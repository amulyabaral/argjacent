#!/usr/bin/env python3
"""Multi-sample cohort co-abundance deconvolution (`argjacent-cohort` / `argjacent cohort`).

Given S >= 2 samples processed by `argjacent`, resolves two classes of ambiguity that cannot
be resolved within a single short-read sample:
  Part A (`cohort_pair_all`):
    Pairs upstream (`left`) and downstream (`right`) flanks across identical long genes when
    single-sample depths tie or share a cassette, via bidirectional weighted Lawson-Hanson NNLS
    (`Y ~ X @ W_fwd`, `X ~ Y @ W_rev`) plus `1 x 1` cohort-anchor transfer and flow-conserved
    `2 x 2` complement subtraction, gated by leave-one-sample-out `Delta_chi2`.
  Part B (`cohort_deconvolve_repeats`):
    Resolves post-repeat exit arms (`[Gene] -> [Shared Repeat >= 400 bp, e.g. IS26] -> [Exits]`)
    by fitting `g_trunk ~ Y_exit @ w` (`w >= 0`) across S >= 4 samples using `gene_depth` and
    `depth_bins` from `flanks.json` directly (and optionally `aw.reads.tsv` when `--keep` was used).
"""
import argparse
import collections
import concurrent.futures
import gzip
import json
import multiprocessing
import os
import statistics
import sys
import time

import numpy as np

from . import __version__ as VERSION
from .gfa import write_gfa
from .pair import (
    DEPTH_BIN,
    divergence_bp,
    extract_trunks,
    pair_gene_single_sample,
    pair_z,
    rc,
    read_fa_dict,
    seq_sim_fast,
)


def _nnls_weighted(A, b, w_sqrt):
    """Exact active-set Lawson-Hanson NNLS: min_{x >= 0} || diag(w_sqrt) (A x - b) ||_2^2.

    For nonnegative A, b, any column j with zero overlap with {i : b[i] > 0} has
    gradient Aw[:, j]^T (bw - Aw x) = -Aw[:, j]^T (Aw x) <= 0 for all x >= 0 and
    can never enter the active set; pruning those columns (and all-zero rows) plus
    using Gram-Cholesky solves with SVD fallback gives identical solutions orders
    of magnitude faster on sparse large-cohort matrices.
    """
    Aw = A * w_sqrt[:, None]
    bw = b * w_sqrt
    n = Aw.shape[1]
    if n == 0:
        return np.zeros(0, dtype=np.float64), float(np.sum(bw ** 2))
    nonneg = bool(np.all(Aw >= 0.0) and np.all(bw >= 0.0))
    b_pos = bw > 0.0
    if nonneg:
        if not np.any(b_pos):
            return np.zeros(n, dtype=np.float64), 0.0
        cols = np.nonzero(np.any(Aw[b_pos] > 0.0, axis=0))[0]
    else:
        cols = np.arange(n)
    if cols.size == 0:
        return np.zeros(n, dtype=np.float64), float(np.sum(bw ** 2))
    Aw_s = Aw[:, cols]
    row_mask = (bw != 0.0) | np.any(Aw_s != 0.0, axis=1)
    Aw_s = Aw_s[row_mask]
    bw_s = bw[row_mask]
    m = cols.size

    G = Aw_s.T @ Aw_s
    c = Aw_s.T @ bw_s
    P_set = np.zeros(m, dtype=bool)
    x_s = np.zeros(m, dtype=np.float64)
    w_grad = c.copy()
    it = 0
    while np.any((~P_set) & (w_grad > 1e-9)) and it < 5 * m + 20:
        it += 1
        j = int(np.argmax(np.where(~P_set, w_grad, -1e30)))
        P_set[j] = True
        while True:
            idx = np.nonzero(P_set)[0]
            if idx.size == 0:
                break
            sp = None
            try:
                L = np.linalg.cholesky(G[np.ix_(idx, idx)])
                diag = np.diag(L)
                if np.min(diag) > 1e-7 * np.max(diag):
                    y = np.linalg.solve(L, c[idx])
                    sp = np.linalg.solve(L.T, y)
            except np.linalg.LinAlgError:
                pass
            if sp is None:
                sp, _, _, _ = np.linalg.lstsq(Aw_s[:, idx], bw_s, rcond=1e-10)
            if np.all(sp > 1e-10):
                x_s[:] = 0.0
                x_s[idx] = sp
                break
            neg = sp <= 1e-10
            alphas = x_s[idx][neg] / (x_s[idx][neg] - sp[neg] + 1e-30)
            alpha = float(np.min(alphas))
            x_new = np.zeros(m, dtype=np.float64)
            x_new[idx] = x_s[idx] + alpha * (sp - x_s[idx])
            x_s = np.maximum(0.0, x_new)
            P_set[idx[x_s[idx] <= 1e-10]] = False
        idx = np.nonzero(P_set)[0]
        w_grad = (c - G[:, idx] @ x_s[idx]) if idx.size else c.copy()
    idx = np.nonzero(P_set)[0]
    res = float(np.sum((bw_s - Aw_s[:, idx] @ x_s[idx]) ** 2)) if idx.size else float(np.sum(bw_s ** 2))
    x = np.zeros(n, dtype=np.float64)
    x[cols] = x_s
    return x, res


class _ExactSeqClusters:
    """Bit-for-bit exact greedy first-match clusterer for `seq_sim_fast(seq, c) >= thr`.

    Stores cluster representatives in a single 0-padded uint8 matrix so each query
    checks all existing clusters via vectorized NumPy comparisons while preserving
    exact first-match insertion order and shift semantics.
    """

    def __init__(self, max_len=500, thr=0.96):
        self.max_len = max(1, int(max_len))
        self.thr = float(thr)
        self.seqs = []
        self.items = []
        self.cap = 64
        self.cnt = 0
        self.mat = np.zeros((self.cap, self.max_len), dtype=np.uint8)
        self.lens = np.zeros(self.cap, dtype=np.int32)

    def _grow_cols(self, new_len):
        nmat = np.zeros((self.cap, new_len), dtype=np.uint8)
        nmat[:, :self.max_len] = self.mat
        self.mat = nmat
        self.max_len = new_len

    def match_or_add(self, seq, payload, add=True):
        L = len(seq)
        if L > self.max_len:
            self._grow_cols(L)
        aa = np.frombuffer(seq.encode("latin-1"), dtype=np.uint8)
        if self.cnt > 0:
            M = self.mat[:self.cnt, :L]
            mL = np.minimum(self.lens[:self.cnt], L)
            zero_hit = np.flatnonzero(mL == 0)
            best = int(zero_hit[0]) if zero_hit.size else None
            direct = np.count_nonzero(M == aa, axis=1) if L > 0 else np.zeros(self.cnt, dtype=np.int64)
            hits = np.flatnonzero((mL > 0) & (direct >= self.thr * mL))
            if hits.size:
                h = int(hits[0])
                if best is None or h < best:
                    best = h
            lim = self.cnt if best is None else best
            if lim > 0 and L >= 40:
                M_lim = M[:lim]
                mL_lim = mL[:lim]
                lens_lim = self.lens[:lim]
                valid40 = mL_lim >= 40
                if np.any(valid40):
                    for sh in (-12, -8, -4, -2, -1, 1, 2, 4, 8, 12):
                        if sh > 0 and L - sh > 20:
                            eff = mL_lim - sh
                            cnts = np.count_nonzero(M_lim[:, sh:L] == aa[:L - sh], axis=1)
                            hs = np.flatnonzero(valid40 & (eff > 20) & (cnts >= self.thr * eff))
                            if hs.size:
                                h = int(hs[0])
                                if best is None or h < best:
                                    best = h
                        elif sh < 0 and L + sh > 20:
                            s = -sh
                            eff = mL_lim - s
                            cnts = np.count_nonzero(M_lim[:, :L - s] == aa[s:L], axis=1)
                            hs = np.flatnonzero(valid40 & (eff > 20) & (cnts >= self.thr * eff))
                            for h_arr in hs:
                                h = int(h_arr)
                                if best is not None and h >= best:
                                    break
                                if lens_lim[h] >= L or seq_sim_fast(seq, self.seqs[h]) >= self.thr:
                                    best = h
                                    break
            if best is not None:
                return best
        if not add:
            return None
        idx = self.cnt
        if idx == self.cap:
            self.cap *= 2
            nmat = np.zeros((self.cap, self.max_len), dtype=np.uint8)
            nmat[:idx] = self.mat
            self.mat = nmat
            nlens = np.zeros(self.cap, dtype=np.int32)
            nlens[:idx] = self.lens
            self.lens = nlens
        if L > 0:
            self.mat[idx, :L] = aa
        self.lens[idx] = L
        self.seqs.append(seq)
        self.items.append(payload)
        self.cnt += 1
        return idx


# ==============================================================================
# Part A: Cohort-wide Left-Right Flank Pairing
# ==============================================================================

def cohort_pair_family(samples_data, fam, W=500, delta_chi2_min=12.0, loo_chi2_min=8.0, include_prox=False):
    """Run cohort-wide left-right pairing for one gene family across S samples."""
    per_sg_single = collections.defaultdict(dict)
    sample_reps = []
    for sd in samples_data:
        sid = sd["sid"]
        doc_genes = sd["doc"].get("genes", {})
        fam_genes = sorted(g for g in doc_genes if g.split("_")[0] == fam and g in sd["genes"])
        if not fam_genes:
            continue
        for g in fam_genes:
            sp = pair_gene_single_sample(
                g, sd["genes"][g], doc_genes[g], W=W,
                pre_phased=sd.get("pre_phased", {}).get(g),
                include_prox=include_prox,
            )
            per_sg_single[sid][g] = list(sp)
        best_g = max(
            fam_genes,
            key=lambda g: (
                max(doc_genes[g].get("left", {}).get("read_depth", 0.0),
                    doc_genes[g].get("right", {}).get("read_depth", 0.0)),
                g,
            ),
        )
        g_ent = doc_genes[best_g]
        l_all = extract_trunks(g_ent.get("left", {}).get("paths", []), "left", W=W)
        r_all = extract_trunks(g_ent.get("right", {}).get("paths", []), "right", W=W)
        if l_all and r_all:
            sample_reps.append((sid, best_g, fam_genes, l_all, r_all, per_sg_single[sid][best_g]))

    if len(sample_reps) < 2:
        return per_sg_single

    U_cl = _ExactSeqClusters(max_len=W, thr=0.96)
    V_cl = _ExactSeqClusters(max_len=W, thr=0.96)

    valid_rows = []
    for sid, best_g, fam_genes, l_all, r_all, sp in sample_reps:
        max_d = max(max(t["depth"] for t in l_all), max(t["depth"] for t in r_all))
        l_maj = [t for t in l_all if t["depth"] >= 0.05 * max_d]
        r_maj = [t for t in r_all if t["depth"] >= 0.05 * max_d]
        l_map = {t["trunk_id"]: U_cl.match_or_add(t["out"], {"out": t["out"], "prox": t["prox"]}) for t in l_all}
        r_map = {t["trunk_id"]: V_cl.match_or_add(t["out"], {"out": t["out"], "prox": t["prox"]}) for t in r_all}
        excl_ok = all(t["exclusive_depth_known"] for t in l_maj + r_maj)
        valid_rows.append((sid, best_g, fam_genes, l_all, r_all, l_map, r_map, excl_ok))

    P, Q = U_cl.cnt, V_cl.cnt
    if P == 0 or Q == 0:
        return per_sg_single

    # Compute co_present and single-sample anchors before NNLS so clusters/families with
    # no co_present >= 2 pairs skip underdetermined NNLS + LOO solves completely.
    anchored = collections.Counter()
    co_present = np.zeros((P, Q), dtype=int)
    for sid, best_g, _, l_all, r_all, l_map, r_map, _ in valid_rows:
        u_here = {l_map[t["trunk_id"]] for t in l_all if t["depth"] >= 8.0 and t["certain_bp"] >= W and t["shared_end"] == 0}
        v_here = {r_map[t["trunk_id"]] for t in r_all if t["depth"] >= 8.0 and t["certain_bp"] >= W and t["shared_end"] == 0}
        for a in u_here:
            for b in v_here:
                co_present[a, b] += 1
        for pr in per_sg_single[sid].get(best_g, []):
            a = l_map.get(pr["left_trunk"])
            b = r_map.get(pr["right_trunk"])
            if a is not None and b is not None:
                anchored[(a, b)] += 1

    need_pair = co_present >= 2
    if not np.any(need_pair):
        return per_sg_single

    reg_rows = [r for r in valid_rows if r[7]]
    Sr = len(reg_rows)
    X = np.zeros((Sr, P), dtype=np.float64)
    Y = np.zeros((Sr, Q), dtype=np.float64)
    for r_i, (_, _, _, l_all, r_all, l_map, r_map, _) in enumerate(reg_rows):
        for t in l_all:
            X[r_i, l_map[t["trunk_id"]]] += t["depth"]
        for t in r_all:
            Y[r_i, r_map[t["trunk_id"]]] += t["depth"]

    W_fwd = np.zeros((P, Q), dtype=np.float64)
    dchi2_fwd = np.zeros((P, Q), dtype=np.float64)
    loo_fwd = np.zeros((P, Q), dtype=np.float64)
    W_rev = np.zeros((Q, P), dtype=np.float64)
    dchi2_rev = np.zeros((P, Q), dtype=np.float64)
    loo_rev = np.zeros((P, Q), dtype=np.float64)

    if Sr >= 2:
        co_excl_mat = (X >= 8.0).T.astype(np.int32) @ (Y >= 8.0).astype(np.int32)
        for b in range(Q):
            if not np.any(need_pair[:, b]):
                continue
            yb = Y[:, b]
            if np.sum(yb > 5.0) < 1:
                continue
            row_scale = np.maximum(yb, np.max(X, axis=1) * 0.15)
            sig = np.sqrt(np.maximum(16.0, yb + (0.22 * row_scale) ** 2))
            w_sqrt = 1.0 / sig
            w_full, res_full = _nnls_weighted(X, yb, w_sqrt)
            W_fwd[:, b] = w_full
            for a in range(P):
                if w_full[a] < 0.25 or not need_pair[a, b]:
                    continue
                X_null = X.copy()
                X_null[:, a] = 0.0
                _, res_null = _nnls_weighted(X_null, yb, w_sqrt)
                dc = max(0.0, res_null - res_full)
                dchi2_fwd[a, b] = dc
                if Sr >= 3 and co_excl_mat[a, b] >= 2 and dc >= loo_chi2_min:
                    best_drop = dc
                    active_s = sorted(
                        (s for s in range(Sr) if X[s, a] > 5.0 or yb[s] > 5.0),
                        key=lambda s: -min(X[s, a], yb[s]),
                    )
                    for drop_s in active_s:
                        keep = [s for s in range(Sr) if s != drop_s]
                        _, rf = _nnls_weighted(X[keep], yb[keep], w_sqrt[keep])
                        _, rn = _nnls_weighted(X_null[keep], yb[keep], w_sqrt[keep])
                        best_drop = min(best_drop, max(0.0, rn - rf))
                        if best_drop < loo_chi2_min:
                            break
                    loo_fwd[a, b] = best_drop

        for a in range(P):
            if not np.any(need_pair[a, :]):
                continue
            xa = X[:, a]
            if np.sum(xa > 5.0) < 1:
                continue
            row_scale = np.maximum(xa, np.max(Y, axis=1) * 0.15)
            sig = np.sqrt(np.maximum(16.0, xa + (0.22 * row_scale) ** 2))
            w_sqrt = 1.0 / sig
            w_full, res_full = _nnls_weighted(Y, xa, w_sqrt)
            W_rev[:, a] = w_full
            for b in range(Q):
                if w_full[b] < 0.25 or not need_pair[a, b]:
                    continue
                Y_null = Y.copy()
                Y_null[:, b] = 0.0
                _, res_null = _nnls_weighted(Y_null, xa, w_sqrt)
                dc = max(0.0, res_null - res_full)
                dchi2_rev[a, b] = dc
                if Sr >= 3 and co_excl_mat[a, b] >= 2 and dc >= loo_chi2_min:
                    best_drop = dc
                    active_s = sorted(
                        (s for s in range(Sr) if Y[s, b] > 5.0 or xa[s] > 5.0),
                        key=lambda s: -min(Y[s, b], xa[s]),
                    )
                    for drop_s in active_s:
                        keep = [s for s in range(Sr) if s != drop_s]
                        _, rf = _nnls_weighted(Y[keep], xa[keep], w_sqrt[keep])
                        _, rn = _nnls_weighted(Y_null[keep], xa[keep], w_sqrt[keep])
                        best_drop = min(best_drop, max(0.0, rn - rf))
                        if best_drop < loo_chi2_min:
                            break
                    loo_rev[a, b] = best_drop

    cohort_links = {}
    for a in range(P):
        for b in range(Q):
            if co_present[a, b] < 2:
                continue
            co_excl = int(np.sum((X[:, a] >= 8.0) & (Y[:, b] >= 8.0))) if Sr >= 2 else 0
            df = dchi2_fwd[a, b] if Sr >= 2 else 0.0
            dr = dchi2_rev[a, b] if Sr >= 2 else 0.0
            lf = loo_fwd[a, b] if Sr >= 3 else 0.0
            lr = loo_rev[a, b] if Sr >= 3 else 0.0
            wf = W_fwd[a, b] if Sr >= 2 else 1.0
            wr = W_rev[b, a] if Sr >= 2 else 1.0
            is_coab = (co_excl >= 2 and max(df, dr) >= delta_chi2_min and min(df, dr) >= 4.0
                       and max(lf, lr) >= loo_chi2_min and (wf >= 0.45 or wr >= 0.45))
            nnls_compat = (co_excl < 2) or (wf >= 0.35 or wr >= 0.35)
            is_anc = (anchored[(a, b)] >= 1 and nnls_compat)
            if is_coab or is_anc:
                src = "cohort_flow" if is_coab else "cohort_anchor"
                sc = max(df, dr) if is_coab else (10.0 * anchored[(a, b)] + max(df, dr))
                cohort_links[(a, b)] = (sc, int(co_present[a, b]), anchored[(a, b)], wf, wr, src)

    if not cohort_links:
        return per_sg_single

    sd_by_id = {x["sid"]: x["doc"].get("genes", {}) for x in samples_data}
    for sid, _, fam_genes, _, _, _, _, _ in valid_rows:
        sd_doc = sd_by_id[sid]
        for g in fam_genes:
            existing = per_sg_single[sid][g]
            used = {(p["left_trunk"], p["right_trunk"]) for p in existing}
            l_paths = sd_doc[g].get("left", {}).get("paths", [])
            r_paths = sd_doc[g].get("right", {}).get("paths", [])
            l_tr = extract_trunks(l_paths, "left", W=W)
            r_tr = extract_trunks(r_paths, "right", W=W)
            if not l_tr or not r_tr:
                continue

            lm = {}
            for t in l_tr:
                ui = U_cl.match_or_add(t["out"], None, add=False)
                if ui is not None:
                    lm[t["trunk_id"]] = ui
            rm = {}
            for t in r_tr:
                vi = V_cl.match_or_add(t["out"], None, add=False)
                if vi is not None:
                    rm[t["trunk_id"]] = vi

            def _append_pair(lt, rt, tier, sc, note):
                pr = {
                    "id": "p%d" % (len(existing) + 1),
                    "left_trunk": lt["trunk_id"],
                    "right_trunk": rt["trunk_id"],
                    "left_id": lt["best_cand"]["id"],
                    "right_id": rt["best_cand"]["id"],
                    "left_certain_bp": lt["certain_bp"],
                    "right_certain_bp": rt["certain_bp"],
                    "left_depth": lt["depth"],
                    "right_depth": rt["depth"],
                    "tier": tier,
                    "score": round(float(sc), 2),
                    "note": note,
                }
                if include_prox:
                    pr["left_prox"] = lt["prox"]
                    pr["right_prox"] = rt["prox"]
                used.add((lt["trunk_id"], rt["trunk_id"]))
                existing.append(pr)

            for lt in l_tr:
                if lt["certain_bp"] < W or lt["shared_end"] != 0 or lt["depth"] < 8.0:
                    continue
                a = lm.get(lt["trunk_id"])
                if a is None:
                    continue
                for rt in r_tr:
                    if rt["certain_bp"] < W or rt["shared_end"] != 0 or rt["depth"] < 8.0:
                        continue
                    b = rm.get(rt["trunk_id"])
                    if b is None or (a, b) not in cohort_links:
                        continue
                    if (lt["trunk_id"], rt["trunk_id"]) in used:
                        continue

                    partner_ls = [t2 for t2 in l_tr if (lm.get(t2["trunk_id"]), b) in cohort_links and t2["depth"] >= 5.0]
                    partner_rs = [t2 for t2 in r_tr if (a, rm.get(t2["trunk_id"])) in cohort_links and t2["depth"] >= 5.0]

                    ok_sample = False
                    if lt["exclusive_depth_known"] and rt["exclusive_depth_known"]:
                        z_direct = pair_z(lt["depth"], rt["depth"])
                        z_many_l = (pair_z(sum(t2["depth"] for t2 in partner_ls), rt["depth"])
                                    if len(partner_ls) > 1 and all(t2["exclusive_depth_known"] for t2 in partner_ls) else 99.0)
                        z_many_r = (pair_z(lt["depth"], sum(t2["depth"] for t2 in partner_rs))
                                    if len(partner_rs) > 1 and all(t2["exclusive_depth_known"] for t2 in partner_rs) else 99.0)
                        if (len(partner_ls) == 1 and len(partner_rs) == 1 and z_direct <= 2.0) or min(z_many_l, z_many_r) <= 1.8:
                            ok_sample = True
                    elif len(partner_ls) == 1 and len(partner_rs) == 1:
                        ok_sample = True
                        if lt["exclusive_depth_known"] and not rt["exclusive_depth_known"] and lt["depth"] > 1.55 * rt["depth"]:
                            ok_sample = False
                        elif rt["exclusive_depth_known"] and not lt["exclusive_depth_known"] and rt["depth"] > 1.55 * lt["depth"]:
                            ok_sample = False

                    if ok_sample:
                        sc, co_n, n_anc, wf, wr, src = cohort_links[(a, b)]
                        _append_pair(lt, rt, src, sc, "co=%d,anc=%d,W=%.2f/%.2f" % (co_n, n_anc, wf, wr))

            l_rd = float(sd_doc[g].get("left", {}).get("read_depth", 0.0))
            r_rd = float(sd_doc[g].get("right", {}).get("read_depth", 0.0))
            if (len(l_tr) == 2 and len(r_tr) == 2 and len(l_paths) < 12 and len(r_paths) < 12
                    and min(l_rd, r_rd) >= 10.0 and pair_z(l_rd, r_rd, unev=0.18) <= 1.6
                    and all(t["certain_bp"] >= W and t["shared_end"] == 0 and t["depth"] >= 8.0 for t in l_tr + r_tr)
                    and len(existing) == 1):
                p0 = existing[0]
                rem_l = [t for t in l_tr if t["trunk_id"] != p0["left_trunk"]]
                rem_r = [t for t in r_tr if t["trunk_id"] != p0["right_trunk"]]
                if len(rem_l) == 1 and len(rem_r) == 1:
                    lt2, rt2 = rem_l[0], rem_r[0]
                    excl2 = lt2["exclusive_depth_known"] and rt2["exclusive_depth_known"]
                    f_l = sum(t["depth"] for t in l_tr) if all(t["exclusive_depth_known"] for t in l_tr) else l_rd
                    f_r = sum(t["depth"] for t in r_tr) if all(t["exclusive_depth_known"] for t in r_tr) else r_rd
                    if pair_z(f_l, f_r, unev=0.18) <= 1.6 and (not excl2 or pair_z(lt2["depth"], rt2["depth"]) <= 1.8):
                        _append_pair(lt2, rt2, "cohort_complement", 9.5,
                                     "2x2_minus_%s_%s" % (p0["left_trunk"], p0["right_trunk"]))
    return per_sg_single


def _pair_family_worker(task):
    fam_samples, fam, W, include_prox = task
    return fam, dict(cohort_pair_family(fam_samples, fam, W=W, include_prox=include_prox))


def cohort_pair_all(samples_data, W=500, include_prox=False, threads=1):
    """Run cohort pairing across all gene families in `samples_data`."""
    fams = sorted({g.split("_")[0] for sd in samples_data for g in sd["doc"].get("genes", {})})
    out = {sd["sid"]: {} for sd in samples_data}
    if threads > 1 and len(fams) > 1:
        tasks = []
        for fam in fams:
            fam_samples = []
            for sd in samples_data:
                doc_genes = sd["doc"].get("genes", {})
                fgenes = {g: doc_genes[g] for g in doc_genes if g.split("_")[0] == fam and g in sd["genes"]}
                if not fgenes:
                    continue
                entry = {
                    "sid": sd["sid"],
                    "genes": {g: sd["genes"][g] for g in fgenes},
                    "doc": {"genes": fgenes},
                }
                if "pre_phased" in sd:
                    entry["pre_phased"] = {g: sd["pre_phased"][g] for g in fgenes if g in sd["pre_phased"]}
                fam_samples.append(entry)
            tasks.append((fam_samples, fam, W, include_prox))
        mp_ctx = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(threads, len(tasks)), mp_context=mp_ctx) as ex:
            results = list(ex.map(_pair_family_worker, tasks))
        for _, res in results:
            for sid, gmap in res.items():
                for g, plist in gmap.items():
                    if plist:
                        out[sid][g] = plist
    else:
        for fam in fams:
            res = cohort_pair_family(samples_data, fam, W=W, include_prox=include_prox)
            for sid, gmap in res.items():
                for g, plist in gmap.items():
                    if plist:
                        out[sid][g] = plist
    return out


# ==============================================================================
# Part B: Cross-Sample Post-Repeat Exit Deconvolution (`[Gene] -> [IS] -> [Exits]`)
# ==============================================================================

def _measure_repeat_from_bins(paths, own_cand, fork_bp, target_bp=5000):
    """Extract `g_depth` and per-exit-cluster `{post_depth, clean_until}` directly from `depth_bins`
    and `gene_depth` written by `argjacent >= 0.5.0` (zero read re-scanning, no `--keep` needed!)."""
    g_vals = [float(p.get("gene_depth", 0.0)) for p in paths if p.get("gene_depth", 0.0) > 0]
    if not g_vals:
        return 0.0, {}
    g_depth = round(max(g_vals), 2)
    b0 = (fork_bp + DEPTH_BIN) // DEPTH_BIN
    max_b = target_bp // DEPTH_BIN

    cstats = {}
    for m_idx, p in own_cand.items():
        bins = p.get("depth_bins")
        if not bins or b0 >= len(bins):
            continue
        post_sl = bins[b0:min(len(bins), b0 + 5)]
        if not post_sl:
            continue
        post_d = float(statistics.median(post_sl))
        clean_until = fork_bp
        for bi in range(b0, min(len(bins), max_b)):
            v = float(bins[bi])
            if v < max(3.0, 0.30 * post_d) or v > 1.85 * max(post_d, 12.0):
                break
            clean_until = min(int(p.get("flank_len", len(p["flank"]))), (bi + 1) * DEPTH_BIN)
        cstats["E%d" % m_idx] = {"post_depth": round(post_d, 2), "clean_until": clean_until}
    return g_depth, cstats


def _measure_repeat_from_reads(gene_seq, synth_paths, side, reads_tsv_path, fork_bp, target_bp=5000):
    """Two-pass k-mer measurement from `aw.reads.tsv(.gz)` when `--keep` was used (recovers even
    minor exit branches dropped by single-sample proposer caps at crowded 4+ way repeat nodes)."""
    if not reads_tsv_path or not os.path.exists(reads_tsv_path) or not synth_paths:
        return 0.0, {}

    K = 31
    RLEN = 150.0
    kmer_counts = collections.Counter()
    tracked = set()

    g_in0 = min(60, max(0, len(gene_seq) // 6))
    g_in1 = max(g_in0 + K + 10, len(gene_seq) - g_in0)
    g_kmers = [gene_seq[i:i + K] for i in range(g_in0, max(g_in0 + 1, g_in1 - K + 1), 4) if i + K <= len(gene_seq)]
    tracked.update(g_kmers)

    post0 = fork_bp + 45
    cand_win_kmers = {}
    for p in synth_paths:
        cid = p["id"]
        flk = p["flank"]
        out = flk if side == "right" else flk[::-1]
        if len(out) < post0 + 120:
            continue
        bins = []
        lim = min(len(out), target_bp)
        for b0 in range(post0, lim - K + 1, 200):
            b1 = min(lim, b0 + 200)
            if b1 - b0 < K + 20:
                break
            b_kms = [out[i:i + K] for i in range(b0, b1 - K + 1, 4)]
            bins.append((b0, b1, b_kms))
            tracked.update(b_kms)
        cand_win_kmers[cid] = bins

    if not tracked:
        return 0.0, {}

    tracked_rc = {rc(x): x for x in tracked}
    op = gzip.open if reads_tsv_path.endswith(".gz") else open
    with op(reads_tsv_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            s = parts[2].rstrip("\n")
            Ls = len(s)
            if Ls < K:
                continue
            hit = False
            for i in range(0, Ls - K + 1, 13):
                sub = s[i:i + K]
                if sub in tracked or sub in tracked_rc:
                    hit = True
                    break
            if not hit:
                continue
            for i in range(0, Ls - K + 1):
                sub = s[i:i + K]
                if sub in tracked:
                    kmer_counts[sub] += 1
                else:
                    fw = tracked_rc.get(sub)
                    if fw is not None:
                        kmer_counts[fw] += 1

    scale = RLEN / (RLEN - K + 1.0)
    g_vals = [kmer_counts[km] * scale for km in g_kmers]
    g_depth = float(np.median(g_vals)) if g_vals else 0.0

    cstats = {}
    for cid, bins in cand_win_kmers.items():
        if not bins:
            continue
        diag = bins[0][2]
        if not diag:
            continue
        post_d = float(np.median([kmer_counts[km] * scale for km in diag]))
        clean_until = fork_bp
        for b0, b1, b_kms in bins:
            vals = [kmer_counts[km] * scale for km in b_kms]
            med_v = float(np.median(vals))
            min_v = float(np.percentile(vals, 10))
            if min_v < 2.0 or med_v < max(4.0, 0.30 * post_d) or med_v > 1.85 * max(post_d, 12.0):
                break
            clean_until = b1
        cstats[cid] = {"post_depth": round(post_d, 2), "clean_until": clean_until}
    return round(g_depth, 2), cstats


def cohort_deconvolve_repeats(samples_data, gene_name, side, min_repeat_bp=400, target_bp=5000,
                              delta_chi2_min=12.0, loo_chi2_min=6.0):
    """Resolve post-repeat exit arms (`[Gene] -> [Shared Repeat >= min_repeat_bp] -> [Exits]`)
    across S >= 4 samples via weighted NNLS (`g_trunk ~ Y_exit @ w`, `w >= 0`).

    Works directly from `depth_bins` + `gene_depth` in `argjacent >= 0.5.0` `flanks.json` (no `--keep`
    required), and automatically uses `reads_tsv` when available to recover proposer-capped minor arms.
    """
    active_samples = []
    exit_clusters = []
    fork_votes = []

    for sd in samples_data:
        gseq = sd["genes"].get(gene_name)
        gent = sd["doc"].get("genes", {}).get(gene_name)
        if not gseq or not gent:
            continue
        paths = gent.get(side, {}).get("paths", [])
        outs = [(p, p["flank"] if side == "right" else p["flank"][::-1])
                for p in paths if len(p["flank"]) >= min_repeat_bp + 250]
        if not outs:
            continue
        has_bins = any(p.get("depth_bins") and p.get("gene_depth", 0.0) > 0 for p, _ in outs)
        has_reads = bool(sd.get("reads_tsv") and os.path.exists(sd["reads_tsv"]))
        if not (has_bins or has_reads):
            continue
        for i in range(len(outs)):
            for j in range(i + 1, len(outs)):
                dv = divergence_bp(outs[i][1], outs[j][1])
                if min_repeat_bp <= dv < min(len(outs[i][1]), len(outs[j][1])) - 50:
                    fork_votes.append(dv)
        active_samples.append((sd, gseq, paths, outs))

    if len(active_samples) < 4 or not fork_votes:
        return {}
    fork_bp = int(min(fork_votes))
    ex_cl = _ExactSeqClusters(max_len=400, thr=0.96)

    for sd, gseq, paths, outs in active_samples:
        for p, out in outs:
            if len(out) < fork_bp + 300:
                continue
            pseq = out[fork_bp:fork_bp + 400]
            prev_cnt = ex_cl.cnt
            ex_cl.match_or_add(pseq, None, add=True)
            if ex_cl.cnt > prev_cnt:
                exit_clusters.append({
                    "id": "E%d" % len(exit_clusters),
                    "flank": p["flank"],
                    "out": out,
                    "post_seq": pseq,
                })

    M = len(exit_clusters)
    if M == 0:
        return {}

    synth_paths = [{"id": cl["id"], "flank": cl["flank"]} for cl in exit_clusters]
    sample_obs = []
    all_have_reads = all(bool(sd.get("reads_tsv") and os.path.exists(sd["reads_tsv"])) for sd, _, _, _ in active_samples)

    for sd, gseq, paths, outs in active_samples:
        own_cand = {}
        for p, out in outs:
            if len(out) >= fork_bp + 250:
                pseq = out[fork_bp:fork_bp + 400]
                m_idx = ex_cl.match_or_add(pseq, None, add=False)
                if m_idx is not None:
                    if m_idx not in own_cand:
                        own_cand[m_idx] = p
                    else:
                        for m2 in range(m_idx + 1, len(exit_clusters)):
                            if m2 not in own_cand and seq_sim_fast(pseq, exit_clusters[m2]["post_seq"]) >= 0.96:
                                own_cand[m2] = p
                                break
        if all_have_reads:
            g_depth, cstats = _measure_repeat_from_reads(gseq, synth_paths, side, sd["reads_tsv"], fork_bp, target_bp=target_bp)
            source_mode = "reads_tsv"
        else:
            g_depth, cstats = _measure_repeat_from_bins(paths, own_cand, fork_bp, target_bp=target_bp)
            source_mode = "depth_bins"
        if g_depth < 5.0 or not cstats:
            continue
        sample_obs.append({
            "sid": sd["sid"], "fork_bp": fork_bp, "g_depth": g_depth,
            "paths": paths, "cstats": cstats, "own_cand": own_cand, "source": source_mode,
        })

    S = len(sample_obs)
    if S < 4:
        return {}

    g_vec = np.array([row["g_depth"] for row in sample_obs], dtype=np.float64)
    Y_mat = np.zeros((S, M), dtype=np.float64)
    for s_i, row in enumerate(sample_obs):
        for m_idx, cl in enumerate(exit_clusters):
            st = row["cstats"].get(cl["id"])
            if st is not None:
                Y_mat[s_i, m_idx] = st["post_depth"]

    sig = np.sqrt(np.maximum(16.0, g_vec + (0.18 * g_vec) ** 2))
    w_sqrt = 1.0 / sig
    w_full, res_full = _nnls_weighted(Y_mat, g_vec, w_sqrt)
    k_pos = int(np.sum(w_full > 0.15))

    # Completeness guard:
    # With `reads_tsv`, all M cohort exits are measured in every sample (`chi2_red <= 2.2`).
    # In JSON-only `depth_bins` mode, `awx` may omit a minor exit in a sample where decoys dominate
    # (`len(own_cand) < M`). Allow excluding up to 1 (or 2 when S >= 6) genuine proposer-dropout
    # samples if the remaining `>= 4` samples fit cleanly with `chi2_red <= 2.2` (never excluding
    # samples that have all M exits measured, and never certifying an excluded sample).
    fit_idx = list(range(S))
    dof = max(1, S - k_pos)
    chi2_red = res_full / dof
    if chi2_red > 2.2:
        if all_have_reads or S < 5:
            return {}
        drop_cands = [i for i in range(S) if len(sample_obs[i]["own_cand"]) < M]
        if not drop_cands:
            return {}
        best_sub = None
        # Try dropping 1 proposer-dropout sample first
        for di in drop_cands:
            sub = [i for i in range(S) if i != di]
            w_s, r_s = _nnls_weighted(Y_mat[sub], g_vec[sub], w_sqrt[sub])
            kp_s = int(np.sum(w_s > 0.15))
            c_s = r_s / max(1, len(sub) - kp_s)
            if c_s <= 2.2 and (best_sub is None or c_s < best_sub[0]):
                best_sub = (c_s, sub, w_s, r_s)
        # If S >= 6 and 2 samples had proposer dropouts, try keeping only complete samples (>= 4)
        if best_sub is None and S >= 6 and len(drop_cands) == 2:
            sub = [i for i in range(S) if i not in drop_cands]
            if len(sub) >= 4:
                w_s, r_s = _nnls_weighted(Y_mat[sub], g_vec[sub], w_sqrt[sub])
                kp_s = int(np.sum(w_s > 0.15))
                c_s = r_s / max(1, len(sub) - kp_s)
                if c_s <= 2.2:
                    best_sub = (c_s, sub, w_s, r_s)
        if best_sub is None:
            return {}
        _, fit_idx, w_full, res_full = best_sub

    Y_fit = Y_mat[fit_idx]
    g_fit = g_vec[fit_idx]
    w_fit = w_sqrt[fit_idx]
    Sf = len(fit_idx)

    validated_exits = {}
    for m in range(M):
        wm = float(w_full[m])
        n_act = int(np.sum(Y_fit[:, m] >= 8.0))
        if n_act < 3 or not (0.65 <= wm <= 1.40):
            continue
        Y_null = Y_fit.copy()
        Y_null[:, m] = 0.0
        _, res_null = _nnls_weighted(Y_null, g_fit, w_fit)
        dchi2 = max(0.0, res_null - res_full)
        if dchi2 < delta_chi2_min:
            continue
        loo_min = dchi2
        drop_order = sorted(range(Sf), key=lambda s: -Y_fit[s, m])
        for drop_s in drop_order:
            keep = [s for s in range(Sf) if s != drop_s]
            _, rf = _nnls_weighted(Y_fit[keep], g_fit[keep], w_fit[keep])
            _, rn = _nnls_weighted(Y_null[keep], g_fit[keep], w_fit[keep])
            loo_min = min(loo_min, max(0.0, rn - rf))
            if loo_min < loo_chi2_min:
                break
        if loo_min >= loo_chi2_min:
            validated_exits[m] = {"w": round(wm, 3), "dChi2": round(dchi2, 2), "loo_dChi2": round(loo_min, 2), "n_samples": n_act}

    if not validated_exits:
        return {}

    resolved = {}
    for s_i, row in enumerate(sample_obs):
        sid = row["sid"]
        active_sum = sum(Y_mat[s_i, m2] * validated_exits[m2]["w"] for m2 in validated_exits if Y_mat[s_i, m2] >= 5.0)
        if pair_z(row["g_depth"], active_sum, unev=0.18) > 2.0:
            continue
        exts = []
        for m, vinfo in sorted(validated_exits.items()):
            cl = exit_clusters[m]
            st = row["cstats"].get(cl["id"])
            if not st or st["post_depth"] < 8.0 or st["clean_until"] <= fork_bp + 100:
                continue
            own_p = row["own_cand"].get(m)
            # Without `reads_tsv`, only certify exits actually proposed in this sample!
            if own_p is None and row["source"] != "reads_tsv":
                continue
            base_cert = int(own_p["context_len"]) if own_p else 0
            flk_seq = own_p["flank"] if own_p else cl["flank"]
            exts.append({
                "cand_id": own_p["id"] if own_p else ("cohort_%s" % cl["id"]),
                "single_certain_bp": base_cert,
                "fork_bp": fork_bp,
                "cohort_certified_bp": max(base_cert, st["clean_until"] - 30),
                "gene_depth": row["g_depth"],
                "exit_post_depth": st["post_depth"],
                "source": row["source"],
                "flank": flk_seq,
                **vinfo,
            })
        if exts:
            resolved[sid] = exts
    return resolved


def _deconvolve_repeat_worker(task):
    gene_samples, gname, side, min_repeat_bp, target_bp = task
    return gname, side, cohort_deconvolve_repeats(
        gene_samples, gname, side, min_repeat_bp=min_repeat_bp, target_bp=target_bp
    )


def cohort_deconvolve_all_repeats(samples_data, min_repeat_bp=400, target_bp=5000, threads=1):
    """Run post-repeat deconvolution across every `(gene, side)` present in `>= 4` samples."""
    gene_counts = collections.Counter(g for sd in samples_data for g in sd["doc"].get("genes", {}))
    out = collections.defaultdict(lambda: collections.defaultdict(dict))
    cand_genes = sorted(g for g, c in gene_counts.items() if c >= 4)
    if threads > 1 and cand_genes:
        tasks = []
        for gname in cand_genes:
            for side in ("left", "right"):
                gene_samples = []
                for sd in samples_data:
                    gseq = sd["genes"].get(gname)
                    gent = sd["doc"].get("genes", {}).get(gname)
                    if not gseq or not gent:
                        continue
                    gene_samples.append({
                        "sid": sd["sid"],
                        "genes": {gname: gseq},
                        "doc": {"genes": {gname: {side: gent.get(side, {})}}},
                        "reads_tsv": sd.get("reads_tsv"),
                    })
                tasks.append((gene_samples, gname, side, min_repeat_bp, target_bp))
        mp_ctx = multiprocessing.get_context("spawn")
        with concurrent.futures.ProcessPoolExecutor(max_workers=min(threads, len(tasks)), mp_context=mp_ctx) as ex:
            results = list(ex.map(_deconvolve_repeat_worker, tasks))
        for gname, side, res in results:
            for sid, exts in res.items():
                out[sid][gname][side] = exts
    else:
        for gname in cand_genes:
            for side in ("left", "right"):
                res = cohort_deconvolve_repeats(samples_data, gname, side, min_repeat_bp=min_repeat_bp, target_bp=target_bp)
                for sid, exts in res.items():
                    out[sid][gname][side] = exts
    return {sid: dict(gmap) for sid, gmap in out.items()}


def _gfa_export_worker(task):
    genes, doc_genes, out_path, pairs_by_gene, cohort_exts_by_gene = task
    write_gfa(
        genes, doc_genes, out_path,
        pairs_by_gene=pairs_by_gene,
        cohort_exts_by_gene=cohort_exts_by_gene,
        verify=True,
    )
    return out_path


# ==============================================================================
# CLI entry point: `argjacent-cohort` / `argjacent cohort`
# ==============================================================================

def _guess_reads_tsv(flanks_json_path):
    """If `PREFIX.tmp/aw.reads.tsv(.gz)` exists next to `PREFIX.flanks.json`, return it."""
    if flanks_json_path.endswith(".flanks.json"):
        pre = flanks_json_path[:-len(".flanks.json")]
        for cand in (pre + ".tmp/aw.reads.tsv", pre + ".tmp/aw.reads.tsv.gz"):
            if os.path.isfile(cand):
                return cand
    return None


def _load_cohort_inputs(args):
    samples = []
    shared_genes = read_fa_dict(args.genes) if args.genes else None
    if args.sample_sheet:
        with open(args.sample_sheet) as f:
            header = f.readline().rstrip("\r\n").split("\t")
            col = {h.strip(): i for i, h in enumerate(header)}
            if "sample" not in col or "flanks_json" not in col:
                sys.stderr.write("argjacent-cohort: --sample-sheet needs at least 'sample' and 'flanks_json' columns\n")
                sys.exit(1)
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                sid = parts[col["sample"]].strip()
                jp = parts[col["flanks_json"]].strip()
                gf = parts[col["genes_fa"]].strip() if "genes_fa" in col and len(parts) > col["genes_fa"] and parts[col["genes_fa"]].strip() else None
                rt = parts[col["reads_tsv"]].strip() if "reads_tsv" in col and len(parts) > col["reads_tsv"] and parts[col["reads_tsv"]].strip() else _guess_reads_tsv(jp)
                gdict = read_fa_dict(gf) if gf else shared_genes
                if not gdict:
                    sys.stderr.write("argjacent-cohort: supply -g genes.fa or a 'genes_fa' column in --sample-sheet\n")
                    sys.exit(1)
                samples.append({"sid": sid, "genes": gdict, "doc": json.load(open(jp)), "reads_tsv": rt, "json_path": jp})
    else:
        if not args.inputs or not shared_genes:
            sys.stderr.write("argjacent-cohort: provide -g genes.fa and -i SAMPLE1.flanks.json SAMPLE2.flanks.json ... (or --sample-sheet)\n")
            sys.exit(1)
        seen_ids = set()
        for idx, item in enumerate(args.inputs, 1):
            if "=" in item and not os.path.exists(item):
                sid, jp = item.split("=", 1)
            else:
                jp = item
                base = os.path.basename(jp)
                sid = base[:-len(".flanks.json")] if base.endswith(".flanks.json") else os.path.splitext(base)[0]
            if sid in seen_ids:
                sid = "%s_%d" % (sid, idx)
            seen_ids.add(sid)
            samples.append({"sid": sid, "genes": shared_genes, "doc": json.load(open(jp)),
                            "reads_tsv": _guess_reads_tsv(jp), "json_path": jp})
    return samples


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="argjacent-cohort",
        description="Multi-sample cohort co-abundance deconvolution for ARGjacent outputs: "
                    "(A) two-sided left-right flank pairing across ambiguous samples, and "
                    "(B) post-repeat exit resolution across >400 bp shared repeats (e.g. IS26).",
    )
    ap.add_argument("-g", "--genes", help="shared genes.fa FASTA (or provide per-sample 'genes_fa' in --sample-sheet)")
    ap.add_argument("-i", "--inputs", nargs="+", metavar="FLANKS_JSON",
                    help="two or more PREFIX.flanks.json files (optionally named as SAMPLE_ID=path/to/run.flanks.json)")
    ap.add_argument("--sample-sheet", metavar="TSV",
                    help="tab-separated table with columns: sample, flanks_json [, genes_fa, reads_tsv]")
    ap.add_argument("-o", "--out", required=True, help="output prefix (writes PREFIX.cohort_pairs.tsv, PREFIX.cohort_repeats.tsv, PREFIX.cohort.json)")
    ap.add_argument("--pair-bp", type=int, default=500, help="certified flank radius W for two-sided pairing (default 500)")
    ap.add_argument("--min-repeat-bp", type=int, default=400, help="minimum shared repeat length for post-repeat deconvolution (default 400)")
    ap.add_argument("--target-bp", type=int, default=5000, help="maximum outward distance for post-repeat cohort certification (default 5000)")
    ap.add_argument("--gfa-dir", metavar="DIR", help="optional directory to write updated per-sample GFA v1.0 graphs with cohort paths")
    ap.add_argument("-t", "--threads", type=int, default=1, help="worker processes for cohort pairing, repeat deconvolution, and GFA export (default 1)")
    ap.add_argument("-v", "--version", action="version", version="argjacent-cohort (argjacent %s)" % VERSION)
    args = ap.parse_args(argv)
    threads = max(1, int(args.threads))

    t0 = time.time()
    samples = _load_cohort_inputs(args)
    if len(samples) < 2:
        sys.stderr.write("argjacent-cohort: at least 2 samples are required (got %d)\n" % len(samples))
        return 1
    sys.stderr.write("argjacent-cohort: loaded %d samples in %.1fs (threads=%d)\n" % (len(samples), time.time() - t0, threads))
    sys.stderr.flush()

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    t_a = time.time()
    pairs_by_sid = cohort_pair_all(samples, W=args.pair_bp, threads=threads)
    sys.stderr.write("argjacent-cohort: Part A (two-sided pairing) finished in %.1fs\n" % (time.time() - t_a))
    sys.stderr.flush()

    t_b = time.time()
    repeats_by_sid = cohort_deconvolve_all_repeats(
        samples, min_repeat_bp=args.min_repeat_bp, target_bp=args.target_bp, threads=threads
    )
    sys.stderr.write("argjacent-cohort: Part B (post-repeat exits) finished in %.1fs\n" % (time.time() - t_b))
    sys.stderr.flush()

    n_pairs = n_cohort_pairs = 0
    with open(args.out + ".cohort_pairs.tsv", "w") as f:
        f.write("sample\tgene\tpair_id\tleft_context\tright_context\tleft_certain_bp\tright_certain_bp\tleft_depth\tright_depth\ttier\tscore\tnote\n")
        for sd in samples:
            sid = sd["sid"]
            for gname, plist in sorted(pairs_by_sid.get(sid, {}).items()):
                for pr in plist:
                    n_pairs += 1
                    if pr["tier"].startswith("cohort_"):
                        n_cohort_pairs += 1
                    f.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.1f\t%.1f\t%s\t%.2f\t%s\n" % (
                        sid, gname, pr["id"], pr["left_id"], pr["right_id"],
                        int(pr["left_certain_bp"]), int(pr["right_certain_bp"]),
                        float(pr["left_depth"]), float(pr["right_depth"]),
                        pr["tier"], float(pr["score"]), pr.get("note", ""),
                    ))

    n_rep = 0
    with open(args.out + ".cohort_repeats.tsv", "w") as f:
        f.write("sample\tgene\tside\tcontext\tfork_bp\tsingle_certain_bp\tcohort_certified_bp\tgene_depth\texit_post_depth\tnnls_weight\tdelta_chi2\tloo_delta_chi2\tn_samples\tsource\n")
        for sd in samples:
            sid = sd["sid"]
            for gname, smap in sorted(repeats_by_sid.get(sid, {}).items()):
                for side, exts in sorted(smap.items()):
                    for ex in exts:
                        n_rep += 1
                        f.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%.1f\t%.1f\t%.3f\t%.2f\t%.2f\t%d\t%s\n" % (
                            sid, gname, side, ex["cand_id"], int(ex["fork_bp"]),
                            int(ex["single_certain_bp"]), int(ex["cohort_certified_bp"]),
                            float(ex["gene_depth"]), float(ex["exit_post_depth"]),
                            float(ex["w"]), float(ex["dChi2"]), float(ex["loo_dChi2"]),
                            int(ex["n_samples"]), ex["source"],
                        ))

    # Strip full flank strings from compact cohort JSON summary? Keep `flank` truncated only in TSV, full in JSON!
    with open(args.out + ".cohort.json", "w") as f:
        json.dump({
            "version": VERSION,
            "parameters": {"pair_bp": args.pair_bp, "min_repeat_bp": args.min_repeat_bp, "target_bp": args.target_bp},
            "samples": [sd["sid"] for sd in samples],
            "pairs": pairs_by_sid,
            "repeats": repeats_by_sid,
        }, f, indent=1)

    if args.gfa_dir:
        t_g = time.time()
        os.makedirs(args.gfa_dir, exist_ok=True)
        if threads > 1 and len(samples) > 1:
            gfa_tasks = [
                (sd["genes"], sd["doc"].get("genes", {}),
                 os.path.join(args.gfa_dir, sd["sid"] + ".gfa"),
                 pairs_by_sid.get(sd["sid"], {}),
                 repeats_by_sid.get(sd["sid"], {}))
                for sd in samples
            ]
            mp_ctx = multiprocessing.get_context("spawn")
            with concurrent.futures.ProcessPoolExecutor(max_workers=min(threads, len(gfa_tasks)), mp_context=mp_ctx) as ex:
                list(ex.map(_gfa_export_worker, gfa_tasks))
        else:
            for sd in samples:
                sid = sd["sid"]
                write_gfa(
                    sd["genes"], sd["doc"].get("genes", {}),
                    os.path.join(args.gfa_dir, sid + ".gfa"),
                    pairs_by_gene=pairs_by_sid.get(sid, {}),
                    cohort_exts_by_gene=repeats_by_sid.get(sid, {}),
                    verify=True,
                )
        sys.stderr.write("argjacent-cohort: wrote %d GFA files in %.1fs\n" % (len(samples), time.time() - t_g))
        sys.stderr.flush()

    print("argjacent-cohort %s: %d samples | %d certified two-sided pairs (%d newly resolved by cohort) | %d post-repeat exits certified" % (
        VERSION, len(samples), n_pairs, n_cohort_pairs, n_rep))
    print("Results: %s.cohort_pairs.tsv  %s.cohort_repeats.tsv  %s.cohort.json" % (args.out, args.out, args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
