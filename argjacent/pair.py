#!/usr/bin/env python3
"""Certified single-sample left-right flank pairing for ARGjacent.

Pairs upstream (`left`) and downstream (`right`) flank trunks across a gene into certified
two-sided `[left + gene + right]` contexts using four strictly ordered evidence tiers:
  1. `span_pair`   : direct paired-end read bridging across short genes (`len(gene) <= 650 bp`)
                     using junction-spanning 21-mers (down to 1-bp overhang past `fork_bp`).
  2. `snv_phase`   : mate-pair haplotype chaining across internal gene SNVs when paired reads
                     cover the gene interior.
  3. `unbranched`  : `1 x 1` topology (after pruning unverified distal ghost forks where
                     `certain_bp <= fork_bp`), plus Kirchhoff flow-partitioned `1 x k` / `k x 1`.
  4. `depth_match` : symmetric log-normal + Poisson delta-method flow bijection (`UNEV = 0.18`),
                     using post-fork exclusive depths from `depth_bins` (when inside `certain_bp`)
                     or boundary depth (`fork_bp < 90`), requiring `z <= 1.8` and permutation
                     `Delta_chi2 >= 9.0` (`>= 3 sigma`) against every alternative assignment.
"""
import collections
import gzip
import itertools
import math
import os
import statistics

UNEV = 0.18       # calibrated single-molecule log-depth standard deviation across ~1-2 kb
DEPTH_BIN = 100   # bin width in verify.py's `depth_bins`
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def pair_z(u, v, unev=UNEV):
    """Symmetric log-ratio z-score between depths u and v under Poisson + UNEV multiplicative drift."""
    u1 = max(float(u), 1.0)
    v1 = max(float(v), 1.0)
    sig_ln = math.sqrt(1.0 / max(u1, 2.0) + 1.0 / max(v1, 2.0) + unev * unev)
    return abs(math.log((u1 + 0.5) / (v1 + 0.5))) / sig_ln


def read_fa_dict(path):
    out, nm, buf = {}, None, []
    with open(path, encoding="latin-1") as f:
        for ln in f:
            ln = ln.strip()
            if ln.startswith(">"):
                if nm is not None:
                    out[nm] = "".join(buf).upper()
                nm, buf = ln[1:].split()[0], []
            elif ln:
                buf.append(ln)
        if nm is not None:
            out[nm] = "".join(buf).upper()
    return out


def prox_seq(flank, side, W):
    if len(flank) <= W:
        return flank
    return flank[-W:] if side == "left" else flank[:W]


def outward_seq(flank, side, W):
    p = prox_seq(flank, side, W)
    return p if side == "right" else p[::-1]


def seq_sim_fast(a, b):
    """Base similarity over shared prefix in outward coordinates, tolerating small indels."""
    L = min(len(a), len(b))
    if L == 0:
        return 1.0
    direct = sum(x == y for x, y in zip(a[:L], b[:L])) / L
    if direct >= 0.97 or L < 40:
        return direct
    best = direct
    for sh in (-12, -8, -4, -2, -1, 1, 2, 4, 8, 12):
        if sh > 0 and L - sh > 20:
            m = sum(x == y for x, y in zip(a[:L - sh], b[sh:L])) / (L - sh)
            if m > best:
                best = m
        elif sh < 0 and L + sh > 20:
            s = -sh
            m = sum(x == y for x, y in zip(a[s:L], b[:L - s])) / (L - s)
            if m > best:
                best = m
    return best


def divergence_bp(out_a, out_b):
    """Outward distance (bp from gene boundary) at which two trunks structurally diverge
    (ignoring isolated 1-bp SNPs by requiring >= 4 mismatches in a 15-bp window)."""
    L = min(len(out_a), len(out_b))
    diffs = [1 if out_a[i] != out_b[i] else 0 for i in range(L)]
    win = sum(diffs[:15])
    if win >= 4:
        for i in range(min(15, L)):
            if diffs[i]:
                return i
    for i in range(15, L):
        win += diffs[i] - diffs[i - 15]
        if win >= 4:
            for j in range(i - 14, i + 1):
                if diffs[j]:
                    return j
    return L


def extract_trunks(paths, side, W=500, sim_thr=0.97):
    """Group candidate paths on `side` into distinct structural proximal trunks up to `W` bp.

    Prunes unverified distal ghost forks (`div >= 60` and `certain_bp <= div`) and computes
    exact post-fork exclusive depths from `depth_bins` whenever the post-fork window lies
    strictly inside `certain_bp`.
    """
    groups = []
    ordered = sorted(
        paths,
        key=lambda p: (-int(p.get("context_len", 0)), -int(p.get("flank_len", 0)), -float(p.get("depth", 0.0)), p.get("id", "")),
    )
    for p in ordered:
        flk = p["flank"]
        if len(flk) < 30:
            continue
        out = outward_seq(flk, side, W)
        placed = False
        for g in groups:
            if seq_sim_fast(out, g["out"]) >= sim_thr:
                g["members"].append(p)
                placed = True
                break
        if not placed:
            groups.append({"out": out, "members": [p]})

    pref = "L" if side == "left" else "R"
    raw_trunks = []
    for g in groups:
        mems = g["members"]
        best = max(mems, key=lambda p: (int(p.get("context_len", 0)), int(p.get("flank_len", 0))))
        cert_mems = [p for p in mems if int(p.get("context_len", 0)) >= W]
        min_sh = min((int(p.get("shared_end", 0)) for p in (cert_mems or mems)), default=0)
        raw_trunks.append({
            "best_cand": best,
            "cand_ids": [p["id"] for p in mems],
            "prox": prox_seq(best["flank"], side, W),
            "out": g["out"],
            "flank_len": int(best.get("flank_len", len(best["flank"]))),
            "certain_bp": max(int(p.get("context_len", 0)) for p in mems),
            "end_bp": max(int(p.get("end_len", p.get("context_len", 0))) for p in mems),
            "shared_end": min_sh,
            "depth": max(float(p.get("depth", 0.0)) for p in mems),
            "boundary_depth": max(float(p.get("depth", 0.0)) for p in mems),
            "depth_bins": best.get("depth_bins"),
            "depth_source": "boundary",
        })

    # Prune ghost forks: unverified candidate branches off an already-certified trunk where verify.py
    # stopped vouching BEFORE or AT the fork point (`div >= 60` and `certain_bp <= div`).
    # Record `pruned_forks` on the surviving parent trunk so 1xk Kirchhoff never treats a partially
    # verified multi-branch side as an unbranched shared arm.
    trunks = []
    for t in raw_trunks:
        t["pruned_forks"] = 0
        is_ghost = False
        if t["certain_bp"] < W:
            for k in trunks:
                if k["certain_bp"] >= W:
                    dv = divergence_bp(t["out"], k["out"])
                    if dv >= 60 and t["certain_bp"] <= dv:
                        k["pruned_forks"] += 1
                        is_ghost = True
                        break
        if not is_ghost:
            t["trunk_id"] = "%s%d" % (pref, len(trunks) + 1)
            trunks.append(t)

    for i, ti in enumerate(trunks):
        mf = W
        for j, tj in enumerate(trunks):
            if i != j:
                mf = min(mf, divergence_bp(ti["out"], tj["out"]))
        ti["min_fork_bp"] = mf
        bins = ti.get("depth_bins")
        # Distal depth near the outer edge of W (used to verify single-trunk flow does not drop before W)
        d_b1 = max(1, W // DEPTH_BIN)
        d_b0 = max(0, d_b1 - 2)
        if bins and d_b0 < len(bins) and d_b1 <= len(bins):
            ti["distal_depth"] = round(float(statistics.median(bins[d_b0:d_b1])), 1)
        else:
            ti["distal_depth"] = ti["depth"]

        if len(trunks) == 1 or mf < 90:
            ti["exclusive_depth_known"] = True
        else:
            # Post-fork exclusive depth from `depth_bins` (`argjacent >= 0.5.0`) when strictly inside `certain_bp`!
            # Check pairwise divergence across ALL other trunks: we need a bin strictly past the LAST fork
            # affecting `ti` within `W`.
            max_f = mf
            for j, tj in enumerate(trunks):
                if i != j:
                    dv = divergence_bp(ti["out"], tj["out"])
                    if dv < W:
                        max_f = max(max_f, dv)
            b0 = (max_f + 120) // DEPTH_BIN
            if bins and b0 < len(bins) and (b0 + 1) * DEPTH_BIN <= ti["certain_bp"]:
                b1 = min(len(bins), b0 + 2, max(b0 + 1, ti["certain_bp"] // DEPTH_BIN))
                sl = bins[b0:b1]
                if sl:
                    ti["depth"] = round(float(statistics.median(sl)), 1)
                    ti["exclusive_depth_known"] = True
                    ti["depth_source"] = "post_fork_bins"
                else:
                    ti["exclusive_depth_known"] = False
            else:
                ti["exclusive_depth_known"] = False
    return trunks


def _iter_read_pairs(reads_source):
    """Yield `(seq1, seq2)` from either `aw.reads.tsv(.gz)` or a `(r1_fq, r2_fq)` tuple."""
    if isinstance(reads_source, (tuple, list)) and len(reads_source) == 2:
        r1, r2 = reads_source
        if not (r1 and r2 and os.path.exists(r1) and os.path.exists(r2)):
            return
        op1 = gzip.open if r1.endswith(".gz") else open
        op2 = gzip.open if r2.endswith(".gz") else open
        with op1(r1, "rt") as f1, op2(r2, "rt") as f2:
            while True:
                h1 = f1.readline()
                if not h1:
                    break
                s1 = f1.readline().rstrip("\n")
                f1.readline(); f1.readline()
                f2.readline()
                s2 = f2.readline().rstrip("\n")
                f2.readline(); f2.readline()
                if s1 and s2:
                    yield s1, s2
    elif isinstance(reads_source, str) and os.path.exists(reads_source):
        op = gzip.open if reads_source.endswith(".gz") else open
        with op(reads_source, "rt") as fh:
            while True:
                l1 = fh.readline()
                if not l1:
                    break
                if l1.startswith("#"):
                    continue
                l2 = fh.readline()
                if not l2:
                    break
                t1 = l1.split("\t")
                t2 = l2.split("\t")
                if len(t1) >= 3 and len(t2) >= 3:
                    yield t1[2].rstrip("\n"), t2[2].rstrip("\n")


def index_sample_reads(genes_dict, doc_genes, reads_source, W=500):
    """Single-pass read-pair indexer across all genes in one sample (`aw.reads.tsv(.gz)` or `(r1_fq, r2_fq)`).
    Returns `{gname: {(li, ri): (tier, note)}}`.
    """
    if not reads_source:
        return {}
    if isinstance(reads_source, str) and not os.path.exists(reads_source):
        return {}
    if isinstance(reads_source, (tuple, list)) and not all(os.path.exists(x) for x in reads_source):
        return {}

    K = 21
    GK = 21
    g_info = {}
    global_seeds = collections.defaultdict(list)

    for gname in sorted(doc_genes):
        gentry = doc_genes[gname]
        gseq = genes_dict.get(gname)
        if not gseq:
            continue
        l_trunks = extract_trunks(gentry.get("left", {}).get("paths", []), "left", W=W)
        r_trunks = extract_trunks(gentry.get("right", {}).get("paths", []), "right", W=W)
        if not l_trunks or not r_trunks:
            continue
        if len(gseq) > 650 and (len(l_trunks) < 2 or len(r_trunks) < 2):
            continue

        l_kmers = collections.defaultdict(set)
        for ti, t in enumerate(l_trunks):
            sq = t["prox"] + gseq[:K - 1]
            for i in range(max(0, len(sq) - W - (K - 1)), len(sq) - K + 1):
                l_kmers[sq[i:i + K]].add(ti)
        r_kmers = collections.defaultdict(set)
        for ti, t in enumerate(r_trunks):
            sq = gseq[-(K - 1):] + t["prox"]
            for i in range(0, min(W + (K - 1), len(sq)) - K + 1):
                r_kmers[sq[i:i + K]].add(ti)

        uniq_L = {km: next(iter(s)) for km, s in l_kmers.items() if len(s) == 1}
        uniq_R = {km: next(iter(s)) for km, s in r_kmers.items() if len(s) == 1}

        g_pos = {}
        for i in range(0, len(gseq) - GK + 1, 5):
            km = gseq[i:i + GK]
            g_pos[km] = -1 if km in g_pos else i
        g_pos = {k: v for k, v in g_pos.items() if v >= 0}

        g_info[gname] = {
            "seq": gseq, "l_trunks": l_trunks, "r_trunks": r_trunks,
            "uniq_L": uniq_L, "uniq_R": uniq_R, "g_pos": g_pos,
            "span_counts": collections.Counter(),
            "gene_pairs": [],
        }
        seed_set = set(uniq_L.keys()) | set(uniq_R.keys())
        for i in range(0, len(gseq) - K + 1, 7):
            seed_set.add(gseq[i:i + K])
        for km in sorted(seed_set):
            global_seeds[km].append(gname)
            global_seeds[rc(km)].append(gname)

    if not g_info:
        return {}

    def scan_read_for_gene(seq, info):
        uniq_L, uniq_R, g_pos, gseq = info["uniq_L"], info["uniq_R"], info["g_pos"], info["seq"]
        for s in (seq, rc(seq)):
            lt = set()
            rt = set()
            for i in range(0, len(s) - K + 1, 3):
                km = s[i:i + K]
                if km in uniq_L:
                    lt.add(uniq_L[km])
                if km in uniq_R:
                    rt.add(uniq_R[km])
            g_off = None
            if len(info["l_trunks"]) >= 2 and len(info["r_trunks"]) >= 2:
                for i in range(0, len(s) - GK + 1, 8):
                    gp = g_pos.get(s[i:i + GK])
                    if gp is not None:
                        g_off = gp - i
                        break
            if lt or rt or g_off is not None:
                snvs = {}
                if g_off is not None:
                    i0 = max(0, -g_off)
                    i1 = min(len(s), len(gseq) - g_off)
                    if i1 - i0 >= 35:
                        mism = 0
                        for i in range(i0, i1):
                            gb = gseq[g_off + i]
                            rb = s[i]
                            if rb != gb and rb in "ACGT":
                                mism += 1
                            snvs[g_off + i] = rb
                        if mism > 0.12 * (i1 - i0):
                            snvs = {}
                return lt, rt, snvs
        return set(), set(), {}

    for s1, s2 in _iter_read_pairs(reads_source):
        cand_genes = set()
        for s in (s1, s2):
            for i in range(0, len(s) - K + 1, 10):
                gs = global_seeds.get(s[i:i + K])
                if gs:
                    cand_genes.update(gs)
        if not cand_genes:
            continue
        for gname in sorted(cand_genes):
            info = g_info[gname]
            lt1, rt1, g1 = scan_read_for_gene(s1, info)
            lt2, rt2, g2 = scan_read_for_gene(s2, info)
            lt = lt1 | lt2
            rt = rt1 | rt2
            if len(lt) == 1 and len(rt) == 1:
                info["span_counts"][(next(iter(lt)), next(iter(rt)))] += 1
            if (lt or rt or g1 or g2) and (g1 or g2):
                merged_g = dict(g1)
                for pos, b in g2.items():
                    if pos in merged_g and merged_g[pos] != b:
                        continue
                    merged_g[pos] = b
                if merged_g:
                    info["gene_pairs"].append((lt, rt, merged_g))

    out_phased = {}
    for gname in sorted(g_info):
        info = g_info[gname]
        phased_tiers = {}
        span_counts = info["span_counts"]
        l_trunks, r_trunks = info["l_trunks"], info["r_trunks"]
        gene_pairs = info["gene_pairs"]
        gseq = info["seq"]

        for (li, ri), cnt in sorted(span_counts.items()):
            other_l = sum(v for (l2, r2), v in span_counts.items() if l2 == li and r2 != ri)
            other_r = sum(v for (l2, r2), v in span_counts.items() if r2 == ri and l2 != li)
            if cnt >= 3 and cnt >= 8 * max(other_l, other_r):
                phased_tiers[(li, ri)] = ("span_pair", "pairs=%d" % cnt)

        if len(l_trunks) >= 2 and len(r_trunks) >= 2 and len(gene_pairs) >= 10:
            pos_bases = collections.defaultdict(collections.Counter)
            for _, _, gm in gene_pairs:
                for pos, b in gm.items():
                    if b in "ACGT":
                        pos_bases[pos][b] += 1
            poly_sites = sorted(
                pos for pos, cnts in pos_bases.items()
                if len([b for b, c in cnts.items() if c >= 4 and c >= 0.10 * sum(cnts.values())]) >= 2
            )
            if poly_sites:
                near5 = [p for p in poly_sites if p <= 450]
                near3 = [p for p in poly_sites if p >= len(gseq) - 450]
                if near5 and near3:
                    for li in range(len(l_trunks)):
                        votes5 = collections.Counter()
                        for lt, _, gm in gene_pairs:
                            if lt == {li} and near5[0] in gm:
                                votes5[gm[near5[0]]] += 1
                        if not votes5:
                            continue
                        b_cur, c_best = sorted(votes5.items(), key=lambda kv: (-kv[1], kv[0]))[0]
                        if c_best < 3 or c_best < 0.85 * sum(votes5.values()):
                            continue
                        chain_sites = [p for p in poly_sites if near5[0] <= p <= near3[-1]]
                        cur_p, cur_b = chain_sites[0], b_cur
                        valid = True
                        for nxt_p in chain_sites[1:]:
                            if nxt_p - cur_p > 420:
                                valid = False
                                break
                            v_step = collections.Counter()
                            for _, _, gm in gene_pairs:
                                if gm.get(cur_p) == cur_b and nxt_p in gm:
                                    v_step[gm[nxt_p]] += 1
                            if not v_step:
                                valid = False
                                break
                            nb, nc = sorted(v_step.items(), key=lambda kv: (-kv[1], kv[0]))[0]
                            if nc < 3 or nc < 0.85 * sum(v_step.values()):
                                valid = False
                                break
                            cur_p, cur_b = nxt_p, nb
                        if not valid:
                            continue
                        v_right = collections.Counter()
                        for _, rt, gm in gene_pairs:
                            if gm.get(cur_p) == cur_b and len(rt) == 1:
                                v_right[next(iter(rt))] += 1
                        if not v_right:
                            continue
                        ri_best, rc_best = sorted(v_right.items(), key=lambda kv: (-kv[1], kv[0]))[0]
                        if rc_best >= 3 and rc_best >= 0.85 * sum(v_right.values()):
                            if (li, ri_best) not in phased_tiers:
                                phased_tiers[(li, ri_best)] = ("snv_phase", "snvs=%d,r_links=%d" % (len(chain_sites), rc_best))
        if phased_tiers:
            out_phased[gname] = phased_tiers
    return out_phased


def _mk_pair(idx, lt, rt, tier, score, note=""):
    return {
        "id": "p%d" % idx,
        "left_trunk": lt["trunk_id"],
        "right_trunk": rt["trunk_id"],
        "left_id": lt["best_cand"]["id"],
        "right_id": rt["best_cand"]["id"],
        "left_certain_bp": lt["certain_bp"],
        "right_certain_bp": rt["certain_bp"],
        "left_depth": lt["depth"],
        "right_depth": rt["depth"],
        "tier": tier,
        "score": round(float(score), 2),
        "note": note,
    }


def pair_gene_single_sample(gene_name, gene_seq, gene_entry, W=500, reads_source=None, pre_phased=None, include_prox=False):
    """Pair left and right proximal trunks for one gene in a single sample."""
    l_paths = gene_entry.get("left", {}).get("paths", [])
    r_paths = gene_entry.get("right", {}).get("paths", [])
    if not l_paths or not r_paths:
        return []

    l_all = extract_trunks(l_paths, "left", W=W)
    r_all = extract_trunks(r_paths, "right", W=W)
    if not l_all or not r_all:
        return []

    def is_cert_trunk(t):
        return t["certain_bp"] >= W and t["shared_end"] == 0

    l_cert = [t for t in l_all if is_cert_trunk(t)]
    r_cert = [t for t in r_all if is_cert_trunk(t)]
    if not l_cert or not r_cert:
        return []

    pairs = []
    used_pairs = set()

    def add_pair(lt, rt, tier, score, note=""):
        key = (lt["trunk_id"], rt["trunk_id"])
        if key in used_pairs:
            return
        used_pairs.add(key)
        pr = _mk_pair(len(pairs) + 1, lt, rt, tier, score, note=note)
        if include_prox:
            pr["left_prox"] = lt["prox"]
            pr["right_prox"] = rt["prox"]
        pairs.append(pr)

    # Tiers 1 & 2: direct paired-end bridging (`span_pair`) & internal SNV phasing (`snv_phase`)
    if pre_phased is not None:
        phased_map = pre_phased
    elif reads_source:
        phased_map = index_sample_reads({gene_name: gene_seq}, {gene_name: gene_entry}, reads_source, W=W).get(gene_name, {})
    else:
        phased_map = {}

    for (li, ri), (tier_name, note) in sorted(phased_map.items()):
        if li < len(l_all) and ri < len(r_all):
            lt, rt = l_all[li], r_all[ri]
            if is_cert_trunk(lt) and is_cert_trunk(rt):
                add_pair(lt, rt, tier_name, 99.0, note=note)

    l_rd = float(gene_entry.get("left", {}).get("read_depth", 0.0))
    r_rd = float(gene_entry.get("right", {}).get("read_depth", 0.0))
    if min(l_rd, r_rd) < 8.0:
        return pairs

    # Tier 3a: `1 x 1` unbranched topology
    if len(l_all) == 1 and len(r_all) == 1 and len(l_cert) == 1 and len(r_cert) == 1:
        lt, rt = l_cert[0], r_cert[0]
        beam_clean = (len(l_paths) < 12 and len(r_paths) < 12)
        u_tol = 0.30 if beam_clean else 0.18
        z_max = 2.2 if beam_clean else 1.8
        if (pair_z(l_rd, r_rd, unev=u_tol) <= z_max
                and pair_z(lt["depth"], rt["depth"], unev=u_tol) <= z_max
                and min(lt["depth"], rt["depth"]) >= 8.0):
            add_pair(lt, rt, "unbranched", 30.0, note="1x1")
        return pairs

    # Whole-gene flow continuity check for multi-trunk depth/Kirchhoff tiers (3b and 4).
    # When all trunks on both sides have known exclusive depths, total flow is sum(depth)
    # (comparing max read_depth would falsely reject a true 1x2 split such as 90x vs 45x+45x).
    excl_l = all(t["exclusive_depth_known"] for t in l_all)
    excl_r = all(t["exclusive_depth_known"] for t in r_all)
    f_l = sum(t["depth"] for t in l_all) if excl_l else max([l_rd] + [t["depth"] for t in l_all])
    f_r = sum(t["depth"] for t in r_all) if excl_r else max([r_rd] + [t["depth"] for t in r_all])
    if pair_z(f_l, f_r, unev=0.18) > 1.8:
        return pairs
    if not (excl_l and excl_r) and pair_z(l_rd, r_rd, unev=0.18) > 1.8:
        return pairs

    # Tier 3b: `1 x k` / `k x 1` Kirchhoff flow partitioning
    if len(l_all) == 1 and len(l_cert) == 1 and len(r_all) >= 2 and len(l_paths) < 12:
        lt = l_cert[0]
        if all(t["exclusive_depth_known"] for t in r_all) and len(r_paths) < 15:
            r_sum = sum(t["depth"] for t in r_all)
            if (len(r_cert) == len(r_all)
                    and lt.get("pruned_forks", 0) == 0
                    and pair_z(lt["depth"], r_sum, unev=0.18) <= 1.6
                    and pair_z(lt.get("distal_depth", lt["depth"]), r_sum, unev=0.18) <= 1.6
                    and all(t["depth"] >= max(8.0, 0.15 * lt["depth"]) for t in r_all)
                    and max(t["depth"] for t in r_all) < 0.88 * lt["depth"]):
                for rt in r_cert:
                    add_pair(lt, rt, "unbranched", 20.0, note="1x%d_kirchhoff" % len(r_all))
                return pairs
            r_sorted = sorted(r_all, key=lambda t: (-t["depth"], t["trunk_id"]))
            r_dom, r_rest = r_sorted[0], sum(t["depth"] for t in r_sorted[1:])
            if is_cert_trunk(r_dom) and pair_z(lt["depth"], r_dom["depth"], unev=0.18) <= 1.5 and r_rest < 0.08 * lt["depth"]:
                add_pair(lt, r_dom, "unbranched", 22.0, note="1xDom")
                return pairs

    if len(r_all) == 1 and len(r_cert) == 1 and len(l_all) >= 2 and len(r_paths) < 12:
        rt = r_cert[0]
        if all(t["exclusive_depth_known"] for t in l_all) and len(l_paths) < 15:
            l_sum = sum(t["depth"] for t in l_all)
            if (len(l_cert) == len(l_all)
                    and rt.get("pruned_forks", 0) == 0
                    and pair_z(rt["depth"], l_sum, unev=0.18) <= 1.6
                    and pair_z(rt.get("distal_depth", rt["depth"]), l_sum, unev=0.18) <= 1.6
                    and all(t["depth"] >= max(8.0, 0.15 * rt["depth"]) for t in l_all)
                    and max(t["depth"] for t in l_all) < 0.88 * rt["depth"]):
                for lt in l_cert:
                    add_pair(lt, rt, "unbranched", 20.0, note="%dx1_kirchhoff" % len(l_all))
                return pairs
            l_sorted = sorted(l_all, key=lambda t: (-t["depth"], t["trunk_id"]))
            l_dom, l_rest = l_sorted[0], sum(t["depth"] for t in l_sorted[1:])
            if is_cert_trunk(l_dom) and pair_z(rt["depth"], l_dom["depth"], unev=0.18) <= 1.5 and l_rest < 0.08 * rt["depth"]:
                add_pair(l_dom, rt, "unbranched", 22.0, note="Domx1")
                return pairs

    # Tier 4: Statistical depth conservation (`depth_match`)
    max_d = max(max(t["depth"] for t in l_all), max(t["depth"] for t in r_all))
    l_maj = [t for t in l_all if t["depth"] >= 0.08 * max_d]
    r_maj = [t for t in r_all if t["depth"] >= 0.08 * max_d]
    if (all(t["exclusive_depth_known"] for t in l_maj + r_maj)
            and len(l_maj) == len(r_maj) and 2 <= len(l_maj) <= 5
            and all(is_cert_trunk(t) for t in l_maj) and all(is_cert_trunk(t) for t in r_maj)):
        P = len(l_maj)
        u = [t["depth"] for t in l_maj]
        v = [t["depth"] for t in r_maj]
        if pair_z(sum(u), sum(v), unev=0.18) <= 1.8:
            perms = []
            for perm in itertools.permutations(range(P)):
                zs = [pair_z(u[a], v[perm[a]], unev=0.18) for a in range(P)]
                perms.append((sum(z * z for z in zs), perm, zs))
            perms.sort(key=lambda x: (round(x[0], 9), x[1]))
            best_cost, best_perm, best_zs = perms[0]
            for a in range(P):
                b = best_perm[a]
                if best_zs[a] > 1.8 or min(u[a], v[b]) < 10.0:
                    continue
                alt_costs = [c for c, pm, _ in perms[1:] if pm[a] != b]
                delta_c = (min(alt_costs) - best_cost) if alt_costs else 0.0
                if delta_c >= 9.0:
                    add_pair(l_maj[a], r_maj[b], "depth_match", delta_c,
                             note="dL=%.0f,dR=%.0f,dChi2=%.1f" % (u[a], v[b], delta_c))
    return pairs


def pair_sample(genes_dict, doc_genes, W=500, reads_source=None, include_prox=False):
    """Run certified single-sample pairing across all genes in one sample (single pass over reads)."""
    phased = index_sample_reads(genes_dict, doc_genes, reads_source, W=W) if reads_source else {}
    out = {}
    for gname in doc_genes:
        gseq = genes_dict.get(gname)
        if not gseq:
            continue
        prs = pair_gene_single_sample(
            gname, gseq, doc_genes[gname], W=W, pre_phased=phased.get(gname, {}), include_prox=include_prox
        )
        if prs:
            out[gname] = prs
    return out


def write_pairs_tsv(pairs_by_gene, out_tsv):
    with open(out_tsv, "w") as f:
        f.write("gene\tpair_id\tleft_context\tright_context\tleft_certain_bp\tright_certain_bp\tleft_depth\tright_depth\ttier\tscore\tnote\n")
        for gname, plist in pairs_by_gene.items():
            for pr in plist:
                f.write("%s\t%s\t%s\t%s\t%d\t%d\t%.1f\t%.1f\t%s\t%.2f\t%s\n" % (
                    gname, pr["id"], pr["left_id"], pr["right_id"],
                    int(pr["left_certain_bp"]), int(pr["right_certain_bp"]),
                    float(pr["left_depth"]), float(pr["right_depth"]),
                    pr["tier"], float(pr["score"]), pr.get("note", ""),
                ))
