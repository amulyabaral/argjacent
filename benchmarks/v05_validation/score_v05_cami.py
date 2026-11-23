#!/usr/bin/env python3
"""Unified genome-truth validation scorer for ARGjacent v0.5.0 on CAMI II / III datasets.

Evaluates all 4 outputs of frozen `argjacent v0.5.0` against exact source genomes:
  1. One-sided certified flanks (`context_len` and `end_len`) at >=100 bp, >=1 kb, >=2 kb, >=5 kb:
       - Strict alignment verification matching `claims_lr3.py` v3.1 (<=3% diff per 3-kb window,
         <=25 diff bases per 300-bp window, no indel >=50 bp from -100 bp inside the gene to
         certain_bp) PLUS molecule-exact (<=6 diffs/300 bp) and lenient checks.
       - Automatic triage of any strict failure (gene-end divergence, draft contig-end truncation,
         or BLAST against all present genomes `cov > 0` for unseeded homologous loci).
       - Independent 31-mer sequence groups and exact 95% Clopper-Pearson confidence intervals.
       - Structural recall (`score_contexts.py` style) for all true contexts with `cov >= min_cov`,
         reported both raw and with circular plasmid rotations collapsed (for CAMI III).
  2. Lossless GFA v1.0 compaction ratio and independent 100% byte-for-byte path round-trip check.
  3. Single-sample two-sided `[left + gene + right]` pairing (`PREFIX.pairs.tsv` at W=500 bp and
     W=1,000 bp) vs Cartesian, Main-to-Main, and Unguarded Greedy baselines, broken down by
     evidence tier and `v0.5.0` `post_fork_bins` usage.
"""
import argparse
import bisect
import collections
import json
import math
import os
import re
import statistics
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "benchmarks", "cami"))

import score_cami as sc
from score_contexts import kindex, struct_len
from argjacent.pair import (
    extract_trunks,
    index_sample_reads,
    pair_gene_single_sample,
    pair_z,
    seq_sim_fast,
)

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
PADG, EDGE, GENE_IN, SUBW, K31 = 300, 30, 100, 300, 31
STRICT_SUB = 25
MOL_SUB = 6
TOL_STRUCT = 40


def cp_ci(k, n, a=0.05):
    if n == 0:
        return (0.0, 1.0)
    def lbin(i, n_, p):
        if p <= 0:
            return 0.0 if i == 0 else -1e300
        if p >= 1:
            return 0.0 if i == n_ else -1e300
        return math.lgamma(n_ + 1) - math.lgamma(i + 1) - math.lgamma(n_ - i + 1) + i * math.log(p) + (n_ - i) * math.log1p(-p)
    def cdf(k_, n_, p):
        return sum(math.exp(lbin(i, n_, p)) for i in range(k_ + 1))
    lo, hi = 0.0, 1.0
    if k > 0:
        x, y = 0.0, 1.0
        for _ in range(55):
            m = (x + y) / 2
            if 1 - cdf(k - 1, n, m) < a / 2:
                x = m
            else:
                y = m
        lo = x
    if k < n:
        x, y = 0.0, 1.0
        for _ in range(55):
            m = (x + y) / 2
            if cdf(k, n, m) > a / 2:
                x = m
            else:
                y = m
        hi = x
    return lo, hi


def family(g):
    return g.split("_")[0]


class Aln:
    __slots__ = ("o0", "o1", "qlen", "qs", "qe", "pos", "w", "cum", "longs", "qid", "hot")

    def __init__(self, o0, o1, qlen, qs, qe, diffs, longs, qid):
        self.o0, self.o1, self.qlen, self.qs, self.qe, self.qid = o0, o1, qlen, qs, qe, qid
        diffs.sort()
        self.pos = [d[0] for d in diffs]
        self.w = [d[1] for d in diffs]
        c, s = [0], 0
        for x in self.w:
            s += x
            c.append(s)
        self.cum = c
        self.longs = sorted(longs)
        self.hot = {}

    def bases(self, x, y):
        i, j = bisect.bisect_left(self.pos, x), bisect.bisect_left(self.pos, y)
        return self.cum[j] - self.cum[i]

    def hot_list(self, cap):
        h = self.hot.get(cap)
        if h is None:
            fr, la = [], []
            n, j = len(self.pos), 0
            for i in range(n):
                if j < i:
                    j = i
                while j + 1 < n and self.pos[j + 1] < self.pos[i] + SUBW:
                    j += 1
                if self.cum[j + 1] - self.cum[i] > cap:
                    fr.append(self.pos[i])
                    la.append(self.pos[j])
            suf = [0] * len(la)
            m = 10 ** 12
            for k in range(len(la) - 1, -1, -1):
                m = min(m, la[k])
                suf[k] = m
            h = self.hot[cap] = (fr, suf)
        return h

    def hot_inside(self, cap, x, y):
        if y - x < SUBW:
            return self.bases(x, y) > cap
        fr, suf = self.hot_list(cap)
        k = bisect.bisect_left(fr, x)
        return k < len(fr) and suf[k] < y

    def long_inside(self, thr, x, y):
        return any(n >= thr and x <= p < y for p, n in self.longs)


def holds_interval(z, x, y, thr, div, sub, fl):
    if not (z.o0 <= x or (x <= -GENE_IN + EDGE and z.o0 <= -GENE_IN + EDGE)):
        return False
    if not (z.o1 >= y or (y >= fl - EDGE and z.o1 >= fl - EDGE)):
        return False
    if z.long_inside(thr, x, y):
        return False
    if z.bases(x, y) > div * (y - x):
        return False
    if sub is not None and z.hot_inside(sub, x, y):
        return False
    return True


def check_claim_strict(alns, claim_len, start_x=-GENE_IN, thr=50, div=0.03, sub=STRICT_SUB):
    """Exact claims_lr3.py v3.1 window check over [start_x, claim_len]."""
    W = 3000
    wins = [(start_x, e) for e in range(start_x + 400, min(claim_len, start_x + W), 100)]
    wins.append((start_x, min(claim_len, start_x + W)))
    wins += [(x, x + W) for x in range(start_x + 100, claim_len - W + 1, 100)]
    if wins[-1][1] < claim_len:
        wins.append((max(start_x, claim_len - W), claim_len))
    for z in alns:
        if all(holds_interval(z, wx, wy, thr, div, sub, claim_len) for wx, wy in wins):
            return True, z.qid
    return False, None


def seq_ident_kmer(a, b, k=15):
    if len(a) < k or len(b) < k:
        return 0.0
    sa = {a[i:i + k] for i in range(len(a) - k + 1)}
    sb = {b[i:i + k] for i in range(len(b) - k + 1)}
    return len(sa & sb) / max(1, min(len(sa), len(sb)))


def match_flank_w(pred_5to3, true_5to3, side, W, strict=False):
    p = pred_5to3[-W:] if side == "left" else pred_5to3[:W]
    t = true_5to3[-W:] if side == "left" else true_5to3[:W]
    if len(p) < W or len(t) < W:
        return False
    if strict:
        # Outward comparison at >= 97% base similarity
        po = p[::-1] if side == "left" else p
        to = t[::-1] if side == "left" else t
        return seq_sim_fast(po, to) >= 0.97
    if seq_ident_kmer(p, t, k=15) >= 0.78:
        return True
    po = p[::-1] if side == "left" else p
    to = t[::-1] if side == "left" else t
    return seq_sim_fast(po, to) >= 0.94


def verify_gfa_file(gfa_path, genes_dict, doc_genes):
    """Independently parse PREFIX.gfa and verify 100% byte-for-byte path reconstruction."""
    if not os.path.isfile(gfa_path):
        return {"exists": False}
    segs = {}
    n_links = n_paths = n_pair_paths = 0
    path_recs = []
    with open(gfa_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("H\t"):
                continue
            parts = line.split("\t")
            rec = parts[0]
            if rec == "S":
                segs[parts[1]] = parts[2]
            elif rec == "L":
                n_links += 1
            elif rec == "P":
                n_paths += 1
                path_recs.append((parts[1], parts[2]))
                if "|pair|" in parts[1]:
                    n_pair_paths += 1

    # Verify every single candidate path in flanks.json matches its GFA P line byte-for-byte
    path_map = {}
    for pname, psegs in path_recs:
        seq_pieces = []
        for tok in psegs.split(","):
            sid, orient = tok[:-1], tok[-1]
            s = segs[sid]
            seq_pieces.append(s if orient == "+" else sc.rc(s))
        path_map[pname] = "".join(seq_pieces)

    def sanitize(s):
        return "".join("_" if ch in " \t,;+-" else ch for ch in s)

    checked = mismatches = 0
    raw_bp = 0
    for gname, gentry in doc_genes.items():
        gseq = genes_dict.get(gname, "")
        if not gseq:
            continue
        gid = sanitize(gname)
        raw_bp += len(gseq)
        for p in gentry.get("left", {}).get("paths", []):
            raw_bp += len(p["flank"])
            checked += 1
            exp = p["flank"] + gseq
            got = path_map.get(f"{gid}|left|{p['id']}")
            if got != exp:
                mismatches += 1
        for p in gentry.get("right", {}).get("paths", []):
            raw_bp += len(p["flank"])
            checked += 1
            exp = gseq + p["flank"]
            got = path_map.get(f"{gid}|right|{p['id']}")
            if got != exp:
                mismatches += 1
    comp_bp = sum(len(s) for s in segs.values())
    return {
        "exists": True,
        "genes": len(doc_genes),
        "segments": len(segs),
        "links": n_links,
        "paths": n_paths,
        "pair_paths": n_pair_paths,
        "checked_paths": checked,
        "roundtrip_mismatches": mismatches,
        "raw_bp": raw_bp,
        "compacted_bp": comp_bp,
        "compression_ratio": round(raw_bp / max(1, comp_bp), 2),
    }


def naive_pairs(gentry, W, mode):
    l_tr = [t for t in extract_trunks(gentry.get("left", {}).get("paths", []), "left", W=W) if t["certain_bp"] >= W and t["shared_end"] == 0]
    r_tr = [t for t in extract_trunks(gentry.get("right", {}).get("paths", []), "right", W=W) if t["certain_bp"] >= W and t["shared_end"] == 0]
    if not l_tr or not r_tr:
        return []
    out = []
    if mode == "cartesian":
        idx = 0
        for lt in l_tr:
            for rt in r_tr:
                idx += 1
                out.append({"id": f"c{idx}", "left_prox": lt["prox"], "right_prox": rt["prox"], "tier": "cartesian",
                            "left_depth": lt["depth"], "right_depth": rt["depth"], "depth_source": "na"})
    elif mode == "main":
        lt = max(l_tr, key=lambda t: (t["certain_bp"], t["depth"]))
        rt = max(r_tr, key=lambda t: (t["certain_bp"], t["depth"]))
        out.append({"id": "m1", "left_prox": lt["prox"], "right_prox": rt["prox"], "tier": "main",
                    "left_depth": lt["depth"], "right_depth": rt["depth"], "depth_source": "na"})
    elif mode == "greedy":
        ls = sorted(l_tr, key=lambda t: -t["depth"])
        rs = sorted(r_tr, key=lambda t: -t["depth"])
        for idx, (lt, rt) in enumerate(zip(ls, rs), 1):
            out.append({"id": f"g{idx}", "left_prox": lt["prox"], "right_prox": rt["prox"], "tier": "greedy",
                        "left_depth": lt["depth"], "right_depth": rt["depth"], "depth_source": "na"})
    return out


def score_pairing_mode(doc_genes, genes, fam_loci, preds_by_gene, W, min_cov=5.0, perm_genomes=None):
    perm_genomes = perm_genomes or set()
    tot = chk = ok = bad = unc = 0
    bad_strict97 = 0
    by_tier = collections.defaultdict(lambda: [0, 0])
    post_fork_used = 0
    bad_rows = []
    seen_fam_pairs = {}

    for gname, plist in sorted(preds_by_gene.items()):
        fam = family(gname)
        floci = fam_loci.get(fam, [])
        for pr in plist:
            tot += 1
            lp, rp = pr["left_prox"], pr["right_prox"]
            if len(lp) < W or len(rp) < W:
                unc += 1
                continue
            if pr.get("depth_source") == "post_fork_bins":
                post_fork_used += 1
            has_both = False
            has_both_97 = False
            has_l_contradict = False
            has_r_contradict = False
            for loc in floci:
                ml = match_flank_w(lp, loc["left_5to3"], "left", W, strict=False) if len(loc["left_5to3"]) >= W else False
                mr = match_flank_w(rp, loc["right_5to3"], "right", W, strict=False) if len(loc["right_5to3"]) >= W else False
                if ml and mr:
                    has_both = True
                    ml97 = match_flank_w(lp, loc["left_5to3"], "left", W, strict=True)
                    mr97 = match_flank_w(rp, loc["right_5to3"], "right", W, strict=True)
                    if ml97 and mr97:
                        has_both_97 = True
                elif ml and len(loc["right_5to3"]) >= W:
                    has_l_contradict = True
                elif mr and len(loc["left_5to3"]) >= W:
                    has_r_contradict = True

            if has_both:
                verdict = "confirmed"
                chk += 1
                ok += 1
                if not has_both_97:
                    bad_strict97 += 1
                by_tier[pr["tier"]][1] += 1
            elif has_l_contradict or has_r_contradict:
                verdict = "chimera"
                chk += 1
                bad += 1
                by_tier[pr["tier"]][0] += 1
                by_tier[pr["tier"]][1] += 1
                bad_rows.append({
                    "gene": gname, "pair_id": pr["id"], "tier": pr["tier"],
                    "dL": pr.get("left_depth", 0), "dR": pr.get("right_depth", 0),
                    "note": pr.get("note", ""), "l_contra": has_l_contradict, "r_contra": has_r_contradict,
                })
            else:
                verdict = "uncheckable"
                unc += 1

            sig = (fam, lp[-120:], rp[:120])
            if verdict in ("confirmed", "chimera"):
                if sig not in seen_fam_pairs or verdict == "confirmed":
                    seen_fam_pairs[sig] = verdict

    fam_chk = len(seen_fam_pairs)
    fam_bad = sum(1 for v in seen_fam_pairs.values() if v == "chimera")

    # Recall over distinct true two-sided contexts (cov >= min_cov, both sides >= W bp)
    true_all = true_hit = true_multi = true_multi_hit = 0
    fam_genes = collections.defaultdict(list)
    for g in doc_genes:
        fam_genes[family(g)].append(g)

    for fam, glist in sorted(fam_genes.items()):
        # Deduplicate rotations and near-identical true two-sided molecules (>= min_cov)
        seen_mol = set()
        clusters = []
        for loc in fam_loci.get(fam, []):
            if loc["cov"] < min_cov or len(loc["left_5to3"]) < W or len(loc["right_5to3"]) < W:
                continue
            seq_norm = re.sub(r"^\d", "", loc["seq"]) if loc["genome"] in perm_genomes else loc["seq"]
            mkey = (loc["genome"], seq_norm, loc["left_5to3"][-60:], loc["right_5to3"][:60])
            if mkey in seen_mol:
                continue
            seen_mol.add(mkey)
            lw, rw = loc["left_5to3"][-W:], loc["right_5to3"][:W]
            if not any(match_flank_w(lw, cl[0], "left", W) and match_flank_w(rw, cl[1], "right", W) for cl in clusters):
                clusters.append((lw, rw))
        if not clusters:
            continue
        is_multi = len(clusters) >= 2
        fam_preds = [pr for g in glist for pr in preds_by_gene.get(g, []) if len(pr["left_prox"]) >= W and len(pr["right_prox"]) >= W]
        for lw, rw in clusters:
            true_all += 1
            if is_multi:
                true_multi += 1
            if any(match_flank_w(pr["left_prox"], lw, "left", W) and match_flank_w(pr["right_prox"], rw, "right", W) for pr in fam_preds):
                true_hit += 1
                if is_multi:
                    true_multi_hit += 1

    lo, hi = cp_ci(bad, chk)
    return {
        "W": W,
        "emitted": tot,
        "checkable": chk,
        "confirmed": ok,
        "chimeras": bad,
        "uncheckable": unc,
        "chimera_rate_pct": round(100.0 * bad / chk, 2) if chk else 0.0,
        "ci95_pct": [round(100.0 * lo, 2), round(100.0 * hi, 2)],
        "fam_checkable": fam_chk,
        "fam_chimeras": fam_bad,
        "non_97_confirmed": bad_strict97,
        "true_all": true_all,
        "true_hit": true_hit,
        "recall_pct": round(100.0 * true_hit / true_all, 1) if true_all else 0.0,
        "true_multi": true_multi,
        "true_multi_hit": true_multi_hit,
        "multi_recall_pct": round(100.0 * true_multi_hit / true_multi, 1) if true_multi else 0.0,
        "by_tier": {k: {"bad": v[0], "chk": v[1]} for k, v in sorted(by_tier.items())},
        "post_fork_bins_pairs": post_fork_used,
        "bad_details": bad_rows,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("genes")
    ap.add_argument("run_prefix")
    ap.add_argument("genomes")
    ap.add_argument("genomes_tsv")
    ap.add_argument("cov_col")
    ap.add_argument("blast_tsv")
    ap.add_argument("out_prefix")
    ap.add_argument("--min-cov", type=float, default=5.0)
    ap.add_argument("--blastdb", default=None)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--blastn", default="blastn")
    ap.add_argument("--reads-tsv", default=None)
    ap.add_argument("--label", default="")
    a = ap.parse_args()

    genes = sc.read_fa(a.genes)
    glen = {g: len(s) for g, s in genes.items()}
    doc = json.load(open(a.run_prefix + ".flanks.json"))
    doc_genes = doc["genes"]
    D = doc["flank_bp"]

    # Coverage + permuted-plasmid flags
    cov, perm_genomes = {}, set()
    with open(a.genomes_tsv) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        cidx = hdr.index(a.cov_col)
        fidx = hdr.index("file") if "file" in hdr else -1
        for line in f:
            p = line.rstrip("\n").split("\t")
            g = p[0]
            cov[g] = float(p[cidx])
            if fidx >= 0 and "permuted" in p[fidx]:
                perm_genomes.add(g)

    loci = sc.loci_from_blast(a.blast_tsv, glen)
    need_seqs = {s for v in loci.values() for s, *_ in v}
    gseq = sc.read_fa(a.genomes, need_seqs)

    # Build true contexts per gene AND pooled per family (for two-sided & strict multi-allele scoring)
    gene_ctx = {}
    fam_loci = collections.defaultdict(list)
    seen_fam_loc = set()
    for g in genes:
        ctxs = []
        for s, pid, qs, qe, ss, se in loci.get(g, []):
            genome, sid = s.split("|", 1)
            cv = cov.get(genome, 0.0)
            if cv <= 0:
                continue
            left_out, right_out, g0, g1 = sc.project(gseq[s], glen[g], qs, qe, ss, se, D + 500)
            gs = gseq[s]
            if ss <= se:
                lpad_out = gs[max(0, g0):min(len(gs), g0 + PADG)][::-1]
                rpad_out = gs[max(0, g1 - PADG):min(len(gs), g1)]
            else:
                lpad_out = sc.rc(gs[max(0, g1 - PADG):min(len(gs), g1)])[::-1]
                rpad_out = sc.rc(gs[max(0, g0):min(len(gs), g0 + PADG)])
            entry = {
                "genome": genome, "seq": sid, "cov": cv, "pid": pid,
                "left": left_out[:D], "right": right_out[:D],
                "left_full": left_out, "right_full": right_out,
                "left_5to3": left_out[::-1], "right_5to3": right_out,
                "padded": {
                    "left": (lpad_out + left_out, len(lpad_out)),
                    "right": (rpad_out + right_out, len(rpad_out)),
                },
                "g0": g0, "g1": g1, "strand": "+" if ss <= se else "-",
            }
            ctxs.append(entry)
            fk = (family(g), genome, sid, g0, g1, entry["strand"])
            if fk not in seen_fam_loc:
                seen_fam_loc.add(fk)
                fam_loci[family(g)].append(entry)
        gene_ctx[g] = ctxs

    # =========================================================================
    # 1. STRICT ALIGNMENT CHECK OF ONE-SIDED CERTIFIED CLAIMS (`claims_lr3.py` v3.1 rules)
    # =========================================================================
    # Build target FASTA of all claimed intervals (`context_len >= 100` and `end_len > context_len`)
    # and query FASTA of true genome gene-flank neighborhoods (cov > 0), run minimap2 once!
    claims = []
    with tempfile.TemporaryDirectory(prefix="argj_score_") as td:
        ref_fa = os.path.join(td, "claims.fa")
        qry_fa = os.path.join(td, "truth.fa")
        paf_file = os.path.join(td, "aln.paf")

        claim_info = {}
        with open(ref_fa, "w") as rf:
            for g in genes:
                gseq_s = genes[g]
                for side in ("left", "right"):
                    paths = doc_genes.get(g, {}).get(side, {}).get("paths", [])
                    for idx, p in enumerate(paths, 1):
                        flk = p["flank"]
                        if not flk:
                            continue
                        clen = min(int(p.get("context_len", 0)), len(flk))
                        elen = min(int(p.get("end_len", clen)), len(flk))
                        pad = min(PADG, len(gseq_s))
                        if clen >= 100:
                            cid = f"{g}|{side}|c{idx}|cert"
                            ref_seq = (flk[-clen:] + gseq_s[:pad]) if side == "left" else (gseq_s[-pad:] + flk[:clen])
                            rf.write(f">{cid}\n{ref_seq}\n")
                            claim_info[cid] = {
                                "gene": g, "fam": family(g), "side": side, "cand": idx, "kind": "certain",
                                "claim": clen, "flank_len": len(flk), "pad": pad, "depth": p.get("depth", 0.0),
                                "outward": (flk[::-1] if side == "left" else flk)[:clen],
                                "seq_5to3": flk[-clen:] if side == "left" else flk[:clen],
                                "stop": p.get("certain_stop", ""), "shared_end": p.get("shared_end", 0),
                                "gene_end": p.get("gene_end", ""),
                            }
                        if elen >= 100 and elen > clen:
                            eid = f"{g}|{side}|c{idx}|end"
                            ref_seq = (flk[-elen:] + gseq_s[:pad]) if side == "left" else (gseq_s[-pad:] + flk[:elen])
                            rf.write(f">{eid}\n{ref_seq}\n")
                            claim_info[eid] = {
                                "gene": g, "fam": family(g), "side": side, "cand": idx, "kind": "end_only",
                                "claim": elen, "flank_len": len(flk), "pad": pad, "depth": p.get("depth", 0.0),
                                "outward": (flk[::-1] if side == "left" else flk)[:elen],
                                "seq_5to3": flk[-elen:] if side == "left" else flk[:elen],
                                "stop": p.get("end_stop", ""), "shared_end": p.get("shared_end", 0),
                                "gene_end": p.get("gene_end", ""),
                            }

        # Write all true genome neighborhoods (cov > 0) grouped by gene family
        qry_meta = {}
        with open(qry_fa, "w") as qf:
            q_idx = 0
            for fam, flist in fam_loci.items():
                for loc in flist:
                    q_idx += 1
                    for side in ("left", "right"):
                        pseq, plen = loc["padded"][side]
                        # Convert outward `pseq` (gene-pad + flank) to 5'->3' orientation
                        seq_5to3 = pseq[::-1] if side == "left" else pseq
                        qid = f"T{q_idx}|{fam}|{side}|{loc['genome']}|{loc['seq']}"
                        qf.write(f">{qid}\n{seq_5to3}\n")
                        qry_meta[qid] = {
                            "fam": fam, "side": side, "genome": loc["genome"], "seq": loc["seq"],
                            "cov": loc["cov"], "flank_len": len(pseq) - plen,
                        }

        alns_by_claim = collections.defaultdict(list)
        if claim_info and qry_meta:
            subprocess.run(
                [a.minimap2, "-c", "-x", "map-hifi", "--eqx", "-N", "2000", "-p", "0.001", "-t", "8", ref_fa, qry_fa],
                stdout=open(paf_file, "w"), stderr=subprocess.DEVNULL, check=True,
            )
            with open(paf_file) as pf:
                for line in pf:
                    p = line.rstrip("\n").split("\t")
                    qid, tname = p[0], p[5]
                    if tname not in claim_info or qid not in qry_meta:
                        continue
                    ci = claim_info[tname]
                    qm = qry_meta[qid]
                    # Only align truth molecules of the SAME gene family and SAME side!
                    if qm["fam"] != ci["fam"] or qm["side"] != ci["side"]:
                        continue
                    qlen, qs, qe, strand, rs, re_ = int(p[1]), int(p[2]), int(p[3]), p[4], int(p[7]), int(p[8])
                    if strand != "+":
                        continue
                    cg = next((t[5:] for t in p[12:] if t.startswith("cg:Z:")), None)
                    if cg is None:
                        continue
                    pad, cl, right = ci["pad"], ci["claim"], ci["side"] == "right"
                    diffs, longs, r = [], [], rs
                    for n_s, op in CIG.findall(cg):
                        n = int(n_s)
                        if op in "=M":
                            r += n
                        elif op == "X":
                            for b in range(r, r + n):
                                diffs.append(((b - pad) if right else (cl - 1 - b), 1))
                            r += n
                        elif op == "D":
                            o = (r - pad) if right else (cl - r - n)
                            diffs.append((o, n))
                            if n >= 50:
                                longs.append((o, n))
                            r += n
                        elif op == "I":
                            o = (r - pad) if right else (cl - r)
                            diffs.append((o, n))
                            if n >= 50:
                                longs.append((o, n))
                    if right:
                        alns_by_claim[tname].append(Aln(rs - pad, re_ - pad, qlen, qs, qe, diffs, longs, qid))
                    else:
                        alns_by_claim[tname].append(Aln(cl - re_, cl - rs, qlen, qs, qe, diffs, longs, qid))

        # Score every claim
        claim_rows = []
        for cid, ci in claim_info.items():
            al = alns_by_claim.get(cid, [])
            cl = ci["claim"]
            ok_strict, hit_q = check_claim_strict(al, cl, start_x=-GENE_IN, thr=50, div=0.03, sub=STRICT_SUB)
            ok_mol, _ = check_claim_strict(al, cl, start_x=-50, thr=50, div=0.03, sub=MOL_SUB) if ok_strict else (False, None)
            ok_lenient, _ = check_claim_strict(al, cl, start_x=-GENE_IN, thr=3000, div=0.15, sub=None) if not ok_strict else (True, hit_q)
            ok_flank_only, _ = check_claim_strict(al, cl, start_x=0, thr=50, div=0.03, sub=STRICT_SUB) if not ok_strict else (True, hit_q)
            # Check if a true contig ends before `cl` while matching strictly up to the contig end
            contig_end_hit = False
            if not ok_strict:
                for z in al:
                    qm = qry_meta[z.qid]
                    if qm["flank_len"] < cl - EDGE and z.o1 >= qm["flank_len"] - EDGE:
                        if qm["flank_len"] >= 100:
                            ok_ce, _ = check_claim_strict([z], qm["flank_len"], start_x=-GENE_IN, thr=50, div=0.03, sub=STRICT_SUB)
                            if ok_ce:
                                contig_end_hit = True
                                break
            if ok_strict:
                verdict = "passed"
                triage = "molecule_exact" if ok_mol else "strain_mosaic"
            elif ok_lenient:
                verdict = "differs"
                triage = "gene_end_diff" if ok_flank_only else "lenient_strain_diff"
            else:
                verdict = "contradicted"
                if ok_flank_only:
                    triage = "gene_end_diff"
                elif contig_end_hit:
                    triage = "past_contig_end"
                else:
                    triage = "unexplained"

            claim_rows.append({
                "cid": cid, "gene": ci["gene"], "fam": ci["fam"], "side": ci["side"], "cand": ci["cand"],
                "kind": ci["kind"], "claim": cl, "flank_len": ci["flank_len"], "depth": ci["depth"],
                "verdict": verdict, "triage": triage, "mol_exact": int(ok_mol),
                "hit_truth": hit_q or "-", "stop": ci["stop"], "shared_end": ci["shared_end"],
                "elsewhere_blast": "", "_out": ci["outward"], "_5to3": ci["seq_5to3"],
            })

        # BLAST any failing claims against ALL present genomes (`cov > 0`) to check for unseeded homologous loci
        failing = [r for r in claim_rows if r["verdict"] != "passed"]
        if failing and a.blastdb:
            bf = os.path.join(td, "fail_q.fa")
            with open(bf, "w") as f:
                for i, r in enumerate(failing):
                    f.write(f">F{i}\n{r['_5to3']}\n")
            bout = subprocess.run(
                [a.blastn, "-query", bf, "-db", a.blastdb, "-perc_identity", "97", "-evalue", "1e-5",
                 "-dust", "no", "-max_target_seqs", "2000", "-num_threads", "8",
                 "-outfmt", "6 qseqid sseqid pident length qlen qstart qend"],
                capture_output=True, text=True, check=True,
            ).stdout
            best_hit = {}
            for line in bout.splitlines():
                q, s, pid, alen, qlen, qs, qe = line.split("\t")
                gn = s.split("|", 1)[0]
                if cov.get(gn, 0.0) <= 0:
                    continue
                qcov = (int(qe) - int(qs) + 1) / max(1, int(qlen))
                pid_f = float(pid)
                prev = best_hit.get(q)
                if prev is None or (qcov, pid_f) > (prev[0], prev[1]):
                    best_hit[q] = (qcov, pid_f, gn, s, round(cov[gn], 1))
            for i, r in enumerate(failing):
                bh = best_hit.get(f"F{i}")
                if bh:
                    r["elsewhere_blast"] = f"{bh[2]}({bh[4]}x):{100*bh[0]:.1f}%cov@{bh[1]:.2f}%id"
                    if bh[0] >= 0.98 and bh[1] >= 99.0 and r["triage"] == "unexplained":
                        r["triage"] = "elsewhere_in_sample_99pct"

    # Independent 31-mer sequence groups over certified (`kind == 'certain'`) claims
    def assign_groups(sub_rows):
        k2c = {}
        parent = list(range(len(sub_rows)))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for i, r in enumerate(sub_rows):
            o = r["_out"]
            for j in range(0, len(o) - K31 + 1):
                km = o[j:j + K31]
                o2 = k2c.get(km)
                if o2 is None:
                    k2c[km] = i
                else:
                    a_, b_ = find(i), find(o2)
                    if a_ != b_:
                        parent[a_] = b_
        return [find(i) for i in range(len(sub_rows))]

    cert_rows = [r for r in claim_rows if r["kind"] == "certain"]
    end_rows = [r for r in claim_rows if r["kind"] == "end_only"]

    strict_summary = {}
    for label_k, pool in [("certain", cert_rows), ("end_only", end_rows)]:
        strict_summary[label_k] = {}
        for lim in (100, 1000, 2000, 5000):
            sub = [r for r in pool if r["claim"] >= lim]
            n = len(sub)
            n_pass = sum(1 for r in sub if r["verdict"] == "passed")
            n_diff = sum(1 for r in sub if r["verdict"] == "differs")
            n_con = sum(1 for r in sub if r["verdict"] == "contradicted")
            n_fail = n_diff + n_con
            n_mol = sum(r["mol_exact"] for r in sub)
            n_true_fault = sum(1 for r in sub if r["verdict"] != "passed" and r["triage"] == "unexplained")
            gids = assign_groups(sub)
            g_fail = collections.defaultdict(int)
            g_fault = collections.defaultdict(int)
            for gid, r in zip(gids, sub):
                g_fail[gid] += int(r["verdict"] != "passed")
                g_fault[gid] += int(r["verdict"] != "passed" and r["triage"] == "unexplained")
            ng = len(g_fail)
            ng_fail = sum(1 for v in g_fail.values() if v > 0)
            ng_fault = sum(1 for v in g_fault.values() if v > 0)
            lo_c, hi_c = cp_ci(n_fail, n)
            lo_g, hi_g = cp_ci(ng_fail, ng)
            strict_summary[label_k][f">={lim}"] = {
                "claims": n, "passed": n_pass, "differs": n_diff, "contradicted": n_con,
                "failed": n_fail, "unexplained_faults": n_true_fault,
                "fail_rate_pct": round(100.0 * n_fail / n, 3) if n else 0.0,
                "claim_ci95_pct": [round(100.0 * lo_c, 2), round(100.0 * hi_c, 2)],
                "groups": ng, "failed_groups": ng_fail, "fault_groups": ng_fault,
                "group_ci95_pct": [round(100.0 * lo_g, 2), round(100.0 * hi_g, 2)],
                "molecule_exact": n_mol,
            }

    # =========================================================================
    # 2. STRUCTURAL RECALL PER TRUE CONTEXT (cov >= min_cov; raw & rotation-collapsed)
    # =========================================================================
    ctx_rows = []
    for g in genes:
        ctxs = gene_ctx[g]
        for side in ("left", "right"):
            paths = doc_genes.get(g, {}).get(side, {}).get("paths", [])
            cands = []
            for i, p in enumerate(paths, 1):
                f = p["flank"][::-1] if side == "left" else p["flank"]
                if not f:
                    continue
                cands.append({"i": i, "f": f, "ctx": int(p.get("context_len", 0)), "end": int(p.get("end_len", p.get("context_len", 0)))})
            T = [c[side] for c in ctxs]

            def shifted(f, j):
                tp, pl = ctxs[j]["padded"][side]
                probe = f[:24]
                if len(probe) < 24:
                    return T[j]
                best, st = None, tp.find(probe, max(0, pl - 60), pl + 84)
                while st != -1:
                    if best is None or abs(st - pl) < abs(best - pl):
                        best = st
                    st = tp.find(probe, st + 1, pl + 84)
                return tp[best:best + D] if best is not None else T[j]

            SL = []
            TL = []
            for c in cands:
                sl, tl = [], []
                for j in range(len(T)):
                    t = shifted(c["f"], j)
                    sl.append(struct_len(c["f"], kindex(t), len(t), 3000))
                    tl.append(len(t))
                SL.append(sl)
                TL.append(tl)

            for j, c in enumerate(ctxs):
                if c["cov"] < a.min_cov:
                    continue
                ach = min(D, len(T[j]))
                nC = range(len(cands))
                s_bp = min(ach, max((SL[ci][j] for ci in nC), default=0))
                cert_bp = min(ach, max((min(cands[ci]["ctx"], SL[ci][j]) for ci in nC if SL[ci][j] >= min(cands[ci]["ctx"], TL[ci][j]) - TOL_STRUCT), default=0))
                end_bp = min(ach, max((min(cands[ci]["end"], SL[ci][j]) for ci in nC if SL[ci][j] >= min(cands[ci]["end"], TL[ci][j]) - TOL_STRUCT), default=0))
                ctx_rows.append({
                    "gene": g, "side": side, "genome": c["genome"], "seq": c["seq"], "cov": round(c["cov"], 2),
                    "achievable": ach, "struct_bp": s_bp, "cert_bp": cert_bp, "end_bp": end_bp,
                    "found": int(s_bp >= min(1000, ach)),
                })

    def summarize_recall(rows):
        n = len(rows)
        reach = lambda k, L: sum(1 for r in rows if r[k] >= min(L, r["achievable"]))
        return {
            "n_contexts": n,
            "found_1k": sum(r["found"] for r in rows),
            "cert_1k": reach("cert_bp", 1000),
            "cert_2k": reach("cert_bp", 2000),
            "cert_5k": reach("cert_bp", 5000),
            "cert_full": reach("cert_bp", D),
            "end_1k": reach("end_bp", 1000),
            "end_5k": reach("end_bp", 5000),
            "mean_cert_bp": int(statistics.mean(r["cert_bp"] for r in rows)) if rows else 0,
            "median_cert_bp": int(statistics.median(r["cert_bp"] for r in rows)) if rows else 0,
            "mean_achievable_bp": int(statistics.mean(r["achievable"] for r in rows)) if rows else 0,
        }

    # Rotation-collapsed recall (important for CAMI III circular plasmids)
    col_grp = collections.defaultdict(list)
    for r in ctx_rows:
        seq_n = re.sub(r"^\d", "", r["seq"]) if r["genome"] in perm_genomes else r["seq"]
        col_grp[(r["gene"], r["side"], r["genome"], seq_n)].append(r)
    col_rows = [
        {
            "gene": k[0], "side": k[1], "genome": k[2], "seq": k[3],
            "achievable": max(x["achievable"] for x in xs),
            "struct_bp": max(x["struct_bp"] for x in xs),
            "cert_bp": max(x["cert_bp"] for x in xs),
            "end_bp": max(x["end_bp"] for x in xs),
            "found": max(x["found"] for x in xs),
        }
        for k, xs in col_grp.items()
    ]

    # =========================================================================
    # 3. GFA v1.0 LOSSLESS COMPACTION VERIFICATION
    # =========================================================================
    gfa_res = verify_gfa_file(a.run_prefix + ".gfa", genes, doc_genes)

    # =========================================================================
    # 4. SINGLE-SAMPLE TWO-SIDED PAIRING (`v0.5.0` WITH `depth_bins` vs BASELINES)
    # =========================================================================
    reads_src = a.reads_tsv
    if not reads_src:
        for cand_r in (a.run_prefix + ".tmp/aw.reads.tsv", a.run_prefix + ".tmp/aw.reads.tsv.gz"):
            if os.path.isfile(cand_r):
                reads_src = cand_r
                break

    pairing_res = {}
    for W in (500, 1000):
        # Load actual emitted pairs at W=500 from flanks.json/pairs.tsv (adding prox seq + depth_source),
        # and run pair_gene_single_sample at W=1000 (with reads_src if present)
        phased = index_sample_reads(genes, doc_genes, reads_src, W=W) if reads_src else {}
        cert_preds = {}
        for gname, gentry in doc_genes.items():
            gseq_s = genes.get(gname)
            if not gseq_s:
                continue
            l_tr = {t["trunk_id"]: t for t in extract_trunks(gentry.get("left", {}).get("paths", []), "left", W=W)}
            r_tr = {t["trunk_id"]: t for t in extract_trunks(gentry.get("right", {}).get("paths", []), "right", W=W)}
            if W == 500 and "pairs" in gentry:
                # Score the ACTUAL pairs emitted by `argjacent v0.5.0`!
                plist = []
                for pr in gentry["pairs"]:
                    lt = l_tr.get(pr["left_trunk"])
                    rt = r_tr.get(pr["right_trunk"])
                    if lt and rt:
                        item = dict(pr)
                        item["left_prox"] = lt["prox"]
                        item["right_prox"] = rt["prox"]
                        item["depth_source"] = "post_fork_bins" if (lt["depth_source"] == "post_fork_bins" or rt["depth_source"] == "post_fork_bins") else "boundary"
                        plist.append(item)
                if plist:
                    cert_preds[gname] = plist
            else:
                prs = pair_gene_single_sample(gname, gseq_s, gentry, W=W, pre_phased=phased.get(gname, {}), include_prox=True)
                for pr in prs:
                    lt = l_tr.get(pr["left_trunk"])
                    rt = r_tr.get(pr["right_trunk"])
                    pr["depth_source"] = "post_fork_bins" if ((lt and lt["depth_source"] == "post_fork_bins") or (rt and rt["depth_source"] == "post_fork_bins")) else "boundary"
                if prs:
                    cert_preds[gname] = prs

        w_out = {"argjacent_v05": score_pairing_mode(doc_genes, genes, fam_loci, cert_preds, W, a.min_cov, perm_genomes)}
        for bname in ("cartesian", "main", "greedy"):
            bpreds = {g: naive_pairs(ge, W, bname) for g, ge in doc_genes.items()}
            w_out[bname] = score_pairing_mode(doc_genes, genes, fam_loci, bpreds, W, a.min_cov, perm_genomes)
        pairing_res[f"W{W}"] = w_out

    # Save detailed tables & JSON summary
    with open(a.out_prefix + ".claims.tsv", "w") as f:
        cols = ["gene", "side", "cand", "kind", "claim", "flank_len", "depth", "verdict", "triage", "mol_exact", "hit_truth", "stop", "shared_end", "elsewhere_blast"]
        f.write("\t".join(cols) + "\n")
        for r in claim_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(a.out_prefix + ".failures.tsv", "w") as f:
        cols = ["gene", "side", "cand", "kind", "claim", "flank_len", "depth", "verdict", "triage", "stop", "shared_end", "elsewhere_blast"]
        f.write("\t".join(cols) + "\n")
        for r in claim_rows:
            if r["verdict"] != "passed":
                f.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(a.out_prefix + ".recall.tsv", "w") as f:
        cols = ["gene", "side", "genome", "seq", "cov", "achievable", "struct_bp", "cert_bp", "end_bp", "found"]
        f.write("\t".join(cols) + "\n")
        for r in ctx_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    summary = {
        "label": a.label or os.path.basename(a.out_prefix),
        "n_genes": len(genes),
        "n_genes_in_json": len(doc_genes),
        "strict_one_sided": strict_summary,
        "recall_raw": summarize_recall(ctx_rows),
        "recall_collapsed": summarize_recall(col_rows),
        "gfa": gfa_res,
        "pairing": pairing_res,
    }
    with open(a.out_prefix + ".v05_score.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Compact human-readable stdout report
    c1k = strict_summary["certain"][">=1000"]
    c100 = strict_summary["certain"][">=100"]
    c5k = strict_summary["certain"][">=5000"]
    rc_s = summary["recall_collapsed"]
    p5 = pairing_res["W500"]["argjacent_v05"]
    p5c = pairing_res["W500"]["cartesian"]
    print(f"=== {summary['label']} ({len(genes)} seed genes) ===")
    print(f"  1-sided strict >=1kb : {c1k['failed']}/{c1k['claims']} failed ({c1k['unexplained_faults']} faults) | groups: {c1k['failed_groups']}/{c1k['groups']} | mol-exact: {c1k['molecule_exact']}/{c1k['claims']}")
    print(f"  1-sided strict >=100 : {c100['failed']}/{c100['claims']} failed ({c100['unexplained_faults']} faults) | >=5kb: {c5k['failed']}/{c5k['claims']} failed")
    print(f"  Recall (>=5x, col.)  : {rc_s['n_contexts']} contexts | found {rc_s['found_1k']} | cert>=1k {rc_s['cert_1k']} | cert>=5k {rc_s['cert_5k']} | full {rc_s['cert_full']} (mean cert {rc_s['mean_cert_bp']} bp / ach {rc_s['mean_achievable_bp']} bp)")
    print(f"  GFA v1.0 lossless    : {gfa_res.get('segments',0)} segs, {gfa_res.get('paths',0)} paths, {gfa_res.get('compression_ratio',0)}x compact, roundtrip_mismatches={gfa_res.get('roundtrip_mismatches','?')}")
    print(f"  2-sided pair (W=500) : ARGjacent {p5['chimeras']}/{p5['checkable']} chimeras ({p5['recall_pct']}% recall, {p5['multi_recall_pct']}% multi, post_fork_bins={p5['post_fork_bins_pairs']}) vs Cartesian {p5c['chimeras']}/{p5c['checkable']}")


if __name__ == "__main__":
    main()
