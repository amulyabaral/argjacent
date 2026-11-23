#!/usr/bin/env python3
"""Score multi-sample `argjacent-cohort` outputs (Part A pairing + Part B post-repeat exits)
on CAMI II / III cohorts against exact source genomes.

Runs the shipped `bin/argjacent-cohort` CLI directly on the retained `flanks.json` files
(JSON-only mode using `v0.5.0` `depth_bins` + `gene_depth`, zero raw reads kept!) and verifies:
  - Part A (`cohort_pairs.tsv`): every emitted two-sided pair across all samples, comparing
    single-sample vs cohort tiers (`cohort_flow`, `cohort_anchor`, `cohort_complement`),
    chimera rate, and recall on multi-context loci.
  - Part B (`cohort_repeats.tsv`): every post-repeat certified exit arm (`cohort_certified_bp`),
    verified strictly (`claims_lr3.py` v3.1 rules) against the exact source genomes present
    in that sample (`cov > 0`)!
"""
import argparse
import collections
import json
import os
import re
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "benchmarks", "cami"))
sys.path.insert(0, os.path.join(REPO, "benchmarks", "v05_validation"))

import score_cami as sc
from argjacent.pair import extract_trunks
from score_v05_cami import (
    Aln, CIG, EDGE, GENE_IN, MOL_SUB, PADG, STRICT_SUB,
    check_claim_strict, cp_ci, family, score_pairing_mode,
)


def load_sample_truth(genes_fa, genomes_fa, genomes_tsv, cov_col, blast_tsv, D=10000):
    genes = sc.read_fa(genes_fa)
    glen = {g: len(s) for g, s in genes.items()}
    cov, perm = {}, set()
    with open(genomes_tsv) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        cidx = hdr.index(cov_col)
        fidx = hdr.index("file") if "file" in hdr else -1
        for line in f:
            p = line.rstrip("\n").split("\t")
            cov[p[0]] = float(p[cidx])
            if fidx >= 0 and "permuted" in p[fidx]:
                perm.add(p[0])
    loci = sc.loci_from_blast(blast_tsv, glen)
    gseq = sc.read_fa(genomes_fa, {s for v in loci.values() for s, *_ in v})
    fam_loci = collections.defaultdict(list)
    seen = set()
    for g in genes:
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
            st = "+" if ss <= se else "-"
            fk = (family(g), genome, sid, g0, g1, st)
            if fk not in seen:
                seen.add(fk)
                fam_loci[family(g)].append({
                    "genome": genome, "seq": sid, "cov": cv,
                    "left_5to3": left_out[::-1], "right_5to3": right_out,
                    "padded": {
                        "left": (lpad_out + left_out, len(lpad_out)),
                        "right": (rpad_out + right_out, len(rpad_out)),
                    },
                })
    return genes, cov, perm, fam_loci


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="JSON list of sample dicts: sid, genes_fa, flanks_json, genomes_fa, genomes_tsv, cov_col, blast_tsv")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--pair-bp", type=int, default=500)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--label", default="cohort")
    a = ap.parse_args()

    manifest = json.load(open(a.manifest))
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)), exist_ok=True)
    sheet_path = a.out_prefix + ".sample_sheet.tsv"
    with open(sheet_path, "w") as f:
        f.write("sample\tflanks_json\tgenes_fa\n")
        for m in manifest:
            f.write(f"{m['sid']}\t{m['flanks_json']}\t{m['genes_fa']}\n")

    cohort_bin = os.path.join(REPO, "bin", "argjacent-cohort")
    subprocess.run(
        [cohort_bin, "--sample-sheet", sheet_path, "-o", a.out_prefix, "--pair-bp", str(a.pair_bp)],
        check=True,
    )
    cohort_doc = json.load(open(a.out_prefix + ".cohort.json"))

    # Score Part A (cohort two-sided pairing) AND Part B (post-repeat exits) across all samples
    tot_chk = tot_bad = tot_ok = tot_unc = 0
    tot_true = tot_hit = tot_multi = tot_multi_hit = 0
    tot_fam_chk = tot_fam_bad = 0
    by_tier = collections.defaultdict(lambda: [0, 0])
    cohort_only_chk = cohort_only_bad = 0
    bad_pairs = []

    rep_total = rep_pass = rep_fail = 0
    rep_gains = []
    rep_rows = []

    for m in manifest:
        sid = m["sid"]
        doc = json.load(open(m["flanks_json"]))
        doc_genes = doc["genes"]
        genes, cov, perm, fam_loci = load_sample_truth(
            m["genes_fa"], m["genomes_fa"], m["genomes_tsv"], m["cov_col"], m["blast_tsv"], D=doc["flank_bp"]
        )

        # Build predicted pairs for this sample from cohort.json (adding prox seq)
        c_preds = {}
        for gname, plist in cohort_doc.get("pairs", {}).get(sid, {}).items():
            gentry = doc_genes.get(gname, {})
            l_tr = {t["trunk_id"]: t for t in extract_trunks(gentry.get("left", {}).get("paths", []), "left", W=a.pair_bp)}
            r_tr = {t["trunk_id"]: t for t in extract_trunks(gentry.get("right", {}).get("paths", []), "right", W=a.pair_bp)}
            out_l = []
            for pr in plist:
                lt = l_tr.get(pr["left_trunk"])
                rt = r_tr.get(pr["right_trunk"])
                if lt and rt:
                    item = dict(pr)
                    item["left_prox"] = lt["prox"]
                    item["right_prox"] = rt["prox"]
                    out_l.append(item)
            if out_l:
                c_preds[gname] = out_l

        ps = score_pairing_mode(doc_genes, genes, fam_loci, c_preds, a.pair_bp, 5.0, perm)
        tot_chk += ps["checkable"]
        tot_bad += ps["chimeras"]
        tot_ok += ps["confirmed"]
        tot_unc += ps["uncheckable"]
        tot_true += ps["true_all"]
        tot_hit += ps["true_hit"]
        tot_multi += ps["true_multi"]
        tot_multi_hit += ps["true_multi_hit"]
        tot_fam_chk += ps["fam_checkable"]
        tot_fam_bad += ps["fam_chimeras"]
        for tname, td in ps["by_tier"].items():
            by_tier[tname][0] += td["bad"]
            by_tier[tname][1] += td["chk"]
            if tname.startswith("cohort_"):
                cohort_only_bad += td["bad"]
                cohort_only_chk += td["chk"]
        for bd in ps["bad_details"]:
            bad_pairs.append(dict(bd, sample=sid))

        # Score Part B (post-repeat exit arms) strictly against exact source genomes!
        s_reps = cohort_doc.get("repeats", {}).get(sid, {})
        if s_reps:
            with tempfile.TemporaryDirectory(prefix="argj_crep_") as td:
                rfa = os.path.join(td, "rep.fa")
                qfa = os.path.join(td, "tru.fa")
                paf = os.path.join(td, "aln.paf")
                rinfo = {}
                with open(rfa, "w") as f:
                    for gname, smap in s_reps.items():
                        gseq_s = genes.get(gname, "")
                        pad = min(PADG, len(gseq_s))
                        for side, exts in smap.items():
                            for ex in exts:
                                flk = ex["flank"]
                                cl = min(int(ex["cohort_certified_bp"]), len(flk))
                                rid = f"{sid}|{gname}|{side}|{ex['cand_id']}"
                                ref_seq = (flk[-cl:] + gseq_s[:pad]) if side == "left" else (gseq_s[-pad:] + flk[:cl])
                                f.write(f">{rid}\n{ref_seq}\n")
                                rinfo[rid] = (gname, family(gname), side, cl, pad, ex)
                qmeta = {}
                with open(qfa, "w") as f:
                    qi = 0
                    for fam, flist in fam_loci.items():
                        for loc in flist:
                            qi += 1
                            for side in ("left", "right"):
                                pseq, plen = loc["padded"][side]
                                seq_5to3 = pseq[::-1] if side == "left" else pseq
                                qid = f"T{qi}|{fam}|{side}|{loc['genome']}"
                                f.write(f">{qid}\n{seq_5to3}\n")
                                qmeta[qid] = (fam, side, len(pseq) - plen)
                alns = collections.defaultdict(list)
                if rinfo and qmeta:
                    subprocess.run(
                        [a.minimap2, "-c", "-x", "map-hifi", "--eqx", "-N", "2000", "-p", "0.001", "-t", "8", rfa, qfa],
                        stdout=open(paf, "w"), stderr=subprocess.DEVNULL, check=True,
                    )
                    for line in open(paf):
                        p = line.rstrip("\n").split("\t")
                        qid, tname = p[0], p[5]
                        if tname not in rinfo or qid not in qmeta:
                            continue
                        _, rfam, rside, cl, pad, _ = rinfo[tname]
                        qfam, qside, _ = qmeta[qid]
                        if rfam != qfam or rside != qside or p[4] != "+":
                            continue
                        qlen, qs, qe, rs, re_ = int(p[1]), int(p[2]), int(p[3]), int(p[7]), int(p[8])
                        cg = next((t[5:] for t in p[12:] if t.startswith("cg:Z:")), None)
                        if cg is None:
                            continue
                        diffs, longs, rpos = [], [], rs
                        right = rside == "right"
                        for n_s, op in CIG.findall(cg):
                            n = int(n_s)
                            if op in "=M":
                                rpos += n
                            elif op == "X":
                                for b in range(rpos, rpos + n):
                                    diffs.append(((b - pad) if right else (cl - 1 - b), 1))
                                rpos += n
                            elif op == "D":
                                o = (rpos - pad) if right else (cl - rpos - n)
                                diffs.append((o, n))
                                if n >= 50:
                                    longs.append((o, n))
                                rpos += n
                            elif op == "I":
                                o = (rpos - pad) if right else (cl - rpos)
                                diffs.append((o, n))
                                if n >= 50:
                                    longs.append((o, n))
                        if right:
                            alns[tname].append(Aln(rs - pad, re_ - pad, qlen, qs, qe, diffs, longs, qid))
                        else:
                            alns[tname].append(Aln(cl - re_, cl - rs, qlen, qs, qe, diffs, longs, qid))
                for rid, (gname, fam, side, cl, pad, ex) in rinfo.items():
                    ok, hq = check_claim_strict(alns.get(rid, []), cl, start_x=-GENE_IN, thr=50, div=0.03, sub=STRICT_SUB)
                    rep_total += 1
                    if ok:
                        rep_pass += 1
                    else:
                        rep_fail += 1
                    gain = cl - int(ex.get("single_certain_bp", 0))
                    rep_gains.append(gain)
                    rep_rows.append({
                        "sample": sid, "gene": gname, "side": side, "cand_id": ex["cand_id"],
                        "fork_bp": ex["fork_bp"], "single_bp": ex["single_certain_bp"],
                        "cohort_bp": cl, "gain_bp": gain, "w": ex["w"], "dChi2": ex["dChi2"],
                        "loo_dChi2": ex["loo_dChi2"], "passed_strict": int(ok), "hit": hq or "-",
                    })

    lo, hi = cp_ci(tot_bad, tot_chk)
    lo_c, hi_c = cp_ci(cohort_only_bad, cohort_only_chk)
    lo_r, hi_r = cp_ci(rep_fail, rep_total)
    res = {
        "label": a.label,
        "n_samples": len(manifest),
        "pair_bp": a.pair_bp,
        "partA_pairing": {
            "checkable": tot_chk, "confirmed": tot_ok, "chimeras": tot_bad, "uncheckable": tot_unc,
            "chimera_rate_pct": round(100.0 * tot_bad / tot_chk, 2) if tot_chk else 0.0,
            "ci95_pct": [round(100.0 * lo, 2), round(100.0 * hi, 2)],
            "fam_checkable": tot_fam_chk, "fam_chimeras": tot_fam_bad,
            "cohort_only_checkable": cohort_only_chk, "cohort_only_chimeras": cohort_only_bad,
            "cohort_only_ci95_pct": [round(100.0 * lo_c, 2), round(100.0 * hi_c, 2)],
            "true_all": tot_true, "true_hit": tot_hit,
            "recall_pct": round(100.0 * tot_hit / tot_true, 1) if tot_true else 0.0,
            "true_multi": tot_multi, "true_multi_hit": tot_multi_hit,
            "multi_recall_pct": round(100.0 * tot_multi_hit / tot_multi, 1) if tot_multi else 0.0,
            "by_tier": {k: {"bad": v[0], "chk": v[1]} for k, v in sorted(by_tier.items())},
            "bad_pairs": bad_pairs,
        },
        "partB_repeats": {
            "total": rep_total, "passed_strict": rep_pass, "failed_strict": rep_fail,
            "ci95_pct": [round(100.0 * lo_r, 2), round(100.0 * hi_r, 2)],
            "mean_gain_bp": round(sum(rep_gains) / len(rep_gains), 1) if rep_gains else 0.0,
            "rows": rep_rows,
        },
    }
    json.dump(res, open(a.out_prefix + ".cohort_score.json", "w"), indent=2)
    print(f"=== COHORT {a.label} (S={len(manifest)}, W={a.pair_bp}) ===")
    print(f"  Part A pairing : {tot_bad}/{tot_chk} chimeras (cohort-added: {cohort_only_bad}/{cohort_only_chk}) | recall {tot_hit}/{tot_true} ({res['partA_pairing']['recall_pct']}%), multi {tot_multi_hit}/{tot_multi} ({res['partA_pairing']['multi_recall_pct']}%)")
    print(f"  By tier        : {res['partA_pairing']['by_tier']}")
    print(f"  Part B repeats : {rep_fail}/{rep_total} failed strict (mean extension gain +{res['partB_repeats']['mean_gain_bp']} bp)")


if __name__ == "__main__":
    main()
