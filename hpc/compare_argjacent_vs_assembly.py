#!/usr/bin/env python3
"""Head-to-head comparison between ARGjacent certified flanks and pre-made metagenome assemblies
(MGnify/ENA ERZ or Logan contigs) — either for a single sample or aggregated across an HPC run.

IMPORTANT SCIENTIFIC FRAMING:
  Pre-made short-read assemblies (MEGAHIT, metaSPAdes, Logan) are NOT ground truth — they are the
  assembly-based baseline. This script quantifies:
    1. Concordance on overlapping certified regions (identity over shared certified span)
    2. Assembly breaks where ARGjacent certifies >= 1 kb (argj_only_1kb)
    3. Multi-context recovery where ARGjacent certifies >= 2 distinct >= 1 kb contexts while the
       assembly collapses to <= 1 contig (argj_multi_vs_asm_single)
    4. Repeat-cliff over-extension where ARGjacent halts certain_bp at a shared repeat
       ('another_locus_on_the_path') while the assembly keeps extending >= 500 bp unverified
       (asm_past_repeat_cliff)
    5. Ready-to-annotate matched FASTA exports strictly sliced to certain_bp (ARGjacent) vs
       extracted contig flanks (assembly), for identical downstream MGE/plasmid/host/vector runs.

Examples:
    # Single sample
    python3 hpc/compare_argjacent_vs_assembly.py \\
        --argj-json sample.flanks.json --asm-json sample.asm_flanks.json \\
        --sample ERR12345 -o sample.cmp

    # Whole global HPC run with manifest metadata
    python3 hpc/compare_argjacent_vs_assembly.py \\
        --manifest global_2023.tsv --out-dir results/global_2023 \\
        -o results/global_2023/comparison/global_cmp --export-fasta
"""
import argparse
import collections
import json
import os
import statistics
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(seq):
    return seq.translate(COMP)[::-1]


def outward_seq(flank_5to3, side):
    """Orient a flank outward away from the gene boundary (index 0 = 1 bp outside gene)."""
    return flank_5to3 if side == "right" else rc(flank_5to3)


def certified_slice_5to3(flank_5to3, side, cert_bp):
    """Slice a 5'->3' flank strictly to its certified bases."""
    if cert_bp <= 0 or not flank_5to3:
        return ""
    c = min(cert_bp, len(flank_5to3))
    return flank_5to3[:c] if side == "right" else flank_5to3[-c:]


def kmer_jaccard_and_identity(out_a, out_b, min_len=200):
    """Compare two outward-oriented flank sequences over their shared prefix length.

    Uses fast exact prefix check + 21-mer containment/agreement (tolerant to small indels)
    and returns (overlap_bp, best_identity_pct, exact_lcp_bp).
    """
    ov = min(len(out_a), len(out_b))
    if ov < min_len:
        return ov, None, 0
    a = out_a[:ov]
    b = out_b[:ov]
    # Exact LCP
    lcp = 0
    while lcp < ov and a[lcp] == b[lcp]:
        lcp += 1
    if lcp == ov:
        return ov, 100.0, lcp
    # Fast direct Hamming identity
    direct_id = 100.0 * sum(1 for x, y in zip(a, b) if x == y) / ov
    if direct_id >= 98.0:
        return ov, direct_id, lcp
    # Indel-tolerant 21-mer shared fraction (catches 1-5 bp indels that shift direct Hamming)
    K = 21
    ka = {a[i:i + K] for i in range(0, ov - K + 1, 3)}
    if not ka:
        return ov, direct_id, lcp
    hits = sum(1 for i in range(0, ov - K + 1, 3) if b[i:i + K] in ka)
    n_tot = max(1, (ov - K) // 3 + 1)
    # Approximate base identity from 21-mer survival: p_kmer ~ id^21 -> id ~ p_kmer^(1/21)
    frac = hits / n_tot
    kmer_est_id = 100.0 * (frac ** (1.0 / K)) if frac > 0 else 0.0
    return ov, max(direct_id, kmer_est_id), lcp


def distinct_context_count(paths, side, min_bp=1000, use_cert=True):
    """Count distinct >= min_bp contexts (deduplicating near-identical branches >= 95% 21-mer overlap)."""
    reps = []
    K = 21
    for p in paths:
        fl = p.get("flank", "")
        bp = int(p.get("context_len", 0)) if use_cert else len(fl)
        if bp < min_bp:
            continue
        out = outward_seq(fl, side)[:bp]
        km = {out[i:i + K] for i in range(0, min(len(out), 1500) - K + 1, 5)}
        if not km:
            continue
        dup = False
        for rkm in reps:
            inter = len(km & rkm)
            denom = max(1, min(len(km), len(rkm)))
            if inter / denom >= 0.85:
                dup = True
                break
        if not dup:
            reps.append(km)
    return len(reps)


def compare_one_sample(sample_id, argj_doc, asm_doc, meta=None, min_bp=1000, cliff_over_bp=500):
    """Compare one sample's ARGjacent flanks.json against its assembly asm_flanks.json."""
    meta = meta or {}
    ag = (argj_doc or {}).get("genes", {}) or {}
    mg = (asm_doc or {}).get("genes", {}) or {}
    all_genes = sorted(set(ag) | set(mg))

    rows = []
    argj_fasta = []
    asm_fasta = []

    biome = meta.get("biome_group", "")
    assembler = meta.get("assembler", "") or (asm_doc or {}).get("source", "")
    study = meta.get("study_accession", "")
    country = meta.get("country", "")

    for g in all_genes:
        for side in ("left", "right"):
            a_paths = ((ag.get(g) or {}).get(side) or {}).get("paths", []) or []
            m_paths = ((mg.get(g) or {}).get(side) or {}).get("paths", []) or []
            if not a_paths and not m_paths:
                continue

            # ARGjacent summary for this (gene, side)
            a_best_cert_p = max(a_paths, key=lambda p: (int(p.get("context_len", 0)), int(p.get("flank_len", 0))), default=None)
            a_cert_bp = int(a_best_cert_p.get("context_len", 0)) if a_best_cert_p else 0
            a_end_bp = max((int(p.get("end_len", 0)) for p in a_paths), default=0)
            a_long_bp = max((int(p.get("flank_len", 0)) for p in a_paths), default=0)
            a_cstop = (a_best_cert_p.get("certain_stop", "") if a_best_cert_p else "")
            a_depth = float(a_best_cert_p.get("depth", 0.0)) if a_best_cert_p else 0.0
            a_n_1kb = distinct_context_count(a_paths, side, min_bp=min_bp, use_cert=True)

            # Assembly summary for this (gene, side)
            m_best_p = max(m_paths, key=lambda p: int(p.get("flank_len", 0)), default=None)
            m_max_bp = int(m_best_p.get("flank_len", 0)) if m_best_p else 0
            m_stop = (m_best_p.get("stop", "") if m_best_p else "not_found")
            m_ctg = (m_best_p.get("contig", "") if m_best_p else "")
            m_n_1kb = distinct_context_count(m_paths, side, min_bp=min_bp, use_cert=False)

            # Export deduplicated certified ARGjacent flanks & assembly flanks >= min_bp for unified downstream annotation
            # (Multiple post-repeat grey candidate exits share the identical certified pre-repeat stem; keep 1 rep per distinct context!)
            seen_a_kms = []
            for p in sorted(a_paths, key=lambda x: (-int(x.get("context_len", 0)), -float(x.get("depth", 0.0)), x.get("id", ""))):
                cb = int(p.get("context_len", 0))
                if cb >= min_bp:
                    sl = certified_slice_5to3(p.get("flank", ""), side, cb)
                    if len(sl) >= min_bp:
                        out_s = outward_seq(sl, side)
                        km = {out_s[i:i + 21] for i in range(0, min(len(out_s), 2000) - 20, 5)}
                        if km and any(len(km & rkm) / max(1, min(len(km), len(rkm))) >= 0.90 for rkm in seen_a_kms):
                            continue
                        seen_a_kms.append(km)
                        cid = p.get("id", "c1")
                        argj_fasta.append((f"{sample_id}|{g}|{side}|{cid}|certain_bp={len(sl)}|biome={biome}", sl))
            seen_m_seqs = set()
            for p in sorted(m_paths, key=lambda x: (-len(x.get("flank", "")), x.get("id", ""))):
                fl = p.get("flank", "")
                if len(fl) >= min_bp and fl not in seen_m_seqs:
                    seen_m_seqs.add(fl)
                    cid = p.get("id", "a1")
                    asm_fasta.append((f"{sample_id}|{g}|{side}|{cid}|asm_bp={len(fl)}|biome={biome}", fl))

            # Best sequence match between any certified ARGjacent flank and any assembly flank
            best_ov = 0
            best_id = None
            best_lcp = 0
            for ap in a_paths:
                cb = int(ap.get("context_len", 0))
                if cb < 200:
                    continue
                a_out = outward_seq(ap.get("flank", ""), side)[:cb]
                for mp in m_paths:
                    m_fl = mp.get("flank", "")
                    if len(m_fl) < 200:
                        continue
                    m_out = outward_seq(m_fl, side)
                    ov, pid, lcp = kmer_jaccard_and_identity(a_out, m_out, min_len=200)
                    if pid is not None and (best_id is None or (ov >= 500 and best_ov < 500) or pid > best_id):
                        best_ov, best_id, best_lcp = ov, pid, lcp

            # Detect assembly extending past an ARGjacent repeat cliff
            hit_repeat_wall = (a_cstop == "another_locus_on_the_path")
            asm_past_cliff = (m_max_bp - a_cert_bp) if (hit_repeat_wall and m_max_bp >= a_cert_bp + cliff_over_bp) else 0

            # Primary diagnostic category (mutually exclusive, ordered by scientific importance)
            if a_n_1kb >= 2 and m_n_1kb <= 1:
                cat = "argj_multi_vs_asm_single"
            elif a_cert_bp >= min_bp and m_max_bp >= min_bp:
                if best_id is not None and best_id < 95.0:
                    cat = "certified_sequence_mismatch"
                elif asm_past_cliff >= cliff_over_bp:
                    cat = "asm_past_repeat_cliff"
                else:
                    cat = "both_concordant_1kb"
            elif a_cert_bp >= min_bp and m_max_bp < min_bp:
                cat = "argj_only_1kb"
            elif a_cert_bp < min_bp and m_max_bp >= min_bp:
                cat = "asm_past_repeat_cliff" if asm_past_cliff >= cliff_over_bp else "asm_only_1kb"
            else:
                cat = "both_short"

            rows.append({
                "sample": sample_id,
                "study_accession": study,
                "biome_group": biome,
                "country": country,
                "assembler": assembler,
                "gene": g,
                "side": side,
                "argj_certain_bp": a_cert_bp,
                "argj_end_bp": a_end_bp,
                "argj_longest_bp": a_long_bp,
                "argj_n_1kb": a_n_1kb,
                "argj_depth": f"{a_depth:.1f}",
                "argj_certain_stop": a_cstop,
                "asm_max_bp": m_max_bp,
                "asm_n_1kb": m_n_1kb,
                "asm_stop": m_stop,
                "asm_contig": m_ctg,
                "overlap_bp": best_ov,
                "overlap_id_pct": f"{best_id:.2f}" if best_id is not None else "",
                "exact_lcp_bp": best_lcp,
                "asm_past_cliff_bp": asm_past_cliff,
                "category": cat,
            })

    return rows, argj_fasta, asm_fasta


def load_manifest_meta(path):
    meta = {}
    with open(path) as fh:
        hdr = [c.strip() for c in fh.readline().rstrip("\r\n").split("\t")]
        col = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            row = {c: (parts[i].strip() if i < len(parts) else "") for c, i in col.items()}
            sid = row.get("sample")
            if sid:
                meta[sid] = row
    return meta


def gene_family(gene_name):
    """Collapse ResFinder allele suffixes (e.g. erm(B)_10_U86375 -> erm(B), tet(Q)_1_L33696 -> tet(Q))."""
    parts = gene_name.split("_")
    if len(parts) >= 3 and parts[1].isdigit():
        return parts[0]
    return gene_name


def summarize_group(label, rows):
    n = len(rows)
    cats = collections.Counter(r["category"] for r in rows)
    a_1kb = sum(1 for r in rows if r["argj_certain_bp"] >= 1000)
    m_1kb = sum(1 for r in rows if r["asm_max_bp"] >= 1000)
    a_5kb = sum(1 for r in rows if r["argj_certain_bp"] >= 5000)
    m_5kb = sum(1 for r in rows if r["asm_max_bp"] >= 5000)
    multi_a = sum(1 for r in rows if r["argj_n_1kb"] >= 2)
    multi_m = sum(1 for r in rows if r["asm_n_1kb"] >= 2)
    past_cliff = sum(1 for r in rows if r["asm_past_cliff_bp"] >= 500)
    ids = [float(r["overlap_id_pct"]) for r in rows if r["overlap_id_pct"] and r["overlap_bp"] >= 1000]
    concord_98 = sum(1 for v in ids if v >= 98.0)
    med_a = statistics.median([r["argj_certain_bp"] for r in rows]) if rows else 0
    med_m = statistics.median([r["asm_max_bp"] for r in rows]) if rows else 0

    # Family-deduplicated counts (1 count per unique (sample, gene_family, side))
    by_fam = collections.defaultdict(list)
    for r in rows:
        by_fam[(r["sample"], gene_family(r["gene"]), r["side"])].append(r)
    fam_a_1kb = sum(1 for rs in by_fam.values() if max(x["argj_certain_bp"] for x in rs) >= 1000)
    fam_m_1kb = sum(1 for rs in by_fam.values() if max(x["asm_max_bp"] for x in rs) >= 1000)
    fam_a_only = sum(1 for rs in by_fam.values() if max(x["argj_certain_bp"] for x in rs) >= 1000 and max(x["asm_max_bp"] for x in rs) < 1000)
    fam_m_only = sum(1 for rs in by_fam.values() if max(x["argj_certain_bp"] for x in rs) < 1000 and max(x["asm_max_bp"] for x in rs) >= 1000)
    fam_cliff = sum(1 for rs in by_fam.values() if max(x["asm_past_cliff_bp"] for x in rs) >= 500)

    return {
        "group": label,
        "samples": len({r["sample"] for r in rows}),
        "gene_sides": n,
        "family_sides": len(by_fam),
        "argj_1kb": a_1kb,
        "asm_1kb": m_1kb,
        "argj_1kb_families": fam_a_1kb,
        "asm_1kb_families": fam_m_1kb,
        "argj_only_1kb_families": fam_a_only,
        "asm_only_1kb_families": fam_m_only,
        "asm_past_cliff_families": fam_cliff,
        "argj_5kb": a_5kb,
        "asm_5kb": m_5kb,
        "argj_multi_ctx": multi_a,
        "asm_multi_ctx": multi_m,
        "argj_multi_vs_asm_single": cats["argj_multi_vs_asm_single"],
        "argj_only_1kb": cats["argj_only_1kb"],
        "asm_only_1kb": cats["asm_only_1kb"],
        "asm_past_repeat_cliff": past_cliff,
        "both_concordant_1kb": cats["both_concordant_1kb"],
        "certified_seq_mismatch": cats["certified_sequence_mismatch"],
        "overlap_1kb_pairs": len(ids),
        "overlap_id_ge98_pct": f"{100.0 * concord_98 / len(ids):.2f}" if ids else "",
        "median_argj_certain_bp": int(med_a),
        "median_asm_bp": int(med_m),
    }


def main():
    ap = argparse.ArgumentParser(description="Compare ARGjacent certified flanks against pre-made metagenome assemblies.")
    ap.add_argument("--argj-json", help="Single-sample ARGjacent .flanks.json")
    ap.add_argument("--asm-json", help="Single-sample extracted .asm_flanks.json")
    ap.add_argument("--sample", default="sample", help="Sample ID for single-sample mode")
    ap.add_argument("--manifest", help="Global manifest TSV (for cohort/global aggregation mode)")
    ap.add_argument("--out-dir", help="HPC output directory containing samples/<sid>/ (for global mode)")
    ap.add_argument("-o", "--out", required=True, help="Output prefix (writes _per_side.tsv, _summary.tsv, optional FASTAs)")
    ap.add_argument("--min-bp", type=int, default=1000, help="Primary comparison flank length threshold (default: 1000)")
    ap.add_argument("--cliff-over-bp", type=int, default=500,
                    help="Flag assembly over-extension past ARGjacent repeat cliff by >= N bp (default: 500)")
    ap.add_argument("--export-fasta", action="store_true",
                    help="Also write matched certified ARGjacent & assembly FASTAs for unified annotation")
    args = ap.parse_args()

    all_rows = []
    all_argj_fa = []
    all_asm_fa = []

    if args.argj_json and args.asm_json:
        ad = json.load(open(args.argj_json)) if os.path.isfile(args.argj_json) else {}
        md = json.load(open(args.asm_json)) if os.path.isfile(args.asm_json) else {}
        rows, afa, mfa = compare_one_sample(args.sample, ad, md, min_bp=args.min_bp, cliff_over_bp=args.cliff_over_bp)
        all_rows.extend(rows)
        all_argj_fa.extend(afa)
        all_asm_fa.extend(mfa)
    elif args.out_dir:
        meta_map = load_manifest_meta(args.manifest) if args.manifest else {}
        sroot = os.path.join(args.out_dir, "samples")
        sample_ids = sorted(meta_map) if meta_map else (sorted(os.listdir(sroot)) if os.path.isdir(sroot) else [])
        n_both = 0
        for sid in sample_ids:
            sdir = os.path.join(sroot, sid)
            aj = os.path.join(sdir, f"{sid}.flanks.json")
            mj = os.path.join(sdir, f"{sid}.asm_flanks.json")
            if not (os.path.isfile(aj) and os.path.isfile(mj)):
                continue
            try:
                ad = json.load(open(aj))
                md = json.load(open(mj))
            except Exception:
                continue
            n_both += 1
            rows, afa, mfa = compare_one_sample(
                sid, ad, md, meta=meta_map.get(sid, {}),
                min_bp=args.min_bp, cliff_over_bp=args.cliff_over_bp,
            )
            all_rows.extend(rows)
            if args.export_fasta:
                all_argj_fa.extend(afa)
                all_asm_fa.extend(mfa)
        sys.stderr.write(f"compare_argjacent_vs_assembly: paired {n_both} samples with both ARGjacent + assembly outputs.\n")
    else:
        ap.error("Provide either (--argj-json and --asm-json) or --out-dir.")

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    side_cols = [
        "sample", "study_accession", "biome_group", "country", "assembler",
        "gene", "side", "argj_certain_bp", "argj_end_bp", "argj_longest_bp",
        "argj_n_1kb", "argj_depth", "argj_certain_stop",
        "asm_max_bp", "asm_n_1kb", "asm_stop", "asm_contig",
        "overlap_bp", "overlap_id_pct", "exact_lcp_bp", "asm_past_cliff_bp", "category",
    ]
    with open(f"{args.out}_per_side.tsv", "w") as out:
        out.write("\t".join(side_cols) + "\n")
        for r in all_rows:
            out.write("\t".join(str(r[c]) for c in side_cols) + "\n")

    summaries = [summarize_group("ALL", all_rows)]
    by_biome = collections.defaultdict(list)
    for r in all_rows:
        if r["biome_group"]:
            by_biome[r["biome_group"]].append(r)
    for b in sorted(by_biome):
        summaries.append(summarize_group(f"biome:{b}", by_biome[b]))

    sum_cols = list(summaries[0].keys())
    with open(f"{args.out}_summary.tsv", "w") as out:
        out.write("\t".join(sum_cols) + "\n")
        for s in summaries:
            out.write("\t".join(str(s[c]) for c in sum_cols) + "\n")

    if args.export_fasta:
        with open(f"{args.out}_argjacent_certified.fa", "w") as out:
            for h, seq in all_argj_fa:
                out.write(f">{h}\n{seq}\n")
        with open(f"{args.out}_assembly_flanks.fa", "w") as out:
            for h, seq in all_asm_fa:
                out.write(f">{h}\n{seq}\n")

    tot = summaries[0]
    print(f"\n=== ARGjacent vs Pre-Made Assembly Comparison ({args.out}) ===")
    print(f"  Samples compared          : {tot['samples']:,}")
    print(f"  Gene sides evaluated      : {tot['gene_sides']:,}")
    print(f"  >= 1 kb recovered         : ARGjacent={tot['argj_1kb']:,} (certified) | Assembly={tot['asm_1kb']:,}")
    print(f"  >= 5 kb recovered         : ARGjacent={tot['argj_5kb']:,} (certified) | Assembly={tot['asm_5kb']:,}")
    print(f"  Multi-context (>=2 >=1kb) : ARGjacent={tot['argj_multi_ctx']:,} | Assembly={tot['asm_multi_ctx']:,} (ARGj multi vs Asm <=1: {tot['argj_multi_vs_asm_single']:,})")
    print(f"  ARGjacent-only >= 1 kb    : {tot['argj_only_1kb']:,} (assembly broke < 1 kb or missed gene)")
    print(f"  Asm past repeat cliff     : {tot['asm_past_repeat_cliff']:,} (assembly extended >={args.cliff_over_bp} bp past ARGj repeat wall)")
    print(f"  Certified overlap >= 1 kb : {tot['overlap_1kb_pairs']:,} pairs (identity >= 98%: {tot['overlap_id_ge98_pct'] or 'n/a'}%)")
    print(f"  Certified seq mismatches  : {tot['certified_seq_mismatch']:,} (<95% identity on >=1 kb overlap)")


if __name__ == "__main__":
    main()
