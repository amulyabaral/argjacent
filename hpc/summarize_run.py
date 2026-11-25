#!/usr/bin/env python3
"""Merge per-sample ARGjacent HPC outputs into cohort-wide tables and a QC metrics report.

Outputs written into --summary-dir (default: <out-dir>/summary):
  - run_metrics.tsv          : one row per sample (status, read pairs, runtime, peak RAM, genes, flanks, pairs)
  - all_samples.summary.tsv  : concatenated per-sample .summary.tsv with leading 'sample' column
  - all_samples.contexts.tsv : concatenated per-sample .contexts.tsv with leading 'sample' column
  - all_samples.pairs.tsv    : concatenated per-sample .pairs.tsv with leading 'sample' column
  - missing_indices.txt      : 1-based Slurm array indices of any incomplete samples (if any)
"""
import argparse
import os
import re
import statistics
import sys


RE_READS = re.compile(r"(\d+)\s+read pairs;\s+kept\s+(\d+)\s+near the genes.*?peak memory\s+([\d.]+)\s+MB")
RE_DONE = re.compile(r"Done in\s+([\d.]+)\s+s\.")
RE_LOG_PROP = re.compile(r"^proposing:\s+([\d.]+)\s+s", re.M)
RE_LOG_JUDG = re.compile(r"^judging:\s+([\d.]+)\s+s", re.M)
RE_LOG_PAIR = re.compile(r"^pairing:\s+(\d+)\s+certified", re.M)
RE_LOG_GFA = re.compile(r'"paths":\s*(\d+)')
RE_AW_PEAK = re.compile(r"peak_rss_mb=(\d+)")


def parse_log(path):
    info = {"wall_s": "", "peak_mb": "", "read_pairs": "", "kept_pairs": ""}
    if not os.path.isfile(path):
        return info
    txt = open(path, errors="replace").read()
    m = RE_DONE.search(txt)
    if m:
        info["wall_s"] = m.group(1)
    else:
        # Estimate total from aw + proposing + judging if Done line wasn't copied into .log
        mp = RE_LOG_PROP.search(txt)
        mj = RE_LOG_JUDG.search(txt)
        m_aw = re.search(r"pass2:.*?\(([0-9.]+)s\)", txt)
        m_p1 = re.search(r"pass1:.*?\(([0-9.]+)s\)", txt)
        tot = 0.0
        for mm in (m_p1, m_aw, mp, mj):
            if mm:
                tot += float(mm.group(1))
        if tot > 0:
            info["wall_s"] = f"{tot:.1f}"
    m2 = RE_READS.search(txt)
    if m2:
        info["read_pairs"], info["kept_pairs"], info["peak_mb"] = m2.group(1), m2.group(2), m2.group(3)
    else:
        m_rss = re.findall(r"rss(?:_mb)?[=:]\s*([\d.]+)", txt, re.I)
        if m_rss:
            info["peak_mb"] = str(max(float(x) for x in m_rss))
    return info


def main():
    ap = argparse.ArgumentParser(description="Merge and QC-summarize an ARGjacent HPC run.")
    ap.add_argument("--samples", required=True, help="Original samples.tsv submitted to the array")
    ap.add_argument("--out-dir", required=True, help="Run output directory containing samples/")
    ap.add_argument("--summary-dir", help="Where to write merged tables (default: <out-dir>/summary)")
    args = ap.parse_args()

    summary_dir = args.summary_dir or os.path.join(args.out_dir, "summary")
    os.makedirs(summary_dir, exist_ok=True)

    sample_ids = []
    with open(args.samples) as fh:
        hdr = [c.strip() for c in fh.readline().rstrip("\r\n").split("\t")]
        scol = hdr.index("sample")
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            sample_ids.append(line.rstrip("\r\n").split("\t")[scol].strip())

    sum_out = open(os.path.join(summary_dir, "all_samples.summary.tsv"), "w")
    ctx_out = open(os.path.join(summary_dir, "all_samples.contexts.tsv"), "w")
    prs_out = open(os.path.join(summary_dir, "all_samples.pairs.tsv"), "w")
    wrote_hdr = {"summary": False, "contexts": False, "pairs": False}

    metrics = []
    missing_idx = []
    tot_c1k = tot_c5k = tot_pairs = tot_sides_found = 0
    wall_vals = []
    ram_vals = []

    for idx, sid in enumerate(sample_ids, 1):
        sdir = os.path.join(args.out_dir, "samples", sid)
        done = os.path.isfile(os.path.join(sdir, f"{sid}.done"))
        p_sum = os.path.join(sdir, f"{sid}.summary.tsv")
        p_ctx = os.path.join(sdir, f"{sid}.contexts.tsv")
        p_prs = os.path.join(sdir, f"{sid}.pairs.tsv")
        p_jsn = os.path.join(sdir, f"{sid}.flanks.json")
        p_log = os.path.join(sdir, f"{sid}.log")

        if not (done and os.path.isfile(p_sum) and os.path.isfile(p_jsn)):
            missing_idx.append(idx)
            metrics.append({
                "array_idx": idx, "sample": sid, "status": "MISSING",
                "read_pairs": "", "kept_pairs": "", "wall_s": "", "peak_mb": "",
                "genes_found": 0, "sides_found": 0, "flanks_1kb": 0, "flanks_5kb": 0,
                "single_pairs": 0,
            })
            continue

        # Parse JSON aw_stats (authoritative read_pairs / kept_pairs / peak_rss)
        import json
        try:
            doc = json.load(open(p_jsn))
            st = doc.get("aw_stats", {}) or {}
        except Exception:
            st = {}
        log_i = parse_log(p_log)
        rp = str(int(st["reads"]) // 2) if st.get("reads", "").isdigit() else log_i["read_pairs"]
        kp = st.get("recruited_pairs", log_i["kept_pairs"])
        pm = st.get("peak_rss_mb_after_pass1" if st.get("map_from_file") == "1" else "peak_rss_mb", log_i["peak_mb"])
        ws = log_i["wall_s"]
        if ws:
            try:
                wall_vals.append(float(ws))
            except ValueError:
                pass
        if pm:
            try:
                ram_vals.append(float(pm))
            except ValueError:
                pass

        # Concatenate summary.tsv
        genes_found = set()
        sides_found = 0
        with open(p_sum) as fh:
            h = fh.readline()
            if h and not wrote_hdr["summary"]:
                sum_out.write("sample\t" + h)
                wrote_hdr["summary"] = True
            for line in fh:
                if not line.strip():
                    continue
                sum_out.write(f"{sid}\t{line}")
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 3 and parts[2] not in ("not_found", "failed"):
                    genes_found.add(parts[0])
                    sides_found += 1

        # Concatenate contexts.tsv
        c1k = c5k = 0
        if os.path.isfile(p_ctx):
            with open(p_ctx) as fh:
                h = fh.readline()
                if h and not wrote_hdr["contexts"]:
                    ctx_out.write("sample\t" + h)
                    wrote_hdr["contexts"] = True
                for line in fh:
                    if not line.strip():
                        continue
                    ctx_out.write(f"{sid}\t{line}")
                    parts = line.rstrip("\r\n").split("\t")
                    if len(parts) >= 5 and parts[4].isdigit():
                        cb = int(parts[4])
                        if cb >= 1000:
                            c1k += 1
                        if cb >= 5000:
                            c5k += 1

        # Concatenate pairs.tsv
        n_pr = 0
        if os.path.isfile(p_prs):
            with open(p_prs) as fh:
                h = fh.readline()
                if h and not wrote_hdr["pairs"]:
                    prs_out.write("sample\t" + h)
                    wrote_hdr["pairs"] = True
                for line in fh:
                    if not line.strip():
                        continue
                    prs_out.write(f"{sid}\t{line}")
                    n_pr += 1

        tot_sides_found += sides_found
        tot_c1k += c1k
        tot_c5k += c5k
        tot_pairs += n_pr
        metrics.append({
            "array_idx": idx, "sample": sid, "status": "OK",
            "read_pairs": rp, "kept_pairs": kp, "wall_s": ws, "peak_mb": pm,
            "genes_found": len(genes_found), "sides_found": sides_found,
            "flanks_1kb": c1k, "flanks_5kb": c5k, "single_pairs": n_pr,
        })

    sum_out.close()
    ctx_out.close()
    prs_out.close()

    mpath = os.path.join(summary_dir, "run_metrics.tsv")
    cols = ["array_idx", "sample", "status", "read_pairs", "kept_pairs", "wall_s", "peak_mb",
            "genes_found", "sides_found", "flanks_1kb", "flanks_5kb", "single_pairs"]
    with open(mpath, "w") as out:
        out.write("\t".join(cols) + "\n")
        for m in metrics:
            out.write("\t".join(str(m[c]) for c in cols) + "\n")

    miss_path = os.path.join(summary_dir, "missing_indices.txt")
    if missing_idx:
        with open(miss_path, "w") as out:
            out.write(",".join(str(i) for i in missing_idx) + "\n")
    elif os.path.exists(miss_path):
        os.remove(miss_path)

    # Cohort results count (if cohort step ran)
    cp_tot = cp_added = cr_tot = 0
    cp_path = os.path.join(args.out_dir, "cohort", "cohort.cohort_pairs.tsv")
    cr_path = os.path.join(args.out_dir, "cohort", "cohort.cohort_repeats.tsv")
    if os.path.isfile(cp_path):
        with open(cp_path) as fh:
            fh.readline()
            for line in fh:
                if not line.strip():
                    continue
                cp_tot += 1
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) >= 10 and parts[9].startswith("cohort_"):
                    cp_added += 1
    if os.path.isfile(cr_path):
        with open(cr_path) as fh:
            fh.readline()
            cr_tot = sum(1 for line in fh if line.strip())

    n_ok = len(sample_ids) - len(missing_idx)
    print(f"\n=== ARGjacent HPC Run Summary ({args.out_dir}) ===")
    print(f"  Samples completed  : {n_ok} / {len(sample_ids)}")
    if missing_idx:
        print(f"  Missing/failed IDs : {','.join(str(i) for i in missing_idx[:20])}{'...' if len(missing_idx) > 20 else ''}")
        print(f"  Tip: rerun missing tasks with `bash hpc/submit_pipeline.sh ... --resubmit-missing`")
    if wall_vals:
        print(f"  Per-sample runtime : median {statistics.median(wall_vals):.1f} s, max {max(wall_vals):.1f} s")
    if ram_vals:
        print(f"  Per-sample peak RAM: median {statistics.median(ram_vals):.0f} MB, max {max(ram_vals):.0f} MB")
    print(f"  Gene sides found   : {tot_sides_found:,}")
    print(f"  Certified flanks   : {tot_c1k:,} (>=1 kb), {tot_c5k:,} (>=5 kb)")
    print(f"  Two-sided pairs    : {tot_pairs:,} single-sample" +
          (f"; {cp_tot:,} cohort total (+{cp_added:,} resolved by cohort Part A)" if os.path.isfile(cp_path) else ""))
    if os.path.isfile(cr_path):
        print(f"  Post-repeat exits  : {cr_tot:,} (Part B: exploratory in chained repeats)")
    print(f"  Summary tables     : {summary_dir}/")


if __name__ == "__main__":
    main()
