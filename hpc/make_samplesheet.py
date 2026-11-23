#!/usr/bin/env python3
"""Scan FASTQ directories and build a validated sample sheet for ARGjacent HPC runs.

Output TSV columns:
    sample   r1   r2   [genes_fa]

Examples:
    python3 hpc/make_samplesheet.py --fastq-dir /data/reads -o samples.tsv
    python3 hpc/make_samplesheet.py --fastq-dir /data/reads --recursive --genes-dir /data/kma_genes -o samples.tsv
    python3 hpc/make_samplesheet.py --check samples.tsv
"""
import argparse
import os
import re
import sys

R1_PATTERNS = [
    (re.compile(r"^(?P<id>.+?)(?:_S\d+)?(?:_L\d+)?_R1(?:_001)?(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$"), "_R1", "_R2"),
    (re.compile(r"^(?P<id>.+?)\.R1(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$"), ".R1", ".R2"),
    (re.compile(r"^(?P<id>.+?)_1(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$"), "_1", "_2"),
    (re.compile(r"^(?P<id>.+?)\.1(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$"), ".1", ".2"),
]

SAFE_ASCII = re.compile(r"^[!-~]+$")


def clean_sample_id(raw):
    sid = re.sub(r"\s+", "_", raw.strip())
    sid = re.sub(r"[^A-Za-z0-9._-]", "_", sid)
    return sid.strip("._-") or "sample"


def discover_pairs(dirs, recursive=False):
    found = []
    seen_r1 = set()
    for root_dir in dirs:
        root_abs = os.path.abspath(root_dir)
        if not os.path.isdir(root_abs):
            sys.exit(f"make_samplesheet: directory not found: {root_dir}")
        walker = os.walk(root_abs) if recursive else [(root_abs, [], sorted(os.listdir(root_abs)))]
        for dp, dns, fns in walker:
            dns.sort()
            for fn in sorted(fns):
                full = os.path.join(dp, fn)
                if not os.path.isfile(full) or full in seen_r1:
                    continue
                for rx, tok1, tok2 in R1_PATTERNS:
                    m = rx.match(fn)
                    if not m:
                        continue
                    # Replace the last occurrence of tok1 before the extension with tok2
                    ext = m.group("ext")
                    stem = fn[:-len(ext)]
                    idx = stem.rfind(tok1)
                    if idx < 0:
                        continue
                    fn2 = stem[:idx] + tok2 + stem[idx + len(tok1):] + ext
                    full2 = os.path.join(dp, fn2)
                    if os.path.isfile(full2):
                        sid = clean_sample_id(m.group("id"))
                        found.append((sid, full, full2))
                        seen_r1.add(full)
                        break
    return found


def find_gene_fa(sid, genes_dir):
    if not genes_dir:
        return ""
    for cand in (
        f"{sid}.fa",
        f"{sid}.fna",
        f"{sid}.fasta",
        f"{sid}.genes.fa",
        f"{sid}_genes.fa",
        os.path.join(sid, "genes.fa"),
        os.path.join(sid, f"{sid}.genes.fa"),
    ):
        p = os.path.join(genes_dir, cand)
        if os.path.isfile(p) and os.path.getsize(p) > 0:
            return os.path.abspath(p)
    return ""


def validate_sheet(path):
    if not os.path.isfile(path):
        sys.exit(f"make_samplesheet: file not found: {path}")
    rows = []
    seen = set()
    errors = []
    with open(path) as fh:
        header_line = fh.readline()
        if not header_line:
            sys.exit(f"make_samplesheet: empty file: {path}")
        cols = [c.strip() for c in header_line.rstrip("\r\n").split("\t")]
        col = {c: i for i, c in enumerate(cols)}
        for req in ("sample", "r1", "r2"):
            if req not in col:
                sys.exit(f"make_samplesheet: missing required column '{req}' in {path} (found: {cols})")
        for lno, line in enumerate(fh, 2):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 3:
                errors.append(f"line {lno}: fewer than 3 tab-separated columns")
                continue
            sid = parts[col["sample"]].strip()
            r1 = parts[col["r1"]].strip()
            r2 = parts[col["r2"]].strip()
            gfa = parts[col["genes_fa"]].strip() if "genes_fa" in col and len(parts) > col["genes_fa"] else ""
            if not sid or not SAFE_ASCII.match(sid):
                errors.append(f"line {lno}: invalid sample ID '{sid}' (must be non-empty printable ASCII without spaces)")
            if sid in seen:
                errors.append(f"line {lno}: duplicate sample ID '{sid}'")
            seen.add(sid)
            for tag, p in (("r1", r1), ("r2", r2)):
                if not os.path.isfile(p):
                    errors.append(f"line {lno} ({sid}): {tag} not found: {p}")
                elif os.path.getsize(p) == 0:
                    errors.append(f"line {lno} ({sid}): {tag} is empty: {p}")
            if os.path.abspath(r1) == os.path.abspath(r2):
                errors.append(f"line {lno} ({sid}): r1 and r2 point to the same file")
            if gfa and not os.path.isfile(gfa):
                errors.append(f"line {lno} ({sid}): genes_fa not found: {gfa}")
            rows.append((sid, r1, r2, gfa))

    if errors:
        sys.stderr.write("Sample sheet validation FAILED:\n")
        for e in errors[:30]:
            sys.stderr.write(f"  - {e}\n")
        if len(errors) > 30:
            sys.stderr.write(f"  ... and {len(errors) - 30} more errors\n")
        sys.exit(1)

    if not rows:
        sys.exit(f"make_samplesheet: no sample rows found in {path}")

    sizes_gb = [(os.path.getsize(r1) + os.path.getsize(r2)) / (1024 ** 3) for _, r1, r2, _ in rows]
    tot_gb = sum(sizes_gb)
    max_gb = max(sizes_gb)
    with_g = sum(1 for _, _, _, g in rows if g)
    rec_mem = 16 if max_gb < 3.0 else (24 if max_gb < 8.0 else (32 if max_gb < 15.0 else 48))
    print(f"OK: {len(rows)} valid samples in {path}")
    print(f"  Total FASTQ size : {tot_gb:.2f} GB (median {sorted(sizes_gb)[len(sizes_gb)//2]:.2f} GB, max {max_gb:.2f} GB/sample)")
    if with_g:
        print(f"  Per-sample genes : {with_g}/{len(rows)} samples have genes_fa")
    print(f"  Suggested Slurm  : --cpus 8 --mem {rec_mem}G --time 02:00:00 (add -w 40 for ultra-deep >15 GB samples)")


def main():
    ap = argparse.ArgumentParser(description="Create or validate a paired-FASTQ sample sheet for ARGjacent Slurm runs.")
    ap.add_argument("--fastq-dir", nargs="+", help="One or more directories containing paired FASTQ files")
    ap.add_argument("--recursive", action="store_true", help="Search subdirectories recursively")
    ap.add_argument("--genes-dir", help="Optional directory containing per-sample FASTA files (<sample>.fa / <sample>.genes.fa)")
    ap.add_argument("--genes", help="Optional single shared genes.fa written into the genes_fa column for every sample")
    ap.add_argument("-o", "--out", default="samples.tsv", help="Output sample sheet TSV path (default: samples.tsv)")
    ap.add_argument("--check", metavar="TSV", help="Validate an existing sample sheet and exit")
    args = ap.parse_args()

    if args.check:
        validate_sheet(args.check)
        return

    if not args.fastq_dir:
        ap.error("provide --fastq-dir DIR [DIR ...] to build a sample sheet, or --check samples.tsv to validate one")

    pairs = discover_pairs(args.fastq_dir, recursive=args.recursive)
    if not pairs:
        sys.exit("make_samplesheet: no paired FASTQ files found matching *_R1/*_R2 or *_1/*_2")

    # Disambiguate duplicate sample IDs cleanly (e.g. from multiple subdirectories)
    counts = {}
    rows = []
    shared_g = os.path.abspath(args.genes) if args.genes else ""
    if shared_g and not os.path.isfile(shared_g):
        sys.exit(f"make_samplesheet: --genes file not found: {args.genes}")

    include_gcol = bool(args.genes_dir or shared_g)
    for sid, r1, r2 in sorted(pairs, key=lambda x: (x[0], x[1])):
        counts[sid] = counts.get(sid, 0) + 1
        final_sid = sid if counts[sid] == 1 else f"{sid}_{counts[sid]}"
        gfa = find_gene_fa(final_sid, args.genes_dir) if args.genes_dir else shared_g
        rows.append((final_sid, r1, r2, gfa))

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w") as fh:
        if include_gcol:
            fh.write("sample\tr1\tr2\tgenes_fa\n")
            for sid, r1, r2, gfa in rows:
                fh.write(f"{sid}\t{r1}\t{r2}\t{gfa}\n")
        else:
            fh.write("sample\tr1\tr2\n")
            for sid, r1, r2, _ in rows:
                fh.write(f"{sid}\t{r1}\t{r2}\n")

    validate_sheet(args.out)


if __name__ == "__main__":
    main()
