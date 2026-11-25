# Running ARGjacent at Scale on a Slurm HPC

This folder gives you a ready-to-run Slurm workflow for running **ARGjacent** across dozens to hundreds of metagenomes, followed automatically by multi-sample **cohort pairing (`argjacent-cohort`)** and merged summary/QC tables.

---

## 1. One-Time Setup on Your HPC Login Node

Clone the repository and compile the three C programs (`aw`, `awx`, `awv`) on Linux:

```bash
# If your HPC uses environment modules, load a C compiler + Python 3.8+ (with numpy) + zlib first, e.g.:
# module load gcc python zlib

git clone https://github.com/amulyabaral/argjacent.git
cd argjacent

# Option A: install inside a conda/venv environment (recommended)
python3 -m pip install -U pip setuptools
python3 -m pip install .
argjacent-selftest          # ~10 s check; should end with "All checks passed."

# Option B: self-contained folder install (no pip install needed outside this folder)
sh install.sh
sh demo/run_demo.sh
```

> **Tip:** If your login node and compute nodes have different CPU architectures, run `pip install .` (or `sh install.sh`) inside an interactive compute session (`srun --cpus-per-task=2 --mem=4G --pty bash`).

---

## 2. Build & Validate Your Sample Sheet

Point `make_samplesheet.py` at the folder containing your paired `.fastq.gz` / `.fq.gz` reads (`*_R1/*_R2` or `*_1/*_2`):

```bash
python3 hpc/make_samplesheet.py \
  --fastq-dir /path/to/trimmed_reads \
  -o samples.tsv
```

- Use `--recursive` if each sample sits in its own subfolder.
- Use `--genes-dir /path/to/per_sample_genes` if each sample already has its own `<sample>.fa` gene list (e.g., from KMA or ResFinder).
- Or validate a hand-made TSV (columns `sample`, `r1`, `r2`, optional `genes_fa`):
  ```bash
  python3 hpc/make_samplesheet.py --check samples.tsv
  ```

---

## 3. Run a Fast 1-Sample Pilot First

Always preview with `--dry-run`, then test **1 sample** (`--limit 1`) on a real compute node before launching the full cohort:

```bash
# Preview exact Slurm commands without submitting
bash hpc/submit_pipeline.sh \
  --samples samples.tsv \
  --genes /path/to/resfinder.fa \
  --out /scratch/$USER/argjacent_pilot \
  --partition YOUR_PARTITION \
  --account YOUR_ACCOUNT \
  --limit 1 --dry-run

# Submit the 1-sample pilot
bash hpc/submit_pipeline.sh \
  --samples samples.tsv \
  --genes /path/to/resfinder.fa \
  --out /scratch/$USER/argjacent_pilot \
  --partition YOUR_PARTITION \
  --account YOUR_ACCOUNT \
  --cpus 8 --mem 24G --time 01:00:00 \
  --limit 1
```

If you prefer on-the-fly KMA pre-screening per sample instead of passing a full gene database directly to `argjacent`, replace `--genes` with `--kma-db /path/to/resfinder_kma_db` (and optionally `--kma-ref-fa /path/to/resfinder.fa`).

---

## 4. Launch the Full Large Run + Automatic Cohort Step

Once your 1-sample pilot succeeds and you check its runtime and peak RAM in `summary/run_metrics.tsv`, launch the full cohort:

```bash
bash hpc/submit_pipeline.sh \
  --samples samples.tsv \
  --genes /path/to/resfinder.fa \
  --out /scratch/$USER/argjacent_large_run \
  --partition YOUR_PARTITION \
  --account YOUR_ACCOUNT \
  --cpus 8 --mem 24G --time 02:00:00 \
  --max-concurrent 20
```

What happens automatically:
1. **Slurm array (`hpc/run_array.sbatch`)** runs up to 20 samples in parallel. Each task:
   - pins `OPENBLAS_NUM_THREADS=1` / `OMP_NUM_THREADS=1` so worker processes never oversubscribe a large node,
   - writes all temporary files (`*.tmp/`) onto fast node-local scratch (`$SLURM_TMPDIR` or `$TMPDIR`) rather than shared Lustre/NFS storage,
   - sets `--max-mem` to 80% of the Slurm memory allocation automatically,
   - copies only the compact final files back and writes a `.done` marker so re-runs skip finished samples.
2. **Dependent cohort & summary job (`hpc/run_cohort.sbatch`)** runs automatically as soon as the array finishes:
   - runs `argjacent-cohort` across all completed samples to resolve tied two-sided pairs (`cohort/cohort.cohort_pairs.tsv`),
   - merges every sample into `summary/all_samples.{summary,contexts,pairs}.tsv` and writes a per-sample QC/runtime/RAM table at `summary/run_metrics.tsv`.

---

## 5. Resuming Failed or Timed-Out Samples

If a few samples in a large run hit walltime or a node hiccup, bump `--time` or `--mem` and resubmit **only** the incomplete indices:

```bash
bash hpc/submit_pipeline.sh \
  --samples samples.tsv \
  --genes /path/to/resfinder.fa \
  --out /scratch/$USER/argjacent_large_run \
  --partition YOUR_PARTITION --account YOUR_ACCOUNT \
  --cpus 8 --mem 48G --time 04:00:00 \
  --resubmit-missing
```

---

## 6. Resource Sizing Guide

| Sample depth (2×150 bp, trimmed `.fq.gz`) | Suggested `--cpus` | Suggested `--mem` | Suggested `--time` | Notes |
|---|---:|---:|---:|---|
| Small / isolate (`< 2 GB` R1+R2, `< 5M` pairs) | 4–8 | `8G` | `00:30:00` | Usually finishes in < 1 min |
| Standard gut metagenome (`2–8 GB`, `15–50M` pairs) | 8 | `24G` | `01:30:00` | Typically 2–6 min, 4–10 GB RAM |
| Deep metagenome (`8–20 GB`, `50–120M` pairs) | 8–16 | `32G–48G` | `03:00:00` | If memory-tight, add `--extra-args "-w 40"` (~40% less RAM) |
| Cohort step (`50–500` samples) | 4 | `16G–32G` | `02:00:00` | Runs directly from `flanks.json` (no raw reads needed) |

---

## 7. Output Layout

```text
OUT_DIR/
├── samples/
│   └── <sample>/
│       ├── <sample>.summary.tsv      # one row per gene side
│       ├── <sample>.contexts.tsv     # all candidate contexts with certainty lengths
│       ├── <sample>.flanks.fa        # candidate flanking DNA sequences
│       ├── <sample>.flanks.json      # full machine-readable record
│       ├── <sample>.pairs.tsv        # certified single-sample two-sided pairs
│       ├── <sample>.gfa              # lossless compacted GFA v1.0 graph (Bandage-ready)
│       ├── <sample>.log              # detailed run log
│       └── <sample>.done             # completion marker
├── cohort/
│   ├── cohort.cohort_pairs.tsv       # single-sample + cohort-resolved two-sided pairs (Part A)
│   ├── cohort.cohort_repeats.tsv     # post-repeat exits (Part B: exploratory in chained repeats)
│   └── cohort.cohort.json            # full cohort JSON output
├── summary/
│   ├── run_metrics.tsv               # per-sample QC, read pairs, runtime, peak RAM, flank counts
│   ├── all_samples.summary.tsv       # merged gene-side summary across all samples
│   ├── all_samples.contexts.tsv      # merged candidate contexts across all samples
│   └── all_samples.pairs.tsv         # merged single-sample two-sided pairs across all samples
└── logs/
    ├── argjacent_<jobid>_<idx>.out
    └── argjacent_cohort_<jobid>.out
```

> **Note on `cohort.cohort_repeats.tsv` (Part B):** Cohort two-sided pairing (`cohort_pairs.tsv`, Part A) achieved **0% added chimeras** across all multi-sample benchmarks (CAMI II Strain Madness, CAMI III Human Gut, and real infant gut HiFi cohorts). Post-repeat exit jumping (`cohort_repeats.tsv`, Part B) is accurate for isolated repeats (e.g., single `IS26`) but should be treated as **exploratory** in dense clinical communities with chained mobile elements.
