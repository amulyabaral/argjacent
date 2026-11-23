#!/usr/bin/env bash
#SBATCH --job-name=argjacent_cohort
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/argjacent_cohort_%j.out
#SBATCH --error=logs/argjacent_cohort_%j.err

# Slurm post-array job:
#   1. Builds the cohort sample sheet from all completed samples in $OUT_DIR/samples/
#   2. Runs `argjacent-cohort` across all completed samples (when >= 2 samples have genes)
#   3. Runs `hpc/summarize_run.py` to merge all per-sample tables and produce run_metrics.tsv

set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ -n "${ENV_SETUP_CMD:-}" ]]; then
    eval "$ENV_SETUP_CMD"
fi

SAMPLES_TSV="${SAMPLES_TSV:?Set SAMPLES_TSV to your sample sheet}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to your output directory}"
GENES_FA="${GENES_FA:-}"
COHORT_PAIR_BP="${COHORT_PAIR_BP:-500}"
COHORT_EXTRA_ARGS="${COHORT_EXTRA_ARGS:-}"
WRITE_COHORT_GFA="${WRITE_COHORT_GFA:-0}"
SKIP_COHORT="${SKIP_COHORT:-0}"

SCRIPT_DIR="${HPC_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY_BIN="${PYTHON_BIN:-python3}"

# Locate argjacent-cohort executable
if [[ -n "${ARGJACENT_COHORT_BIN:-}" ]]; then
    COHORT_BIN="$ARGJACENT_COHORT_BIN"
elif command -v argjacent-cohort >/dev/null 2>&1; then
    COHORT_BIN="$(command -v argjacent-cohort)"
elif [[ -n "${REPO_ROOT:-}" && -x "${REPO_ROOT}/bin/argjacent-cohort" ]]; then
    COHORT_BIN="${REPO_ROOT}/bin/argjacent-cohort"
elif [[ -x "${SCRIPT_DIR}/../bin/argjacent-cohort" ]]; then
    COHORT_BIN="${SCRIPT_DIR}/../bin/argjacent-cohort"
else
    echo "ERROR: argjacent-cohort executable not found." >&2
    exit 1
fi

COHORT_DIR="${OUT_DIR}/cohort"
SUMMARY_DIR="${OUT_DIR}/summary"
mkdir -p "$COHORT_DIR" "$SUMMARY_DIR"

echo "=== ARGjacent HPC cohort & summary job (host=$(hostname)) ==="

# 1. Build cohort sample sheet from completed samples
COHORT_SHEET="${COHORT_DIR}/cohort_samplesheet.tsv"
N_READY="$("$PY_BIN" - "$SAMPLES_TSV" "$OUT_DIR" "$COHORT_SHEET" "$GENES_FA" <<'PY'
import json, os, sys
samples_tsv, out_dir, out_sheet, shared_genes = sys.argv[1:5]
rows = []
with open(samples_tsv) as fh:
    hdr = [c.strip() for c in fh.readline().rstrip("\r\n").split("\t")]
    col = {c: i for i, c in enumerate(hdr)}
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\r\n").split("\t")
        sid = parts[col["sample"]].strip()
        sdir = os.path.join(out_dir, "samples", sid)
        jpath = os.path.join(sdir, f"{sid}.flanks.json")
        if not (os.path.isfile(jpath) and os.path.getsize(jpath) > 0):
            continue
        try:
            doc = json.load(open(jpath))
        except Exception:
            continue
        if not doc.get("genes"):
            continue
        gpath = ""
        gfile = os.path.join(sdir, f"{sid}.genes_path.txt")
        if os.path.isfile(gfile):
            cand = open(gfile).read().strip()
            if cand and os.path.isfile(cand):
                gpath = cand
        if not gpath and "genes_fa" in col and len(parts) > col["genes_fa"] and parts[col["genes_fa"]].strip():
            gpath = parts[col["genes_fa"]].strip()
        if not gpath:
            gpath = shared_genes
        rtsv = ""
        for cand in (os.path.join(sdir, f"{sid}.tmp", "aw.reads.tsv.gz"),
                     os.path.join(sdir, f"{sid}.tmp", "aw.reads.tsv")):
            if os.path.isfile(cand):
                rtsv = cand
                break
        rows.append((sid, os.path.abspath(jpath), os.path.abspath(gpath) if gpath else "", os.path.abspath(rtsv) if rtsv else ""))

with open(out_sheet, "w") as out:
    out.write("sample\tflanks_json\tgenes_fa\treads_tsv\n")
    for r in rows:
        out.write("\t".join(r) + "\n")
print(len(rows))
PY
)"

echo "  Completed samples with non-empty genes: ${N_READY}"

# 2. Run argjacent-cohort if requested and >= 2 samples are ready
if [[ "$SKIP_COHORT" != "1" && "$N_READY" -ge 2 ]]; then
    CMD=("$COHORT_BIN" --sample-sheet "$COHORT_SHEET" -o "${COHORT_DIR}/cohort" --pair-bp "$COHORT_PAIR_BP")
    if [[ -n "$GENES_FA" ]]; then
        CMD+=(-g "$GENES_FA")
    fi
    if [[ "$WRITE_COHORT_GFA" == "1" ]]; then
        CMD+=(--gfa-dir "${COHORT_DIR}/cohort_gfa")
    fi
    if [[ -n "$COHORT_EXTRA_ARGS" ]]; then
        # shellcheck disable=SC2206
        EXTRA_ARR=($COHORT_EXTRA_ARGS)
        CMD+=("${EXTRA_ARR[@]}")
    fi
    echo "  Running: ${CMD[*]}"
    START_TS=$(date +%s)
    "${CMD[@]}"
    echo "  argjacent-cohort finished in $(( $(date +%s) - START_TS )) s."
elif [[ "$SKIP_COHORT" == "1" ]]; then
    echo "  SKIP_COHORT=1 set; skipping argjacent-cohort."
else
    echo "  Fewer than 2 completed samples with genes (${N_READY}); skipping argjacent-cohort."
fi

# 3. Always run summary merger & QC report
"$PY_BIN" "${SCRIPT_DIR}/summarize_run.py" \
    --samples "$SAMPLES_TSV" \
    --out-dir "$OUT_DIR" \
    --summary-dir "$SUMMARY_DIR"

echo "=== Cohort & summary job completed ==="
