#!/usr/bin/env bash
# One-command Slurm launcher for large ARGjacent runs (array job + dependent cohort & summary job).
#
# Quick pilot (first 1 sample, dry-run first):
#   bash hpc/submit_pipeline.sh --samples samples.tsv --genes resfinder.fa --out results/pilot --limit 1 --dry-run
#
# Full large run (e.g. 200 samples, max 20 concurrent, auto cohort deconvolution):
#   bash hpc/submit_pipeline.sh --samples samples.tsv --genes resfinder.fa --out results/cohort \
#       --partition standard --account myproj --cpus 8 --mem 24G --time 02:00:00 --max-concurrent 20
#
# Resubmit only missing/failed samples after a partial run:
#   bash hpc/submit_pipeline.sh --samples samples.tsv --genes resfinder.fa --out results/cohort --resubmit-missing

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SAMPLES_TSV=""
GENES_FA=""
KMA_DB=""
KMA_REF_FA=""
OUT_DIR=""
PARTITION=""
ACCOUNT=""
QOS=""
CPUS=8
MEM="24G"
TIME="02:00:00"
COHORT_CPUS=4
COHORT_MEM="24G"
COHORT_TIME="02:00:00"
MAX_CONCURRENT=20
LIMIT=0
EXTRA_ARGS=""
COHORT_PAIR_BP=500
COHORT_EXTRA_ARGS=""
ENV_SETUP_CMD=""
NO_COHORT=0
WRITE_COHORT_GFA=0
KEEP_TMP=0
FORCE=0
RESUBMIT_MISSING=0
DRY_RUN=0

usage() {
    cat <<'EOF'
Usage: bash hpc/submit_pipeline.sh --samples samples.tsv --out OUT_DIR [options]

Required:
  --samples FILE          Tab-separated sample sheet (columns: sample, r1, r2 [, genes_fa])
  --out DIR               Output directory (writes samples/, cohort/, summary/, logs/)

Gene source (provide one, unless samples.tsv has a 'genes_fa' column):
  --genes FILE            Shared target genes FASTA (e.g. ResFinder/CARD/custom genes)
  --kma-db PREFIX         Optional: run per-sample KMA pre-screen against this KMA database
  --kma-ref-fa FILE       Optional FASTA matching --kma-db (uses reference alleles instead of KMA consensus)

Slurm resources:
  --partition NAME        Slurm partition (-p)
  --account NAME          Slurm account/project (-A)
  --qos NAME              Slurm QoS (--qos)
  --cpus N                CPUs per sample task (default: 8)
  --mem SIZE              RAM per sample task, e.g. 16G, 24G, 48G (default: 24G)
  --time HH:MM:SS         Walltime per sample task (default: 02:00:00)
  --max-concurrent N      Max simultaneous array tasks (default: 20)
  --cohort-cpus N         CPUs for the post-array cohort+summary job (default: 4)
  --cohort-mem SIZE       RAM for the post-array cohort+summary job (default: 24G)
  --cohort-time HH:MM:SS  Walltime for the post-array cohort+summary job (default: 02:00:00)

ARGjacent & workflow options:
  --limit N               Submit only the first N samples (great for 1-sample or 5-sample HPC pilot tests)
  --extra-args "ARGS"     Extra flags passed to `argjacent` (e.g. "-D 10000 -w 40 --pair-bp 500")
  --pair-bp BP            Flank radius for cohort pairing (default: 500)
  --write-cohort-gfa      Also emit cohort-updated per-sample GFA graphs
  --keep-tmp              Save gzipped aw.reads.tsv per sample (for repeat read-level mode)
  --no-cohort             Skip argjacent-cohort (still merges summary tables)
  --env-cmd "CMD"         Shell command run inside jobs first (e.g. "module load gcc python; source ~/venv/bin/activate")
  --force                 Re-run samples even if .done marker exists
  --resubmit-missing      Submit ONLY sample array indices without a completed .done marker
  --dry-run               Validate everything and print exact sbatch commands without submitting
  -h, --help              Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samples)          SAMPLES_TSV="$2"; shift 2 ;;
        --genes)            GENES_FA="$2"; shift 2 ;;
        --kma-db)           KMA_DB="$2"; shift 2 ;;
        --kma-ref-fa)       KMA_REF_FA="$2"; shift 2 ;;
        --out)              OUT_DIR="$2"; shift 2 ;;
        --partition)        PARTITION="$2"; shift 2 ;;
        --account)          ACCOUNT="$2"; shift 2 ;;
        --qos)              QOS="$2"; shift 2 ;;
        --cpus)             CPUS="$2"; shift 2 ;;
        --mem)              MEM="$2"; shift 2 ;;
        --time)             TIME="$2"; shift 2 ;;
        --cohort-cpus)      COHORT_CPUS="$2"; shift 2 ;;
        --cohort-mem)       COHORT_MEM="$2"; shift 2 ;;
        --cohort-time)      COHORT_TIME="$2"; shift 2 ;;
        --max-concurrent)   MAX_CONCURRENT="$2"; shift 2 ;;
        --limit)            LIMIT="$2"; shift 2 ;;
        --extra-args)       EXTRA_ARGS="$2"; shift 2 ;;
        --pair-bp)          COHORT_PAIR_BP="$2"; shift 2 ;;
        --cohort-extra)     COHORT_EXTRA_ARGS="$2"; shift 2 ;;
        --env-cmd)          ENV_SETUP_CMD="$2"; shift 2 ;;
        --no-cohort)        NO_COHORT=1; shift ;;
        --write-cohort-gfa) WRITE_COHORT_GFA=1; shift ;;
        --keep-tmp)         KEEP_TMP=1; shift ;;
        --force)            FORCE=1; shift ;;
        --resubmit-missing) RESUBMIT_MISSING=1; shift ;;
        --dry-run)          DRY_RUN=1; shift ;;
        -h|--help)          usage; exit 0 ;;
        *)                  echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$SAMPLES_TSV" || -z "$OUT_DIR" ]]; then
    echo "ERROR: --samples and --out are required." >&2
    usage
    exit 1
fi

SAMPLES_ABS="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$SAMPLES_TSV")"
OUT_ABS="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$OUT_DIR")"
GENES_ABS=""
if [[ -n "$GENES_FA" ]]; then
    GENES_ABS="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$GENES_FA")"
    if [[ ! -s "$GENES_ABS" ]]; then
        echo "ERROR: --genes file not found or empty: $GENES_ABS" >&2
        exit 1
    fi
fi

# Validate sample sheet before touching Slurm
python3 "${SCRIPT_DIR}/make_samplesheet.py" --check "$SAMPLES_ABS"

# Check that a gene source exists
HAS_GCOL="$(python3 -c 'import sys; h=open(sys.argv[1]).readline().rstrip("\r\n").split("\t"); print("1" if "genes_fa" in [c.strip() for c in h] else "0")' "$SAMPLES_ABS")"
if [[ -z "$GENES_ABS" && -z "$KMA_DB" && "$HAS_GCOL" != "1" ]]; then
    echo "ERROR: provide --genes FILE, --kma-db PREFIX, or a 'genes_fa' column in $SAMPLES_ABS" >&2
    exit 1
fi

# Count total valid sample rows
N_TOTAL="$(python3 - "$SAMPLES_ABS" <<'PY'
import sys
with open(sys.argv[1]) as fh:
    fh.readline()
    print(sum(1 for l in fh if l.strip() and not l.startswith("#")))
PY
)"

if [[ "$LIMIT" -gt 0 && "$LIMIT" -lt "$N_TOTAL" ]]; then
    N_RUN="$LIMIT"
else
    N_RUN="$N_TOTAL"
fi

# Determine Slurm --array specification
if [[ "$RESUBMIT_MISSING" == "1" ]]; then
    ARRAY_SPEC="$(python3 - "$SAMPLES_ABS" "$OUT_ABS" "$N_RUN" <<'PY'
import os, sys
tsv, out_dir, limit = sys.argv[1], sys.argv[2], int(sys.argv[3])
miss = []
with open(tsv) as fh:
    hdr = [c.strip() for c in fh.readline().rstrip("\r\n").split("\t")]
    scol = hdr.index("sample")
    idx = 0
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        idx += 1
        if idx > limit:
            break
        sid = line.rstrip("\r\n").split("\t")[scol].strip()
        sdir = os.path.join(out_dir, "samples", sid)
        if not (os.path.isfile(os.path.join(sdir, f"{sid}.done")) and os.path.isfile(os.path.join(sdir, f"{sid}.flanks.json"))):
            miss.append(str(idx))
print(",".join(miss))
PY
)"
    if [[ -z "$ARRAY_SPEC" ]]; then
        echo "All ${N_RUN} samples in ${OUT_ABS} are already complete! Nothing to resubmit."
        exit 0
    fi
    ARRAY_SPEC="${ARRAY_SPEC}%${MAX_CONCURRENT}"
else
    ARRAY_SPEC="1-${N_RUN}%${MAX_CONCURRENT}"
fi

LOG_DIR="${OUT_ABS}/logs"
mkdir -p "$OUT_ABS" "$LOG_DIR"

COMMON_SBATCH=()
[[ -n "$PARTITION" ]] && COMMON_SBATCH+=(--partition="$PARTITION")
[[ -n "$ACCOUNT" ]]   && COMMON_SBATCH+=(--account="$ACCOUNT")
[[ -n "$QOS" ]]       && COMMON_SBATCH+=(--qos="$QOS")

# Build export variables cleanly (use environment inheritance so commas/spaces in EXTRA_ARGS or ENV_SETUP_CMD never break --export)
export REPO_ROOT="$REPO_ROOT"
export HPC_DIR="$SCRIPT_DIR"
export SAMPLES_TSV="$SAMPLES_ABS"
export OUT_DIR="$OUT_ABS"
export GENES_FA="$GENES_ABS"
export KMA_DB="$KMA_DB"
export KMA_REF_FA="$KMA_REF_FA"
export EXTRA_ARGS="$EXTRA_ARGS"
export COHORT_PAIR_BP="$COHORT_PAIR_BP"
export COHORT_EXTRA_ARGS="$COHORT_EXTRA_ARGS"
export WRITE_COHORT_GFA="$WRITE_COHORT_GFA"
export SKIP_COHORT="$NO_COHORT"
export KEEP_TMP="$KEEP_TMP"
export FORCE="$FORCE"
export ENV_SETUP_CMD="$ENV_SETUP_CMD"

ARRAY_CMD=(
    sbatch --parsable
    ${COMMON_SBATCH[@]+"${COMMON_SBATCH[@]}"}
    --job-name="argjacent"
    --array="$ARRAY_SPEC"
    --cpus-per-task="$CPUS"
    --mem="$MEM"
    --time="$TIME"
    --output="${LOG_DIR}/argjacent_%A_%a.out"
    --error="${LOG_DIR}/argjacent_%A_%a.err"
    --export=ALL
    "${SCRIPT_DIR}/run_array.sbatch"
)

COHORT_CMD_BASE=(
    sbatch --parsable
    ${COMMON_SBATCH[@]+"${COMMON_SBATCH[@]}"}
    --job-name="argjacent_cohort"
    --cpus-per-task="$COHORT_CPUS"
    --mem="$COHORT_MEM"
    --time="$COHORT_TIME"
    --output="${LOG_DIR}/argjacent_cohort_%j.out"
    --error="${LOG_DIR}/argjacent_cohort_%j.err"
    --export=ALL
)

echo "=== ARGjacent Slurm Submission Plan ==="
echo "  Samples to run : ${N_RUN} / ${N_TOTAL} (array=${ARRAY_SPEC})"
echo "  Output dir     : ${OUT_ABS}"
echo "  Gene source    : ${GENES_ABS:-${KMA_DB:-(per-sample genes_fa column)}}"
echo "  Array resources: ${CPUS} CPUs, ${MEM} RAM, ${TIME} walltime"
echo "  Cohort job     : $( [[ "$NO_COHORT" == "1" ]] && echo "summary only (--no-cohort)" || echo "argjacent-cohort + summary (${COHORT_CPUS} CPUs, ${COHORT_MEM} RAM, ${COHORT_TIME})" )"

if [[ "$DRY_RUN" == "1" ]]; then
    echo ""
    echo "[DRY RUN] Array command:"
    printf "  %q" "${ARRAY_CMD[@]}"; echo ""
    echo "[DRY RUN] Dependent cohort+summary command:"
    printf "  %q" "${COHORT_CMD_BASE[@]}" "--dependency=afterany:<ARRAY_JOB_ID>" "${SCRIPT_DIR}/run_cohort.sbatch"; echo ""
    exit 0
fi

if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: 'sbatch' not found on PATH. Run this on your Slurm HPC login node (or use --dry-run to preview)." >&2
    exit 1
fi

ARR_JOB_ID="$("${ARRAY_CMD[@]}")"
COHORT_JOB_ID="$("${COHORT_CMD_BASE[@]}" --dependency="afterany:${ARR_JOB_ID}" "${SCRIPT_DIR}/run_cohort.sbatch")"

echo ""
echo "Submitted!"
echo "  Array job ID          : ${ARR_JOB_ID}"
echo "  Cohort + summary job  : ${COHORT_JOB_ID} (runs automatically after array ${ARR_JOB_ID} finishes)"
echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f ${LOG_DIR}/argjacent_${ARR_JOB_ID}_1.out"
echo "  cat ${OUT_ABS}/summary/run_metrics.tsv   # once cohort job finishes"
