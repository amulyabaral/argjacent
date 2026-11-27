#!/usr/bin/env bash
#SBATCH --job-name=argj_annotate
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/argj_annotate_%j.out
#SBATCH --error=logs/argj_annotate_%j.err

# Slurm job for unified 5-layer annotation & comparison of ARGjacent certified flanks vs.
# pre-made metagenome assembly flanks produced by `hpc/run_global_post.sbatch`.
#
# Usage (standalone or chained via --dependency=afterok:<POST_JOB_ID>):
#   sbatch --export=ALL,OUT_DIR=/scratch/$USER/argj_global,MANIFEST_TSV=global_2023.tsv,DB_ENV=/scratch/$USER/argj_dbs/db_paths.env \
#       hpc/run_annotate.sbatch

set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

if [[ -n "${ENV_SETUP_CMD:-}" ]]; then
    eval "$ENV_SETUP_CMD"
fi

if [[ -n "${DB_ENV:-}" && -f "${DB_ENV}" ]]; then
    # shellcheck disable=SC1090
    source "${DB_ENV}"
fi

OUT_DIR="${OUT_DIR:?Set OUT_DIR to your global benchmark output directory}"
MANIFEST_TSV="${MANIFEST_TSV:-}"
ARGJ_FA="${ARGJ_FA:-${OUT_DIR}/comparison/global_cmp_argjacent_certified.fa}"
ASM_FA="${ASM_FA:-${OUT_DIR}/comparison/global_cmp_assembly_flanks.fa}"
ANN_OUT="${ANN_OUT:-${OUT_DIR}/annotation}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
MIN_LEN="${MIN_LEN:-1000}"
CLUSTER_HGT="${CLUSTER_HGT:-1}"

SCRIPT_DIR="${HPC_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY_BIN="${PYTHON_BIN:-python3}"

CMD=("$PY_BIN" "${SCRIPT_DIR}/annotate_flanks.py" -o "$ANN_OUT" -t "$THREADS" --min-len "$MIN_LEN")

[[ -s "$ARGJ_FA" ]] && CMD+=(--argj-fa "$ARGJ_FA")
[[ -s "$ASM_FA" ]]  && CMD+=(--asm-fa "$ASM_FA")
[[ -n "$MANIFEST_TSV" && -s "$MANIFEST_TSV" ]] && CMD+=(--manifest "$MANIFEST_TSV")

[[ -n "${UNIVEC_DB:-}" ]]        && CMD+=(--univec-db "$UNIVEC_DB")
[[ -n "${ARG_DB:-}" ]]           && CMD+=(--arg-db "$ARG_DB")
[[ -n "${ISFINDER_DB:-}" ]]      && CMD+=(--isfinder-db "$ISFINDER_DB")
[[ -n "${MOBILEOG_DB:-}" ]]      && CMD+=(--mobileog-db "$MOBILEOG_DB")
[[ -n "${PLASMIDFINDER_DB:-}" ]] && CMD+=(--plasmidfinder-db "$PLASMIDFINDER_DB")
[[ -n "${GENOMAD_DB:-}" ]]       && CMD+=(--genomad-db "$GENOMAD_DB")
[[ -n "${KRAKEN2_DB:-}" ]]       && CMD+=(--kraken2-db "$KRAKEN2_DB")
[[ -n "${MMSEQS_TAX_DB:-}" ]]    && CMD+=(--mmseqs-tax-db "$MMSEQS_TAX_DB")
[[ "$CLUSTER_HGT" == "1" ]]      && CMD+=(--cluster-hgt)

if [[ -n "${ANNOTATE_EXTRA_ARGS:-}" ]]; then
    # shellcheck disable=SC2206
    EXTRA_ARR=($ANNOTATE_EXTRA_ARGS)
    CMD+=("${EXTRA_ARR[@]}")
fi

echo "=== ARGjacent Unified Annotation Job (host=$(hostname) threads=${THREADS}) ==="
echo "  Running: ${CMD[*]}"
START_TS=$(date +%s)
"${CMD[@]}"
echo "=== Annotation job completed in $(( $(date +%s) - START_TS )) s ==="
