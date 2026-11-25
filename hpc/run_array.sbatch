#!/usr/bin/env bash
#SBATCH --job-name=argjacent
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/argjacent_%A_%a.out
#SBATCH --error=logs/argjacent_%A_%a.err

# Slurm array worker for single-sample ARGjacent runs.
# Typically launched via `bash hpc/submit_pipeline.sh ...`, or directly with:
#   sbatch --array=1-N%20 \
#     --export=ALL,SAMPLES_TSV=samples.tsv,GENES_FA=resfinder.fa,OUT_DIR=results \
#     hpc/run_array.sbatch
#
# Key HPC safeguards built in:
#   1. Pins OpenBLAS/OMP/MKL threads to 1 so Python worker processes never oversubscribe cores.
#   2. Stages all temporary I/O (*.tmp/) onto local node scratch ($SLURM_TMPDIR / $TMPDIR),
#      avoiding Lustre/GPFS/NFS metadata storms, and cleans scratch on exit.
#   3. Atomically copies final outputs back and writes a .done marker for safe resume.
#   4. Derives --max-mem from Slurm's allocation automatically (80% headroom).
#   5. Supports shared GENES_FA, per-sample genes_fa column, or on-the-fly KMA pre-screening (KMA_DB).

set -euo pipefail

# Prevent BLAS/OpenMP oversubscription inside argjacent's ProcessPoolExecutor workers
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# Optional user environment hook (e.g. `module load gcc python` or `conda activate argjacent`)
if [[ -n "${ENV_SETUP_CMD:-}" ]]; then
    eval "$ENV_SETUP_CMD"
fi

SAMPLES_TSV="${SAMPLES_TSV:?Set SAMPLES_TSV to your tab-separated sample sheet}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to your output directory}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?Set SLURM_ARRAY_TASK_ID (1-based sample index)}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
GENES_FA="${GENES_FA:-}"
KMA_DB="${KMA_DB:-}"
KMA_REF_FA="${KMA_REF_FA:-}"
KMA_MIN_ID="${KMA_MIN_ID:-90}"
KMA_MIN_COV="${KMA_MIN_COV:-60}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
FORCE="${FORCE:-0}"
KEEP_TMP="${KEEP_TMP:-0}"

# Locate argjacent executable
if [[ -n "${ARGJACENT_BIN:-}" ]]; then
    ARGJACENT="$ARGJACENT_BIN"
elif command -v argjacent >/dev/null 2>&1; then
    ARGJACENT="$(command -v argjacent)"
elif [[ -n "${REPO_ROOT:-}" && -x "${REPO_ROOT}/bin/argjacent" ]]; then
    ARGJACENT="${REPO_ROOT}/bin/argjacent"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    if [[ -x "${SCRIPT_DIR}/../bin/argjacent" ]]; then
        ARGJACENT="${SCRIPT_DIR}/../bin/argjacent"
    else
        echo "ERROR: argjacent executable not found. Install via pip or 'sh install.sh', or set ARGJACENT_BIN." >&2
        exit 1
    fi
fi

# Locate Python 3 (for parsing TSV / KMA hits)
PY_BIN="${PYTHON_BIN:-python3}"
if ! command -v "$PY_BIN" >/dev/null 2>&1; then
    echo "ERROR: python3 not found on PATH." >&2
    exit 1
fi

# Extract row TASK_ID (1-based, skipping header and blank/comment lines)
ROW_DATA="$("$PY_BIN" - "$SAMPLES_TSV" "$TASK_ID" <<'PY'
import sys
tsv, tid = sys.argv[1], int(sys.argv[2])
with open(tsv) as fh:
    hdr = [c.strip() for c in fh.readline().rstrip("\r\n").split("\t")]
    col = {c: i for i, c in enumerate(hdr)}
    idx = 0
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        idx += 1
        if idx == tid:
            parts = line.rstrip("\r\n").split("\t")
            s = parts[col["sample"]].strip()
            r1 = parts[col["r1"]].strip()
            r2 = parts[col["r2"]].strip()
            gf = parts[col["genes_fa"]].strip() if "genes_fa" in col and len(parts) > col["genes_fa"] else ""
            print(f"{s}\t{r1}\t{r2}\t{gf}")
            sys.exit(0)
sys.exit(2)
PY
)" || { echo "ERROR: sample index ${TASK_ID} not found in ${SAMPLES_TSV}" >&2; exit 1; }

IFS=$'\t' read -r SAMPLE R1 R2 ROW_GENES <<< "$ROW_DATA"
SAMPLE_OUT_DIR="${OUT_DIR}/samples/${SAMPLE}"
DONE_MARK="${SAMPLE_OUT_DIR}/${SAMPLE}.done"

echo "=== ARGjacent HPC task ${TASK_ID}: sample=${SAMPLE} host=$(hostname) threads=${THREADS} ==="
echo "  R1: ${R1}"
echo "  R2: ${R2}"
echo "  Out: ${SAMPLE_OUT_DIR}"

if [[ "$FORCE" != "1" && -f "$DONE_MARK" && -s "${SAMPLE_OUT_DIR}/${SAMPLE}.flanks.json" && -s "${SAMPLE_OUT_DIR}/${SAMPLE}.summary.tsv" ]]; then
    echo "  Already completed (${DONE_MARK} exists). Skipping (set FORCE=1 to rerun)."
    exit 0
fi

# Validate input FASTQs
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing FASTQ input for ${SAMPLE}: R1=${R1} R2=${R2}" >&2
    exit 1
fi

# Derive --max-mem in GB (80% of Slurm allocation, leaving 20% headroom for OS + Python workers)
if [[ -z "${MAX_MEM_GB:-}" ]]; then
    MEM_MB=""
    if [[ -n "${SLURM_MEM_PER_NODE:-}" && "${SLURM_MEM_PER_NODE}" =~ ^[0-9]+$ && "${SLURM_MEM_PER_NODE}" -gt 0 ]]; then
        MEM_MB="$SLURM_MEM_PER_NODE"
    elif [[ -n "${SLURM_MEM_PER_CPU:-}" && "${SLURM_MEM_PER_CPU}" =~ ^[0-9]+$ && "${SLURM_MEM_PER_CPU}" -gt 0 ]]; then
        MEM_MB=$((SLURM_MEM_PER_CPU * THREADS))
    fi
    if [[ -n "$MEM_MB" ]]; then
        MAX_MEM_GB="$("$PY_BIN" -c "print(f'{max(1.0, ($MEM_MB / 1024.0) * 0.80):.1f}')")"
    fi
fi

# Create node-local scratch workspace (cleaned automatically on exit/SIGTERM)
SCRATCH_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
WORK_DIR="$(mktemp -d "${SCRATCH_BASE}/argj_${SAMPLE}_XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

mkdir -p "$SAMPLE_OUT_DIR"

# Determine gene FASTA source: per-sample column > KMA_DB pre-screen > shared GENES_FA
ACTIVE_GENES=""
if [[ -n "$ROW_GENES" ]]; then
    ACTIVE_GENES="$ROW_GENES"
elif [[ -n "$KMA_DB" ]]; then
    if ! command -v kma >/dev/null 2>&1; then
        echo "ERROR: KMA_DB was set (${KMA_DB}) but 'kma' is not on PATH." >&2
        exit 1
    fi
    echo "  [Pre-screen] Running KMA against ${KMA_DB} (id>=${KMA_MIN_ID}%, cov>=${KMA_MIN_COV}%) ..."
    KMA_PRE="${WORK_DIR}/kma_out"
    kma -ipe "$R1" "$R2" -o "$KMA_PRE" -t_db "$KMA_DB" -t "$THREADS" -1t1 -cge -ef -nf >/dev/null 2>&1 || \
        kma -ipe "$R1" "$R2" -o "$KMA_PRE" -t_db "$KMA_DB" -t "$THREADS" -1t1 -cge -ef -nf
    DETECTED_FA="${WORK_DIR}/${SAMPLE}.kma_genes.fa"
    N_GENES="$("$PY_BIN" - "$KMA_PRE.res" "$KMA_PRE.fsa" "$DETECTED_FA" "$KMA_MIN_ID" "$KMA_MIN_COV" "$KMA_REF_FA" <<'PY'
import sys, os
res_p, fsa_p, out_p, min_id, min_cov, ref_p = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
kept = []
if os.path.isfile(res_p):
    with open(res_p) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 6:
                continue
            name = parts[0].split()[0]
            tid, tcov = float(parts[4]), float(parts[5])
            if tid >= min_id and tcov >= min_cov:
                kept.append(name)
kept_set = set(kept)
seqs = {}
src = ref_p if (ref_p and os.path.isfile(ref_p)) else fsa_p
if kept_set and os.path.isfile(src):
    cur, buf = None, []
    with open(src) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur in kept_set:
                    seqs[cur] = "".join(buf).upper().replace("-", "")
                cur = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur in kept_set:
            seqs[cur] = "".join(buf).upper().replace("-", "")
n = 0
with open(out_p, "w") as out:
    for name in kept:
        s = seqs.get(name, "")
        if len(s) >= 31:
            out.write(f">{name}\n{s}\n")
            n += 1
print(n)
PY
)"
    echo "  [Pre-screen] KMA detected ${N_GENES} resistance gene(s)."
    if [[ "$N_GENES" -eq 0 ]]; then
        # Clean zero-gene completion so the array job succeeds and summary marks 0 genes cleanly
        printf "gene\tside\tstatus\tcontexts\tcertain_bp\tlongest_bp\tcontexts_certain_1kb\tread_depth\tgene_end\tcertain_end_bp\n" > "${SAMPLE_OUT_DIR}/${SAMPLE}.summary.tsv"
        printf "gene\tside\tcontext\tlength_bp\tcertain_bp\tread_depth\tstop\tcertain_stop\tother_genes\tgene_end\tmain\tcertain_end_bp\tshared_end_bp\n" > "${SAMPLE_OUT_DIR}/${SAMPLE}.contexts.tsv"
        : > "${SAMPLE_OUT_DIR}/${SAMPLE}.flanks.fa"
        printf "gene\tpair_id\tleft_context\tright_context\tleft_certain_bp\tright_certain_bp\tleft_depth\tright_depth\ttier\tscore\tnote\n" > "${SAMPLE_OUT_DIR}/${SAMPLE}.pairs.tsv"
        ARGJ_VER="$("$ARGJACENT" --version 2>/dev/null | awk '{print $2}' || echo "0.5.1")"
        printf '{"version": "%s", "flank_bp": 10000, "aw_stats": {}, "parameters": {}, "genes": {}}\n' "${ARGJ_VER:-0.5.1}" > "${SAMPLE_OUT_DIR}/${SAMPLE}.flanks.json"
        printf "argjacent hpc: 0 genes detected by KMA pre-screen\n" > "${SAMPLE_OUT_DIR}/${SAMPLE}.log"
        date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE_MARK"
        exit 0
    fi
    cp "$DETECTED_FA" "${SAMPLE_OUT_DIR}/${SAMPLE}.genes.fa"
    ACTIVE_GENES="${SAMPLE_OUT_DIR}/${SAMPLE}.genes.fa"
elif [[ -n "$GENES_FA" ]]; then
    ACTIVE_GENES="$GENES_FA"
else
    echo "ERROR: No genes FASTA provided (set GENES_FA, KMA_DB, or genes_fa column in samples.tsv)." >&2
    exit 1
fi

if [[ ! -s "$ACTIVE_GENES" ]]; then
    echo "ERROR: Gene FASTA missing or empty: ${ACTIVE_GENES}" >&2
    exit 1
fi

# Assemble command
LOCAL_PRE="${WORK_DIR}/${SAMPLE}"
CMD=("$ARGJACENT" -g "$ACTIVE_GENES" -1 "$R1" -2 "$R2" -o "$LOCAL_PRE" -t "$THREADS")
if [[ -n "${MAX_MEM_GB:-}" ]]; then
    CMD+=(--max-mem "$MAX_MEM_GB")
fi
if [[ "$KEEP_TMP" == "1" ]]; then
    CMD+=(--keep)
fi
if [[ -n "$EXTRA_ARGS" ]]; then
    # shellcheck disable=SC2206
    EXTRA_ARR=($EXTRA_ARGS)
    CMD+=("${EXTRA_ARR[@]}")
fi

echo "  Running: ${CMD[*]}"
START_TS=$(date +%s)
"${CMD[@]}"
ELAPSED=$(( $(date +%s) - START_TS ))
echo "  Finished ARGjacent in ${ELAPSED} s. Staging final outputs to ${SAMPLE_OUT_DIR} ..."

# Copy compact final results atomically from local scratch back to shared storage
for ext in .summary.tsv .contexts.tsv .flanks.fa .flanks.json .pairs.tsv .gfa .log .near_1.fq .near_2.fq; do
    if [[ -f "${LOCAL_PRE}${ext}" ]]; then
        cp "${LOCAL_PRE}${ext}" "${SAMPLE_OUT_DIR}/.${SAMPLE}${ext}.part"
        mv -f "${SAMPLE_OUT_DIR}/.${SAMPLE}${ext}.part" "${SAMPLE_OUT_DIR}/${SAMPLE}${ext}"
    fi
done

if [[ "$KEEP_TMP" == "1" && -f "${LOCAL_PRE}.tmp/aw.reads.tsv" ]]; then
    mkdir -p "${SAMPLE_OUT_DIR}/${SAMPLE}.tmp"
    gzip -c "${LOCAL_PRE}.tmp/aw.reads.tsv" > "${SAMPLE_OUT_DIR}/${SAMPLE}.tmp/aw.reads.tsv.gz"
fi

# Record exact genes FASTA used (helpful for argjacent-cohort)
printf "%s\n" "$ACTIVE_GENES" > "${SAMPLE_OUT_DIR}/${SAMPLE}.genes_path.txt"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE_MARK"
echo "=== Completed ${SAMPLE} (${ELAPSED} s) ==="
