#!/usr/bin/env bash
#SBATCH --job-name=argj_global
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/argj_global_%A_%a.out
#SBATCH --error=logs/argj_global_%A_%a.err

# Ephemeral-scratch Slurm array worker for Global Metagenome Benchmarking:
#   1. Fetches paired-end FASTQs (from ENA/SRA URLs in global_manifest.tsv) into node-local $SLURM_TMPDIR
#   2. Optionally trims reads with fastp and deletes raw FASTQs immediately
#   3. Runs per-sample KMA pre-screen (or uses shared GENES_FA) -> produces <sample>.genes.fa
#   4. Runs ARGjacent (±10 kb certified flanks + GFA + pairs + depth_bins) and deletes FASTQs immediately
#   5. Fetches the matched pre-made ERZ (or Logan) assembly into $SLURM_TMPDIR, extracts ±10 kb flanks
#      around the exact same <sample>.genes.fa, runs per-sample comparison, and deletes the assembly
#   6. Atomically stages only compact outputs (~1-15 MB/sample) back to $OUT_DIR/samples/<sample>/

set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

if [[ -n "${ENV_SETUP_CMD:-}" ]]; then
    eval "$ENV_SETUP_CMD"
fi

MANIFEST_TSV="${MANIFEST_TSV:?Set MANIFEST_TSV to manifest produced by hpc/fetch_global_manifest.py}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to output directory}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?Set SLURM_ARRAY_TASK_ID (1-based row index)}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
GENES_FA="${GENES_FA:-}"
KMA_DB="${KMA_DB:-}"
KMA_REF_FA="${KMA_REF_FA:-}"
KMA_MIN_ID="${KMA_MIN_ID:-90}"
KMA_MIN_COV="${KMA_MIN_COV:-60}"
FLANK_BP="${FLANK_BP:-10000}"
SKIP_FASTP="${SKIP_FASTP:-0}"
SKIP_ASSEMBLY="${SKIP_ASSEMBLY:-0}"
USE_LOGAN="${USE_LOGAN:-0}"
ASM_ALIGNER="${ASM_ALIGNER:-auto}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
FORCE="${FORCE:-0}"

SCRIPT_DIR="${HPC_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY_BIN="${PYTHON_BIN:-python3}"

# Locate argjacent executable
if [[ -n "${ARGJACENT_BIN:-}" ]]; then
    ARGJACENT="$ARGJACENT_BIN"
elif command -v argjacent >/dev/null 2>&1; then
    ARGJACENT="$(command -v argjacent)"
elif [[ -n "${REPO_ROOT:-}" && -x "${REPO_ROOT}/bin/argjacent" ]]; then
    ARGJACENT="${REPO_ROOT}/bin/argjacent"
elif [[ -x "${SCRIPT_DIR}/../bin/argjacent" ]]; then
    ARGJACENT="${SCRIPT_DIR}/../bin/argjacent"
else
    echo "ERROR: argjacent executable not found." >&2
    exit 1
fi

# Read manifest row TASK_ID (supports both remote URLs and local file paths for offline testing)
ROW_DATA="$("$PY_BIN" - "$MANIFEST_TSV" "$TASK_ID" <<'PY'
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
            def get(k, alt=""):
                i = col.get(k)
                if i is not None and i < len(parts) and parts[i].strip():
                    return parts[i].strip()
                return alt
            sid = get("sample")
            r1 = get("r1_url") or get("r1")
            r2 = get("r2_url") or get("r2")
            erz = get("erz_accession")
            asm = get("assembly_url") or get("assembly")
            logan = get("logan_url")
            gf = get("genes_fa")
            print("\t".join([sid, r1, r2, erz, asm, logan, gf]))
            sys.exit(0)
sys.exit(2)
PY
)" || { echo "ERROR: task index ${TASK_ID} not found in ${MANIFEST_TSV}" >&2; exit 1; }

IFS=$'\t' read -r SAMPLE R1_SRC R2_SRC ERZ_ACC ASM_SRC LOGAN_SRC ROW_GENES <<< "$ROW_DATA"
SAMPLE_OUT_DIR="${OUT_DIR}/samples/${SAMPLE}"
DONE_MARK="${SAMPLE_OUT_DIR}/${SAMPLE}.done"

echo "=== ARGjacent Global Benchmark Task ${TASK_ID}: sample=${SAMPLE} erz=${ERZ_ACC:-none} host=$(hostname) threads=${THREADS} ==="

if [[ "$FORCE" != "1" && -f "$DONE_MARK" && -s "${SAMPLE_OUT_DIR}/${SAMPLE}.flanks.json" ]]; then
    echo "  Already completed (${DONE_MARK} exists). Skipping."
    exit 0
fi

# Derive --max-mem in GB (80% of Slurm allocation)
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

SCRATCH_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
WORK_DIR="$(mktemp -d "${SCRATCH_BASE}/argj_glob_${SAMPLE}_XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT
mkdir -p "$SAMPLE_OUT_DIR" "${WORK_DIR}/raw" "${WORK_DIR}/trim" "${WORK_DIR}/asm"

fetch_file() {
    local src="$1" dst="$2"
    if [[ -f "$src" ]]; then
        ln -sf "$(cd "$(dirname "$src")" && pwd)/$(basename "$src")" "$dst"
        return 0
    fi
    if command -v curl >/dev/null 2>&1; then
        if ! curl -fL -C - --retry 2 --retry-delay 5 --connect-timeout 15 --max-time 3600 -sS -o "$dst" "$src"; then
            if [[ "$src" == https://ftp.sra.ebi.ac.uk/* ]]; then
                local alt="http://${src#https://}"
                curl -fL -C - --retry 5 --retry-delay 8 --retry-all-errors --connect-timeout 30 --max-time 3600 -sS -o "$dst" "$alt"
            else
                return 1
            fi
        fi
    elif command -v wget >/dev/null 2>&1; then
        wget -q -c --tries=5 --waitretry=5 --timeout=30 -O "$dst" "$src"
    else
        echo "ERROR: neither curl nor wget found to fetch ${src}" >&2
        return 1
    fi
}

# 1. Fetch paired-end FASTQs into ephemeral node-local scratch
T0=$(date +%s)
RAW_R1="${WORK_DIR}/raw/${SAMPLE}_1.fastq.gz"
RAW_R2="${WORK_DIR}/raw/${SAMPLE}_2.fastq.gz"
echo "  [1/5] Fetching paired FASTQs..."
fetch_file "$R1_SRC" "$RAW_R1"
fetch_file "$R2_SRC" "$RAW_R2"
T_DL=$(( $(date +%s) - T0 ))
echo "        Downloaded R1+R2 in ${T_DL} s."

# 2. Optional fastp trimming (and immediate deletion of raw FASTQs to halve scratch usage)
if [[ "$SKIP_FASTP" != "1" ]] && command -v fastp >/dev/null 2>&1; then
    echo "  [2/5] Trimming with fastp..."
    T1=$(date +%s)
    USE_R1="${WORK_DIR}/trim/${SAMPLE}_1.fq.gz"
    USE_R2="${WORK_DIR}/trim/${SAMPLE}_2.fq.gz"
    fastp -i "$RAW_R1" -I "$RAW_R2" -o "$USE_R1" -O "$USE_R2" \
        --detect_adapter_for_pe -q 20 -l 80 -w "$THREADS" \
        -j "${WORK_DIR}/${SAMPLE}.fastp.json" -h /dev/null >/dev/null 2>&1
    rm -f "$RAW_R1" "$RAW_R2"
    echo "        fastp done in $(( $(date +%s) - T1 )) s."
else
    echo "  [2/5] Skipping fastp (SKIP_FASTP=${SKIP_FASTP} or fastp not installed)."
    USE_R1="$RAW_R1"
    USE_R2="$RAW_R2"
fi

# 3. Determine target gene FASTA (per-row > KMA pre-screen > shared GENES_FA)
ACTIVE_GENES=""
N_GENES=-1
if [[ -n "$ROW_GENES" ]]; then
    ACTIVE_GENES="$ROW_GENES"
elif [[ -n "$KMA_DB" ]]; then
    echo "  [3/5] Running KMA pre-screen against ${KMA_DB}..."
    T2=$(date +%s)
    KMA_PRE="${WORK_DIR}/kma_out"
    kma -ipe "$USE_R1" "$USE_R2" -o "$KMA_PRE" -t_db "$KMA_DB" -t "$THREADS" -1t1 -cge -ef -nf >/dev/null 2>&1 || \
        kma -ipe "$USE_R1" "$USE_R2" -o "$KMA_PRE" -t_db "$KMA_DB" -t "$THREADS" -1t1 -cge -ef -nf
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
            name, tid, tcov = parts[0].split()[0], float(parts[4]), float(parts[5])
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
                cur, buf = line[1:].strip().split()[0], []
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
    echo "        KMA detected ${N_GENES} genes in $(( $(date +%s) - T2 )) s."
    if [[ "$N_GENES" -eq 0 ]]; then
        printf "gene\tside\tstatus\tcontexts\tcertain_bp\tlongest_bp\tcontexts_certain_1kb\tread_depth\tgene_end\tcertain_end_bp\n" > "${SAMPLE_OUT_DIR}/${SAMPLE}.summary.tsv"
        printf "gene\tside\tcontext\tlength_bp\tcertain_bp\tread_depth\tstop\tcertain_stop\tother_genes\tgene_end\tmain\tcertain_end_bp\tshared_end_bp\n" > "${SAMPLE_OUT_DIR}/${SAMPLE}.contexts.tsv"
        : > "${SAMPLE_OUT_DIR}/${SAMPLE}.flanks.fa"
        printf "gene\tpair_id\tleft_context\tright_context\tleft_certain_bp\tright_certain_bp\tleft_depth\tright_depth\ttier\tscore\tnote\n" > "${SAMPLE_OUT_DIR}/${SAMPLE}.pairs.tsv"
        printf '{"version": "0.5.1", "flank_bp": %s, "aw_stats": {}, "parameters": {}, "genes": {}}\n' "$FLANK_BP" > "${SAMPLE_OUT_DIR}/${SAMPLE}.flanks.json"
        printf '{"version": "asm_flanks_v1", "source": "%s", "flank_bp": %s, "genes": {}}\n' "${ERZ_ACC:-none}" "$FLANK_BP" > "${SAMPLE_OUT_DIR}/${SAMPLE}.asm_flanks.json"
        printf "argjacent global hpc: 0 genes detected by KMA pre-screen\n" > "${SAMPLE_OUT_DIR}/${SAMPLE}.log"
        date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE_MARK"
        exit 0
    fi
    cp "$DETECTED_FA" "${SAMPLE_OUT_DIR}/${SAMPLE}.genes.fa"
    ACTIVE_GENES="${SAMPLE_OUT_DIR}/${SAMPLE}.genes.fa"
elif [[ -n "$GENES_FA" ]]; then
    ACTIVE_GENES="$GENES_FA"
else
    echo "ERROR: No gene source provided (set GENES_FA, KMA_DB, or genes_fa column)." >&2
    exit 1
fi

# 4. Run ARGjacent and immediately purge FASTQs from scratch
LOCAL_PRE="${WORK_DIR}/${SAMPLE}"
CMD=("$ARGJACENT" -g "$ACTIVE_GENES" -1 "$USE_R1" -2 "$USE_R2" -o "$LOCAL_PRE" -D "$FLANK_BP" -t "$THREADS")
[[ -n "${MAX_MEM_GB:-}" ]] && CMD+=(--max-mem "$MAX_MEM_GB")
if [[ -n "$EXTRA_ARGS" ]]; then
    # shellcheck disable=SC2206
    EXTRA_ARR=($EXTRA_ARGS)
    CMD+=("${EXTRA_ARR[@]}")
fi

echo "  [4/5] Running ARGjacent: ${CMD[*]}"
T3=$(date +%s)
"${CMD[@]}"
T_ARGJ=$(( $(date +%s) - T3 ))
echo "        ARGjacent finished in ${T_ARGJ} s. Deleting local FASTQs."
rm -rf "${WORK_DIR}/raw" "${WORK_DIR}/trim"

# 5. Fetch pre-made assembly (ERZ or Logan), extract +/-D bp flanks, compare, and delete assembly
ASM_TARGET="$ASM_SRC"
ASM_LABEL="${ERZ_ACC:-assembly}"
if [[ "$USE_LOGAN" == "1" && -n "$LOGAN_SRC" ]]; then
    ASM_TARGET="$LOGAN_SRC"
    ASM_LABEL="logan:${SAMPLE}"
fi

if [[ "$SKIP_ASSEMBLY" != "1" && -n "$ASM_TARGET" ]]; then
    echo "  [5/5] Fetching pre-made assembly (${ASM_LABEL}) and extracting ±${FLANK_BP} bp flanks..."
    T4=$(date +%s)
    ASM_EXT=".fasta.gz"
    case "$ASM_TARGET" in
        *.zst) ASM_EXT=".fa.zst" ;;
        *.fa|*.fna|*.fasta) ASM_EXT=".fa" ;;
    esac
    LOCAL_ASM="${WORK_DIR}/asm/${SAMPLE}${ASM_EXT}"
    fetch_file "$ASM_TARGET" "$LOCAL_ASM"
    "$PY_BIN" "${SCRIPT_DIR}/extract_assembly_flanks.py" \
        -g "$ACTIVE_GENES" -c "$LOCAL_ASM" -o "$LOCAL_PRE" \
        -D "$FLANK_BP" --aligner "$ASM_ALIGNER" -t "$THREADS" --source-label "$ASM_LABEL"
    rm -rf "${WORK_DIR}/asm"
    "$PY_BIN" "${SCRIPT_DIR}/compare_argjacent_vs_assembly.py" \
        --argj-json "${LOCAL_PRE}.flanks.json" --asm-json "${LOCAL_PRE}.asm_flanks.json" \
        --sample "$SAMPLE" -o "${LOCAL_PRE}.cmp" >/dev/null
    echo "        Assembly extract + compare done in $(( $(date +%s) - T4 )) s."
else
    echo "  [5/5] Skipping pre-made assembly comparison (SKIP_ASSEMBLY=${SKIP_ASSEMBLY}, url='${ASM_TARGET}')."
fi

# Atomically stage compact outputs back to shared filesystem
for ext in .summary.tsv .contexts.tsv .flanks.fa .flanks.json .pairs.tsv .gfa .log .fastp.json \
           .asm_flanks.json .asm_contexts.tsv .asm_flanks.fa .cmp_per_side.tsv .cmp_summary.tsv; do
    if [[ -f "${LOCAL_PRE}${ext}" ]]; then
        cp "${LOCAL_PRE}${ext}" "${SAMPLE_OUT_DIR}/.${SAMPLE}${ext}.part"
        mv -f "${SAMPLE_OUT_DIR}/.${SAMPLE}${ext}.part" "${SAMPLE_OUT_DIR}/${SAMPLE}${ext}"
    fi
done

printf "%s\n" "$ACTIVE_GENES" > "${SAMPLE_OUT_DIR}/${SAMPLE}.genes_path.txt"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$DONE_MARK"
echo "=== Completed ${SAMPLE} (dl=${T_DL}s argj=${T_ARGJ}s) ==="
