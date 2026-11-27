#!/usr/bin/env bash
#SBATCH --job-name=argj_global_post
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/argj_global_post_%j.out
#SBATCH --error=logs/argj_global_post_%j.err

# Post-array aggregation job for Global Metagenome Benchmarking:
#   1. Runs `hpc/summarize_run.py` to merge per-sample tables and write summary/run_metrics.tsv
#   2. Runs per-BioProject `argjacent-cohort` across every study with >= COHORT_MIN_SAMPLES completed runs
#   3. Runs global `hpc/compare_argjacent_vs_assembly.py --export-fasta` to build biome-stratified
#      head-to-head comparison tables + matched certified/assembly FASTAs for unified MGE/plasmid annotation

set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ -n "${ENV_SETUP_CMD:-}" ]]; then
    eval "$ENV_SETUP_CMD"
fi

MANIFEST_TSV="${MANIFEST_TSV:?Set MANIFEST_TSV to global manifest TSV}"
OUT_DIR="${OUT_DIR:?Set OUT_DIR to output directory}"
GENES_FA="${GENES_FA:-}"
COHORT_PAIR_BP="${COHORT_PAIR_BP:-500}"
COHORT_MIN_SAMPLES="${COHORT_MIN_SAMPLES:-4}"
SKIP_COHORT="${SKIP_COHORT:-0}"
SKIP_ASSEMBLY="${SKIP_ASSEMBLY:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

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
    COHORT_BIN=""
fi

mkdir -p "${OUT_DIR}/summary" "${OUT_DIR}/cohorts" "${OUT_DIR}/comparison"

echo "=== ARGjacent Global Benchmark Post-Processing (host=$(hostname)) ==="

# 1. Merge standard per-sample tables + runtime/RAM metrics
"$PY_BIN" "${SCRIPT_DIR}/summarize_run.py" \
    --samples "$MANIFEST_TSV" \
    --out-dir "$OUT_DIR" \
    --summary-dir "${OUT_DIR}/summary"

# 2. Per-BioProject cohort deconvolution (only group samples from the same BioProject!)
if [[ "$SKIP_COHORT" != "1" && -n "$COHORT_BIN" ]]; then
    echo ""
    echo "-- Building per-BioProject cohort sheets (min_samples=${COHORT_MIN_SAMPLES}) --"
    STUDY_LIST="$("$PY_BIN" - "$MANIFEST_TSV" "$OUT_DIR" "$GENES_FA" "$COHORT_MIN_SAMPLES" <<'PY'
import collections, json, os, sys
man, out_dir, shared_g, min_n = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
by_study = collections.defaultdict(list)
with open(man) as fh:
    hdr = [c.strip() for c in fh.readline().rstrip("\r\n").split("\t")]
    col = {c: i for i, c in enumerate(hdr)}
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        p = line.rstrip("\r\n").split("\t")
        sid = p[col["sample"]].strip()
        st = p[col["study_accession"]].strip() if "study_accession" in col and col["study_accession"] < len(p) else "ALL"
        st = st or "ALL"
        sdir = os.path.join(out_dir, "samples", sid)
        jp = os.path.join(sdir, f"{sid}.flanks.json")
        if not (os.path.isfile(jp) and os.path.getsize(jp) > 0):
            continue
        try:
            doc = json.load(open(jp))
        except Exception:
            continue
        if not doc.get("genes"):
            continue
        gp = ""
        gfile = os.path.join(sdir, f"{sid}.genes_path.txt")
        if os.path.isfile(gfile):
            cand = open(gfile).read().strip()
            if cand and os.path.isfile(cand):
                gp = cand
        if not gp:
            gp = shared_g
        by_study[st].append((sid, os.path.abspath(jp), os.path.abspath(gp) if gp else "", ""))

for st, rows in sorted(by_study.items()):
    if len(rows) < min_n:
        continue
    cdir = os.path.join(out_dir, "cohorts", st)
    os.makedirs(cdir, exist_ok=True)
    sheet = os.path.join(cdir, "cohort_samplesheet.tsv")
    with open(sheet, "w") as out:
        out.write("sample\tflanks_json\tgenes_fa\treads_tsv\n")
        for r in rows:
            out.write("\t".join(r) + "\n")
    print(f"{st}\t{len(rows)}\t{sheet}")
PY
)"
    if [[ -n "$STUDY_LIST" ]]; then
        while IFS=$'\t' read -r STUDY N_S SHEET; do
            echo "  Running argjacent-cohort on BioProject ${STUDY} (${N_S} samples)..."
            CMD=("$COHORT_BIN" --sample-sheet "$SHEET" -o "${OUT_DIR}/cohorts/${STUDY}/${STUDY}" --pair-bp "$COHORT_PAIR_BP" -t "$THREADS")
            [[ -n "$GENES_FA" ]] && CMD+=(-g "$GENES_FA")
            "${CMD[@]}"
        done <<< "$STUDY_LIST"
    else
        echo "  No BioProjects reached >= ${COHORT_MIN_SAMPLES} completed samples with detected ARGs."
    fi
fi

# 3. Global head-to-head comparison against pre-made assemblies + unified FASTA exports
if [[ "$SKIP_ASSEMBLY" != "1" ]]; then
    echo ""
    "$PY_BIN" "${SCRIPT_DIR}/compare_argjacent_vs_assembly.py" \
        --manifest "$MANIFEST_TSV" \
        --out-dir "$OUT_DIR" \
        -o "${OUT_DIR}/comparison/global_cmp" \
        --export-fasta
fi

echo "=== Global Benchmark Post-Processing Completed ==="
