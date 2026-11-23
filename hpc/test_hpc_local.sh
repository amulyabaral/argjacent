#!/usr/bin/env bash
# Local end-to-end smoke test for the Slurm HPC scripts (runs without needing Slurm installed).
# Simulates a 3-sample paired-read cohort, tests:
#   1. hpc/make_samplesheet.py discovery & validation
#   2. hpc/submit_pipeline.sh --dry-run
#   3. hpc/run_array.sbatch across tasks 1..3 + resume skip check
#   4. hpc/run_cohort.sbatch (argjacent-cohort + summarize_run.py)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/argj_hpc_test_XXXXXX")"
trap 'rm -rf "$TEST_DIR"' EXIT

echo "== Local HPC pipeline smoke test in ${TEST_DIR} =="

# 1. Generate 3 small paired-read metagenome samples carrying 2 plasmids with varying abundance
mkdir -p "${TEST_DIR}/reads" "${TEST_DIR}/out"
python3 - "${TEST_DIR}" <<'PY'
import gzip, os, random, sys

root = sys.argv[1]
rng = random.Random(42)
def rseq(n): return "".join(rng.choice("ACGT") for _ in range(n))
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]

# 1 kb gene shared by two plasmids (P1 and P2) with distinct flanks
gene = rseq(1000)
p1 = rseq(3000) + gene + rseq(3000)
p2 = rseq(3000) + gene + rseq(3000)
with open(os.path.join(root, "genes.fa"), "w") as f:
    f.write(f">testARG_1\n{gene}\n")

# Vary depths across 4 samples (including a 50x/50x tie in S1!) so cohort Part A resolves the tie
abunds = [("S1", 50, 50), ("S2", 85, 25), ("S3", 25, 85), ("S4", 60, 0)]
for sid, d1, d2 in abunds:
    r1_p = os.path.join(root, "reads", f"{sid}_R1.fastq.gz")
    r2_p = os.path.join(root, "reads", f"{sid}_R2.fastq.gz")
    q = "I" * 150
    rid = 0
    with gzip.open(r1_p, "wt") as f1, gzip.open(r2_p, "wt") as f2:
        for seq, cov in ((p1, d1), (p2, d2)):
            if cov <= 0:
                continue
            n_pairs = int(len(seq) * cov / 300)
            for _ in range(n_pairs):
                ins = rng.randint(360, 440)
                st = rng.randint(0, len(seq) - ins)
                frag = seq[st:st + ins]
                if rng.random() < 0.5:
                    a, b = frag[:150], rc(frag[-150:])
                else:
                    a, b = rc(frag[-150:]), frag[:150]
                rid += 1
                f1.write(f"@r{rid}/1\n{a}\n+\n{q}\n")
                f2.write(f"@r{rid}/2\n{b}\n+\n{q}\n")
PY

# 2. Sample sheet creation & check
echo "-- [1/4] make_samplesheet.py"
python3 "${SCRIPT_DIR}/make_samplesheet.py" --fastq-dir "${TEST_DIR}/reads" -o "${TEST_DIR}/samples.tsv"

# 3. submit_pipeline.sh --dry-run
echo "-- [2/4] submit_pipeline.sh --dry-run"
bash "${SCRIPT_DIR}/submit_pipeline.sh" \
    --samples "${TEST_DIR}/samples.tsv" \
    --genes "${TEST_DIR}/genes.fa" \
    --out "${TEST_DIR}/out" \
    --cpus 4 --mem 8G --time 00:30:00 \
    --dry-run

# 4. Run Slurm array script directly for tasks 1..4 (and test resume skip + --resubmit-missing)
echo "-- [3/4] run_array.sbatch (tasks 1..4 + resume skip)"
for tid in 1 2 3 4; do
    SLURM_ARRAY_TASK_ID="$tid" \
    SLURM_CPUS_PER_TASK=4 \
    SLURM_MEM_PER_NODE=4096 \
    SLURM_TMPDIR="${TEST_DIR}/scratch_${tid}" \
    REPO_ROOT="$REPO_ROOT" \
    SAMPLES_TSV="${TEST_DIR}/samples.tsv" \
    GENES_FA="${TEST_DIR}/genes.fa" \
    OUT_DIR="${TEST_DIR}/out" \
    EXTRA_ARGS="-D 2500" \
    bash -c 'mkdir -p "$SLURM_TMPDIR" && bash "'"${SCRIPT_DIR}/run_array.sbatch"'"'
done

# Verify --resubmit-missing detects 0 missing when all 4 are done
bash "${SCRIPT_DIR}/submit_pipeline.sh" \
    --samples "${TEST_DIR}/samples.tsv" \
    --genes "${TEST_DIR}/genes.fa" \
    --out "${TEST_DIR}/out" \
    --resubmit-missing --dry-run | grep -q "already complete" || { echo "FAIL: --resubmit-missing check failed"; exit 1; }

# Verify resume skip works
OUT_SKIP="$(
    SLURM_ARRAY_TASK_ID=1 SLURM_CPUS_PER_TASK=4 SLURM_MEM_PER_NODE=4096 \
    REPO_ROOT="$REPO_ROOT" SAMPLES_TSV="${TEST_DIR}/samples.tsv" \
    GENES_FA="${TEST_DIR}/genes.fa" OUT_DIR="${TEST_DIR}/out" \
    bash "${SCRIPT_DIR}/run_array.sbatch"
)"
echo "$OUT_SKIP" | grep -q "Already completed" || { echo "FAIL: resume skip did not trigger"; exit 1; }

# 5. Run Slurm cohort & summary script directly
echo "-- [4/4] run_cohort.sbatch (argjacent-cohort + summarize_run.py)"
HPC_DIR="$SCRIPT_DIR" \
REPO_ROOT="$REPO_ROOT" \
SAMPLES_TSV="${TEST_DIR}/samples.tsv" \
GENES_FA="${TEST_DIR}/genes.fa" \
OUT_DIR="${TEST_DIR}/out" \
WRITE_COHORT_GFA=1 \
bash "${SCRIPT_DIR}/run_cohort.sbatch"

# Verify outputs
for f in \
    "${TEST_DIR}/out/summary/run_metrics.tsv" \
    "${TEST_DIR}/out/summary/all_samples.summary.tsv" \
    "${TEST_DIR}/out/summary/all_samples.contexts.tsv" \
    "${TEST_DIR}/out/summary/all_samples.pairs.tsv" \
    "${TEST_DIR}/out/cohort/cohort.cohort_pairs.tsv" \
    "${TEST_DIR}/out/cohort/cohort.cohort.json"; do
    [[ -s "$f" ]] || { echo "FAIL: missing or empty output $f"; exit 1; }
done

echo ""
echo "--- run_metrics.tsv ---"
cat "${TEST_DIR}/out/summary/run_metrics.tsv"
echo "--- cohort.cohort_pairs.tsv ---"
cat "${TEST_DIR}/out/cohort/cohort.cohort_pairs.tsv"
echo ""
echo "ALL HPC SMOKE CHECKS PASSED."
