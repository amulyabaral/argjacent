#!/usr/bin/env bash
# End-to-end local smoke test for the Global Metagenome Benchmark HPC scripts
# (runs completely offline without Slurm, blastn, or minimap2).
#
# Tests:
#   1. Two synthetic samples with BOTH + and - strand assembly contigs (.fasta.gz and .fa.zst):
#      - Gene A: two co-existing plasmids (60x + 40x), assembly collapses to 1 contig -> argj_multi_vs_asm_single
#      - Gene B: clean left flank + 5-copy 1.3 kb repeat at 1.2 kb on right; assembly stitches chimeric
#                extension past repeat wall -> both_concordant_1kb (left, 100% id) & asm_past_repeat_cliff (right)
#   2. submit_global_benchmark.sh --dry-run
#   3. run_global_array.sbatch (tasks 1..2 + resume skip check)
#   4. run_global_post.sbatch (summarize_run.py + per-BioProject cohort + global comparison + FASTA export)
#   5. run_annotate.sbatch / annotate_flanks.py (UniVec reagent vector, ISfinder exact left/right distance,
#      co-localized ARG, ORF + mobileOG plasmid classification, MGE masking, cross-biome HGT clustering)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/argj_global_test_XXXXXX")"
trap 'rm -rf "$TEST_DIR"' EXIT

echo "== Local Global Metagenome Benchmark smoke test in ${TEST_DIR} =="
mkdir -p "${TEST_DIR}/reads" "${TEST_DIR}/asms" "${TEST_DIR}/out" "${TEST_DIR}/dbs"

# 1. Build synthetic paired reads + simulated pre-made assemblies (gzip and zstd) + annotation DBs + manifest
python3 - "${TEST_DIR}" <<'PY'
import gzip, os, random, shutil, subprocess, sys

root = sys.argv[1]
rng = random.Random(123)
COMP = str.maketrans("ACGT", "TGCA")
CODONS = [
    "GCT", "GCC", "GCA", "CGT", "CGC", "AAT", "AAC", "GAT", "GAC", "TGT", "TGC",
    "CAA", "CAG", "GAA", "GAG", "GGT", "GGC", "CAT", "CAC", "ATT", "ATC", "CTG",
    "AAA", "AAG", "TTT", "TTC", "CCT", "CCG", "TCT", "TCC", "ACT", "ACC", "TGG",
    "TAT", "TAC", "GTT", "GTC",
]
CODON_MAP = {
    "ATG": "M", "GCT": "A", "GCC": "A", "GCA": "A", "CGT": "R", "CGC": "R", "AAT": "N", "AAC": "N",
    "GAT": "D", "GAC": "D", "TGT": "C", "TGC": "C", "CAA": "Q", "CAG": "Q", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "CAT": "H", "CAC": "H", "ATT": "I", "ATC": "I", "CTG": "L", "AAA": "K",
    "AAG": "K", "TTT": "F", "TTC": "F", "CCT": "P", "CCG": "P", "TCT": "S", "TCC": "S", "ACT": "T",
    "ACC": "T", "TGG": "W", "TAT": "Y", "TAC": "Y", "GTT": "V", "GTC": "V",
}
def rseq(n): return "".join(rng.choice("ACGT") for _ in range(n))
def rc(s): return s.translate(COMP)[::-1]
def make_orf(n_aa=100):
    cods = ["ATG"] + [rng.choice(CODONS) for _ in range(n_aa - 1)]
    return "".join(cods) + "TAA", "".join(CODON_MAP[c] for c in cods)

# Synthetic reference elements for Layer 1-3 annotation verification
is26_seq = rseq(650)                    # IS26_IS6_IS6
coarg_seq = rseq(500)                   # sul1_1_U12338
vec_seq = rseq(550)                     # pMB1_ColE1_vector
tra_dna, tra_aa = make_orf(100)         # mobileOG_001|traI|transfer|relaxase

with open(os.path.join(root, "dbs", "univec.fa"), "w") as f:
    f.write(f">pMB1_ColE1_vector cloning vector\n{vec_seq}\n")
with open(os.path.join(root, "dbs", "isfinder.fa"), "w") as f:
    f.write(f">IS26_IS6_IS6 insertion sequence\n{is26_seq}\n")
with open(os.path.join(root, "dbs", "coarg.fa"), "w") as f:
    f.write(f">sul1_1_U12338 sulfonamide resistance\n{coarg_seq}\n")
with open(os.path.join(root, "dbs", "mobileog.faa"), "w") as f:
    f.write(f">mobileOG_001|traI|transfer|relaxase|Plasmid\n{tra_aa}\n")

with open(os.path.join(root, "dbs", "db_paths.env"), "w") as f:
    f.write(f'export UNIVEC_DB="{root}/dbs/univec.fa"\n')
    f.write(f'export ISFINDER_DB="{root}/dbs/isfinder.fa"\n')
    f.write(f'export ARG_DB="{root}/dbs/coarg.fa"\n')
    f.write(f'export MOBILEOG_DB="{root}/dbs/mobileog.faa"\n')

# Gene A: 1 kb gene present on TWO distinct plasmids:
#   - plasmid1: carries traI plasmid ORF + IS26 at exact distance 400 bp on LEFT and sul1 at exact distance 300 bp on RIGHT
#   - plasmid2: carries pMB1 cloning vector (reagent contamination!) — missed by collapsed assembly!
geneA = rseq(1000)
p1_L = rseq(1600) + tra_dna + rseq(447) + is26_seq + rseq(400)  # len=3400; IS26 is 400 bp upstream of geneA!
p1_R = rseq(300) + coarg_seq + rseq(2600)                       # len=3400; sul1 is 300 bp downstream of geneA!
p2_L = rseq(2450) + vec_seq + rseq(400)                         # len=3400; pMB1 vector on plasmid2!
p2_R = rseq(3400)
plasmid1 = p1_L + geneA + p1_R
plasmid2 = p2_L + geneA + p2_R

# Gene B: 1 kb gene with clean left flank and a 5-copy 900 bp repeat starting at +800 bp on right
# (repeat exits at 1700 bp -> ARGjacent certifies ~1670 bp; chimeric assembly extends 3200 bp past cliff!)
geneB = rseq(1000)
rep = rseq(900)
b_L = rseq(3400)
b_R_pre = rseq(800)
b_R_true_exit = rseq(2000)
b_R_chim_exit = rseq(2000)
genomeB = b_L + geneB + b_R_pre + rep + b_R_true_exit
# Four decoy repeat copies at 40x each -> total 200x inside repeat, forces ARGjacent to halt at repeat exit
decoys = [rseq(1500) + rep + rseq(1500) for _ in range(4)]

with open(os.path.join(root, "genes.fa"), "w") as f:
    f.write(f">geneA_multi\n{geneA}\n>geneB_repeat\n{geneB}\n")

def sim_reads(molecules, r1_path, r2_path):
    q = "I" * 150
    rid = 0
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for seq, cov in molecules:
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

samples = [
    ("G1", [(plasmid1, 65), (plasmid2, 45), (genomeB, 45)] + [(d, 45) for d in decoys], "human_gut", "metaSPAdes 3.15.3"),
    ("G2", [(plasmid1, 50), (plasmid2, 80), (genomeB, 45)] + [(d, 45) for d in decoys], "wastewater", "MEGAHIT 1.2.9"),
]

# Simulated pre-made assemblies:
#   - Gene A collapsed to 1 contig (only plasmid1, broken at 700 bp on right, placed on MINUS strand in G1!)
#   - Gene B assembled cleanly on left, but CHIMERICALLY extends 2,000 bp past the repeat wall on right!
asm_ctg_A = p1_L + geneA + p1_R[:700]
asm_ctg_B = b_L + geneB + b_R_pre + rep + b_R_chim_exit

manifest_rows = []
for idx, (sid, mols, biome, asmr) in enumerate(samples):
    r1 = os.path.join(root, "reads", f"{sid}_1.fastq.gz")
    r2 = os.path.join(root, "reads", f"{sid}_2.fastq.gz")
    sim_reads(mols, r1, r2)

    # Put G1 contigs on reverse-complement (-) strand in .fasta.gz, and G2 on (+) strand (test .fa.zst if zstd present)
    fa_txt = (
        f">ctgA_{sid}\n{rc(asm_ctg_A) if idx == 0 else asm_ctg_A}\n"
        f">ctgB_{sid}\n{rc(asm_ctg_B) if idx == 0 else asm_ctg_B}\n"
    )
    if idx == 1 and shutil.which("zstd"):
        asm_p = os.path.join(root, "asms", f"{sid}.contigs.fa.zst")
        p = subprocess.Popen(["zstd", "-q", "-o", asm_p], stdin=subprocess.PIPE, text=True)
        p.communicate(fa_txt)
        assert p.returncode == 0
    else:
        asm_p = os.path.join(root, "asms", f"{sid}.fasta.gz")
        with gzip.open(asm_p, "wt") as out:
            out.write(fa_txt)

    manifest_rows.append([
        sid, sid, f"ERZ99900{idx+1}", "PRJEB_TEST1", f"ERS99900{idx+1}",
        biome, f"{biome} metagenome", "", "Denmark", "2023-06-01", "2023-10-01",
        "ILLUMINA", "NovaSeq 6000", "20000", "3000000", "0.00", "0.00",
        r1, r2, asm_p, asmr, "", "",
    ])

cols = [
    "sample", "run_accession", "erz_accession", "study_accession", "sample_accession",
    "biome_group", "scientific_name", "environment_biome", "country",
    "collection_date", "first_public", "instrument_platform", "instrument_model",
    "read_count", "base_count", "est_gbp", "fastq_gb",
    "r1_url", "r2_url", "assembly_url", "assembler", "logan_url", "logan_mb",
]
with open(os.path.join(root, "manifest.tsv"), "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in manifest_rows:
        out.write("\t".join(r) + "\n")
PY

echo "-- [1/4] submit_global_benchmark.sh --dry-run (with auto-chained --db-env)"
bash "${SCRIPT_DIR}/submit_global_benchmark.sh" \
    --manifest "${TEST_DIR}/manifest.tsv" \
    --genes "${TEST_DIR}/genes.fa" \
    --db-env "${TEST_DIR}/dbs/db_paths.env" \
    --out "${TEST_DIR}/out" \
    --cpus 4 --mem 8G --time 00:30:00 \
    --dry-run | grep -q "run_annotate.sbatch" || { echo "FAIL: --db-env did not chain run_annotate.sbatch"; exit 1; }

echo "-- [2/4] run_global_array.sbatch (tasks 1..2 + resume skip)"
for tid in 1 2; do
    SLURM_ARRAY_TASK_ID="$tid" \
    SLURM_CPUS_PER_TASK=4 \
    SLURM_MEM_PER_NODE=4096 \
    SLURM_TMPDIR="${TEST_DIR}/scratch_${tid}" \
    REPO_ROOT="$REPO_ROOT" \
    HPC_DIR="$SCRIPT_DIR" \
    MANIFEST_TSV="${TEST_DIR}/manifest.tsv" \
    GENES_FA="${TEST_DIR}/genes.fa" \
    OUT_DIR="${TEST_DIR}/out" \
    FLANK_BP=3200 \
    SKIP_FASTP=1 \
    bash -c 'mkdir -p "$SLURM_TMPDIR" && bash "'"${SCRIPT_DIR}/run_global_array.sbatch"'"'
done

OUT_SKIP="$(
    SLURM_ARRAY_TASK_ID=1 SLURM_CPUS_PER_TASK=4 SLURM_MEM_PER_NODE=4096 \
    REPO_ROOT="$REPO_ROOT" HPC_DIR="$SCRIPT_DIR" \
    MANIFEST_TSV="${TEST_DIR}/manifest.tsv" GENES_FA="${TEST_DIR}/genes.fa" \
    OUT_DIR="${TEST_DIR}/out" bash "${SCRIPT_DIR}/run_global_array.sbatch"
)"
echo "$OUT_SKIP" | grep -q "Already completed" || { echo "FAIL: resume skip failed"; exit 1; }

echo "-- [3/4] run_global_post.sbatch (summary + per-BioProject cohort + global comparison)"
HPC_DIR="$SCRIPT_DIR" \
REPO_ROOT="$REPO_ROOT" \
MANIFEST_TSV="${TEST_DIR}/manifest.tsv" \
GENES_FA="${TEST_DIR}/genes.fa" \
OUT_DIR="${TEST_DIR}/out" \
COHORT_MIN_SAMPLES=2 \
bash "${SCRIPT_DIR}/run_global_post.sbatch"

echo "-- [4/4] run_annotate.sbatch (5-layer annotation + ARGjacent vs Assembly biological delta)"
HPC_DIR="$SCRIPT_DIR" \
OUT_DIR="${TEST_DIR}/out" \
MANIFEST_TSV="${TEST_DIR}/manifest.tsv" \
DB_ENV="${TEST_DIR}/dbs/db_paths.env" \
THREADS=4 \
ANNOTATE_EXTRA_ARGS="--allow-python-fallback" \
bash "${SCRIPT_DIR}/run_annotate.sbatch"

# Verify fail-fast on HPC when neither blastn nor minimap2 is present and --allow-python-fallback is omitted
if ! command -v blastn >/dev/null 2>&1 && ! command -v minimap2 >/dev/null 2>&1; then
    if python3 "${SCRIPT_DIR}/annotate_flanks.py" \
        --fasta "${TEST_DIR}/out/comparison/global_cmp_argjacent_certified.fa" \
        --univec-db "${TEST_DIR}/dbs/univec.fa" \
        -o "${TEST_DIR}/out/fail_test" 2>/dev/null; then
        echo "FAIL: annotate_flanks.py should fail fast without --allow-python-fallback when aligners are missing"
        exit 1
    fi
fi

# Validate expected comparison + annotation assertions!
python3 - "${TEST_DIR}/out" <<'PY'
import os, sys
out = sys.argv[1]
for rel in (
    "summary/run_metrics.tsv",
    "cohorts/PRJEB_TEST1/PRJEB_TEST1.cohort_pairs.tsv",
    "comparison/global_cmp_per_side.tsv",
    "comparison/global_cmp_summary.tsv",
    "comparison/global_cmp_argjacent_certified.fa",
    "comparison/global_cmp_assembly_flanks.fa",
    "annotation/argjacent_annotations.tsv",
    "annotation/assembly_annotations.tsv",
    "annotation/argjacent_hgt_clusters.tsv",
    "annotation/reagent_vector_audit.tsv",
    "annotation/annotation_comparison_per_side.tsv",
    "annotation/annotation_comparison_summary.tsv",
):
    p = os.path.join(out, rel)
    assert os.path.isfile(p) and os.path.getsize(p) > 0, f"Missing/empty output: {p}"

rows = []
with open(os.path.join(out, "comparison/global_cmp_per_side.tsv")) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        rows.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))

print("\n--- global_cmp_per_side.tsv ---")
for r in rows:
    print(
        f"  {r['sample']} {r['gene']:12s} {r['side']:5s} "
        f"argj_cert={r['argj_certain_bp']:>4s} (n1k={r['argj_n_1kb']}) "
        f"asm={r['asm_max_bp']:>4s} (n1k={r['asm_n_1kb']}) "
        f"ov_id={r['overlap_id_pct'] or '-':>6s}% past_cliff={r['asm_past_cliff_bp']:>4s} -> {r['category']}"
    )

by_key = {(r["sample"], r["gene"], r["side"]): r for r in rows}
for sid in ("G1", "G2"):
    assert by_key[(sid, "geneA_multi", "left")]["category"] == "argj_multi_vs_asm_single"
    assert float(by_key[(sid, "geneA_multi", "left")]["overlap_id_pct"]) >= 99.5
    assert by_key[(sid, "geneA_multi", "right")]["category"] == "argj_multi_vs_asm_single"
    assert by_key[(sid, "geneB_repeat", "left")]["category"] == "both_concordant_1kb"
    assert float(by_key[(sid, "geneB_repeat", "left")]["overlap_id_pct"]) >= 99.5
    assert by_key[(sid, "geneB_repeat", "right")]["category"] == "asm_past_repeat_cliff"
    assert int(by_key[(sid, "geneB_repeat", "right")]["asm_past_cliff_bp"]) >= 500

# Verify 5-layer annotations (reagent vector, exact left/right distances, plasmid ORF, MGE masking, HGT clusters)
def load_tsv(rel):
    with open(os.path.join(out, rel)) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]

a_ann = load_tsv("annotation/argjacent_annotations.tsv")
m_ann = load_tsv("annotation/assembly_annotations.tsv")
delta = load_tsv("annotation/annotation_comparison_per_side.tsv")
hgt_a = load_tsv("annotation/argjacent_hgt_clusters.tsv")

# 1. ARGjacent detects pMB1 reagent vector on plasmid2 left in BOTH G1 and G2; collapsed assembly misses it!
a_vecs = [r for r in a_ann if r["is_reagent_vector"] == "1"]
m_vecs = [r for r in m_ann if r["is_reagent_vector"] == "1"]
assert len(a_vecs) == 2, f"Expected 2 reagent vector flanks in ARGjacent, got {len(a_vecs)}"
assert len(m_vecs) == 0, f"Expected 0 reagent vector flanks in collapsed assembly, got {len(m_vecs)}"
assert all(r["element_class"] == "reagent_vector" for r in a_vecs)
# Reagent vectors must be strictly excluded from HGT clustering!
assert all(r["hgt_cluster_id"] == "" for r in a_vecs), "Reagent vector leaked into HGT clustering!"

# 2. Exact left-side IS26 distance (400 bp) & plasmid traI ORF & MGE masking on plasmid1 left!
p1_lefts = [r for r in a_ann if r["gene"] == "geneA_multi" and r["side"] == "left" and r["has_is"] == "1"]
assert len(p1_lefts) == 2, f"Expected IS26 on plasmid1 left in G1 & G2, got {len(p1_lefts)}"
for r in p1_lefts:
    assert abs(int(r["nearest_is_bp"]) - 400) <= 5, f"Wrong left IS distance: {r['nearest_is_bp']}"
    assert r["is_families"] == "IS6", f"Wrong IS family: {r['is_families']}"
    assert r["element_class"] == "plasmid", f"Expected plasmid class via traI, got {r['element_class']}"
    assert int(r["unmasked_bp"]) < int(r["flank_len"]) - 600, "MGE/IS masking did not mask IS26+traI!"

# 3. Exact right-side co-ARG (sul1) distance (300 bp) on plasmid1 right (missed by 700 bp assembly break!)
p1_rights = [r for r in a_ann if r["gene"] == "geneA_multi" and r["side"] == "right" and int(r["n_coargs"]) > 0]
assert len(p1_rights) == 2, f"Expected sul1 on plasmid1 right in G1 & G2, got {len(p1_rights)}"
for r in p1_rights:
    assert abs(int(r["nearest_coarg_bp"]) - 300) <= 5, f"Wrong right co-ARG distance: {r['nearest_coarg_bp']}"

# 4. Cross-biome HGT clustering detects plasmid1 shared between G1 (human_gut) and G2 (wastewater)!
xb_clusters = [h for h in hgt_a if h["is_cross_biome"] == "1"]
assert len(xb_clusters) >= 2, f"Expected >=2 cross-biome HGT clusters between human_gut & wastewater, got {len(xb_clusters)}"

# 5. Per-side delta table explicitly flags vector_argj_only + multi_class_argj_collapsed_in_asm!
d_by = {(r["sample"], r["gene"], r["side"]): r for r in delta}
for sid in ("G1", "G2"):
    fl = d_by[(sid, "geneA_multi", "left")]["delta_flags"]
    assert "vector_argj_only" in fl, f"Missing vector_argj_only flag: {fl}"
    assert "multi_class_argj_collapsed_in_asm" in fl, f"Missing multi_class collapse flag: {fl}"

print("\n--- annotation_comparison_per_side.tsv ---")
for r in delta:
    print(
        f"  {r['sample']} {r['gene']:12s} {r['side']:5s} "
        f"ARGj=[{r['argj_class']:22s} IS={r['argj_nearest_is_bp'] or '-':>3s} coARG={r['argj_max_coargs']}] "
        f"Asm=[{r['asm_class']:22s} IS={r['asm_nearest_is_bp'] or '-':>3s} coARG={r['asm_max_coargs']}] "
        f"-> {r['delta_flags']}"
    )

print("\nALL GLOBAL METAGENOME BENCHMARK + 5-LAYER ANNOTATION SMOKE CHECKS PASSED!")
PY
