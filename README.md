# ARGjacent

**Certified genomic contexts of antibiotic resistance genes (or arbitrary target genes) directly from paired short-read metagenomes, without whole-metagenome assembly.**

ARGjacent reconstructs the upstream and downstream flanking DNA (default 10 kb per side) around every copy of every target gene directly from raw paired-end short reads. When a gene occurs in multiple genomic neighborhoods within the same sample (distinct plasmids, chromosomes, strains, or mobile elements), ARGjacent separates each structural branch and certifies base by base how far the reads themselves prove that molecules carrying the gene continue along that flank (`certain_bp`), stopping cleanly before unresolvable repeats rather than stitching chimeras.

![ARGjacent mechanistic workflow](docs/graphical-abstract.png)

---

## How It Works

1. **Minimizer anchor graph and directed outward walk (`aw`, streaming passes):**
   A saturating Bloom pre-pass (`pass0`) identifies canonical minimizers (`k=21, w=11` CLI default, roughly one anchor per 6 bp) seen at least twice across the sample before link construction (`pass1`), discarding singleton sequencing errors while preserving `2x` low-coverage links and recording directed within-read edges plus dead-end `FR` mate bridges (`-i 420`). Starting from anchors inside each target gene, ARGjacent walks strictly outward (no backtracking) up to distance `D` (default 10 kb per side). At every branch, the gene's read share splits in proportion to link counts, and branches carrying less than `--min-share` (default `0.2%`) are pruned so repeats stay bounded. Pass 2 streams the reads once more and recruits every read pair in which at least one mate hits a reached anchor.
2. **Variant-aware local extension and sub-21-bp stall rescue (`awx`):**
   Starting from each gene end (using both the database reference end and the sample's own inferred terminal 300 bp, guarded against terminal paralog hijacks), ARGjacent extends outward one base at a time from reads overlapping by at least 31 bp (with multi-tier `FR` mate, sub-21-bp [`6-30 bp`] mate-stepping bridge, 1-bp-step backward `K=31`, and 16-mer tip rescue when coverage drops to `2.5-10x`). Single-base SNPs, small indels, and short variant patches that reconverge within 31 to 60 bp remain on the same path, and reads differing only at known variant sites stay active (up to 4 sites per read) so they can phase downstream junctions. Sustained structural divergence (plasmid vs chromosome junctions, mobile element boundaries) forks into separate candidate arms (up to `--max-paths`, default 16 per side), explored in order of read share.
3. **Dual strict and tolerant verification (`awv` + `verify.py`):**
   All recruited read pairs are mapped back onto every candidate in a single fused strict (exact) and tolerant (SNP/small-indel aware) pass. A position along a flank is certified (`certain_bp`) only if molecules originating in the gene demonstrably continue past all incoming foreign-locus joins, verified by spanning paired-end templates, an own-read budget allowing 20% coverage unevenness, read depth exceeding incoming foreign joins, or unbranched mate-bridged low-depth continuity (`Rule U`, `<= 20x` with zero foreign joins and single-copy depth ceilings), withholding a 30 bp safety margin at the cutoff. When another locus joins inside the gene's final 100 to 300 bp, `certain_end_bp` separately reports how far the flank is certified for genes sharing that terminal 100 bp.
4. **GFA compaction, two-sided pairing, and cohort deconvolution (`gfa.py`, `pair.py`, `cohort.py`):**
   Shared candidate prefixes and suffixes are compacted into a round-trip lossless blunt-end (`0M`) GFA v1.0 graph split at exact certainty boundaries (`gene`, `certain`, `end_only`, `candidate`) with 100-bp depth bins (`depth_bins`). Left and right flanks are linked (`[left + gene + right]`) across four strictly ordered tiers: direct read-pair spans across short genes (`<= 650 bp`), phased internal gene SNPs, unbranched (`1x1`) or flow-conserved (`1xk`) splits, and post-fork depth matching requiring 3-sigma separation from competing pairs (abstaining otherwise). Across multiple samples (`argjacent-cohort`), bidirectional weighted non-negative least squares (NNLS) with leave-one-sample-out `Delta chi^2` gating resolves tied-depth pairs (`Part A`, `S >= 2`) and explores post-repeat exits across `>= 400 bp` repeats (`Part B`, `S >= 4`) directly from `flanks.json`.

---

## Validation Summary

Evaluated with frozen `v0.5.0` code against matched single-molecule **PacBio HiFi** long reads from the same DNA samples and official **CAMI II/III** benchmarks:

- **Unseen real PacBio HiFi metagenomes (21 samples):** **9 / 2,757 (`0.33%`, 95% CI `0.15-0.62%`)** primary flanks `>= 1 kb` failed strict long-read checks (`0 / 1,388` in Block F native DNA; `9 / 1,369` across 9 additional real gut/mock metagenomes, of which 8 were single 72-80 bp tandem-repeat indels and 1 was a 74-bp structural overshoot).
- **Head-to-head on 3 real PacBio HiFi metagenomes (`H1`, `B-A2`, `H4`):** **ARGjacent `0 / 834 (0.0%)`** in 9.1 min, vs **ARGextender `11 / 79 (13.9%)`** in 18.2 min, **MEGAHIT `31 / 138 (22.5%)`** in 3 h 43 min, and **ARGContextProfiler / Sarand `50.5-52.0%`**.
- **CAMI II/III benchmarks (16 samples: Strain Madness, CAMI III Human Gut, Marine):** `1 / 384` unexplained `>= 1 kb` fault groups (`0.26%`) and `0 / 292` failures at `>= 5 kb` on Strain Madness.
- **Two-sided pairing (31 real PacBio HiFi metagenomes):** **0 / 151 (`0.00%`, 95% CI `0.00-2.41%`)** chimeras at `±500 bp` and **0 / 105 (`0.00%`)** at `±1,000 bp`, vs `16.8%` for Cartesian pairing. Across 18 multi-sample cohort samples, `argjacent-cohort` Part A added **0 / 29** chimeras at `±500 bp` while doubling multi-context recall on Strain Madness (`7.6% -> 14.1%`). Part B post-repeat extension had `0 / 108` false exits on isolated synthetic `IS26` repeats but `23.6%` failures (`35 / 148`) in chained Strain Madness repeats and should be treated as exploratory in dense repeat chains.
- **Lossless GFA v1.0 compaction (51 metagenomes, 29,720 paths):** `1.60x` sequence compaction with `0` round-trip mismatches.
- **Ultra-low-coverage recovery (`v0.6.0`) & 86-metagenome global ENA `metaSPAdes` benchmark (`561.2 Gbp`, 5 biomes, 40 BioProjects):** `0 / 288` adversarial `2-6x` chimera-trap faults; `2.67x` more true `>= 1 kb` flanks than `MEGAHIT` on *E. coli* 536 at `2.5x` (`16` vs `6`, `0` errors, `98x` faster) and `8.0x` more `< 12x` certified flanks in paired real metagenomes (`2 -> 16`). Across all 86 global metagenomes (`median 12.6 s`, `2.51 GB` RAM), `ARGjacent v0.6.0` recovered **743 certified flanks `>= 1 kb` (`699` sides; `1,636` candidate sides `>= 1 kb`)** vs **622 flanks (`608` sides)** in pre-computed ENA `metaSPAdes` assemblies (`+19.5%` flanks, **3.45x more multi-context ARGs [`38` vs `11`]**, **1.63x more IS linkages [`147` vs `90`]**, **1.53x more co-ARG cassettes [`185` vs `121`]**, **98.93% [`370/374`] overlap identity concordance**, `13` cohort-certified post-repeat exits in `0.4 s`, `105 / 743 [14.1%]` `UniVec` vector backbones flagged, and `7` vector-free intercontinental human-livestock HGT clusters).

Full validation protocols, per-sample tables, and forensic checks are in [benchmarks/v05_validation/RESULTS.md](benchmarks/v05_validation/RESULTS.md) and [docs/heldout-validation.md](docs/heldout-validation.md).

---

## Install

Requires Python 3.8+, `numpy >= 1.20`, a C compiler (`gcc` or `clang`), and `zlib` headers (`xcode-select --install` on macOS; `apt install gcc zlib1g-dev` on Debian/Ubuntu).

```bash
git clone https://github.com/amulyabaral/argjacent.git
cd argjacent
python3 -m pip install -U pip setuptools
python3 -m pip install .
argjacent-selftest
```

To build locally inside the repository folder without pip:

```bash
sh install.sh
sh demo/run_demo.sh
```

A Conda build recipe is also provided in `recipe/`.

---

## Usage

### Single sample

```bash
argjacent -g genes.fa -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o results/sample
```

Inputs: paired 4-line FASTQ (plain or `.gz`, adapter-trimmed) and a FASTA of target genes. Passing a full database (e.g., ResFinder or CARD) works directly; on deep metagenomes, pre-screening detected alleles first (for example with KMA) keeps local extension faster and leaner. Key options: `-D` flank length per side (default `10000`), `-t` threads, `--max-mem` memory limit in GB, `-k` anchor k-mer size (default `21`), `-w` minimizer window (default `11`), `--max-paths` candidate branches per side (default `16`), `--pair-bp` pairing radius (default `500`), `--save-map` / `--map` to reuse the Pass-1 anchor graph across gene sets. Run `argjacent -h` for all options.

### Multi-sample cohort deconvolution

```bash
argjacent-cohort \
  -g genes.fa \
  -i results/s1.flanks.json results/s2.flanks.json results/s3.flanks.json \
  -o results/cohort \
  --gfa-dir results/cohort_gfa \
  -t 8
```

You can also invoke this as `argjacent cohort ...`. For cohorts with sample-specific gene FASTA files, pass `--sample-sheet samples.tsv` (columns: `sample`, `flanks_json`, optional `genes_fa`, `reads_tsv`). Gene FASTA headers must match the keys inside `flanks.json`; `--family-mode` (default `auto`) collapses ResFinder allele copies (`<allele>_<copy>_<acc>`) without colliding RefSeq, CARD, or custom gene names (`exact` disables grouping).

### Large runs on a Slurm HPC

Ready-to-run Slurm array and dependent cohort scripts live in [`hpc/`](hpc/README.md) (`make_samplesheet.py`, `submit_pipeline.sh`, `run_array.sbatch`, `run_cohort.sbatch`, `summarize_run.py`).

---

## Outputs

| File | Description |
|---|---|
| `PREFIX.summary.tsv` | One row per `(gene, side)`: `status`, candidate count, `certain_bp`, `longest_bp`, `read_depth`, `gene_end`, and `certain_end_bp`. |
| `PREFIX.contexts.tsv` | One row per candidate branch (`c1`, `c2`, ...): `length_bp`, `certain_bp`, `read_depth`, `stop`, `certain_stop`, `other_genes`, `main`, `certain_end_bp`, `shared_end_bp`. |
| `PREFIX.flanks.fa` | FASTA of all candidate flanking sequences. |
| `PREFIX.pairs.tsv` | Certified single-sample two-sided `[left + gene + right]` pairings with evidence tier and score. |
| `PREFIX.gfa` | Lossless compacted gene-centered GFA v1.0 graph with `TP:Z:`, `DP:f:`, and Bandage `CL:Z:` colors. |
| `PREFIX.flanks.json` | Full machine-readable record including `joins`, `depth_bins`, `gene_depth`, and `pairs` (input to `argjacent-cohort`). |

Example `PREFIX.summary.tsv` rows (full 10 kb flank certifies to `9970` bp after the 30 bp safety margin; a repeat 3 kb downstream stops `certain_bp` at the repeat while retaining all 4 candidate branches in `PREFIX.flanks.fa` and `PREFIX.gfa`, which opens directly in Bandage with certainty colors applied):

```text
gene                side   status    contexts  certain_bp  longest_bp  contexts_certain_1kb  read_depth  gene_end  certain_end_bp
blaCTX-M-15_1_...   left   complete  1         9970        10000       1                     42          sample    9970
blaCTX-M-15_1_...   right  partial   4         4268        10000       4                     41          sample    4268
```

---

## License and Citation

Released under the [MIT License](LICENSE). See [CITATION.cff](CITATION.cff) for citation metadata and [CHANGELOG.md](CHANGELOG.md) for release notes.
