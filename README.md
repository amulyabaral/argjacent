# ARGjacent

**Certified genomic contexts of antibiotic resistance genes (and arbitrary target genes) directly from paired short-read metagenomes — without whole-metagenome assembly.**

ARGjacent (`v0.5.0`) recovers the upstream and downstream flanking DNA (e.g., ±10 kb) around every copy of every target gene in a metagenome or isolate in **two streaming passes over the reads**. When a gene sits in multiple genetic neighborhoods within the same sample (different plasmids, chromosomes, strains, or mobile elements), ARGjacent outputs one candidate flank per neighborhood and **certifies base-by-base how far the reads themselves prove that molecules carrying the gene continue along that flank (`certain_bp`)** — stopping cleanly before unresolvable repeats instead of stitching chimeras.

Key capabilities in `v0.5.0`:
- **Assembly-free, memory-bounded flank reconstruction (`argjacent`):** C minimizer anchor-walk + local multi-branch proposer + dual strict/tolerant read verifier with explicit read-pair, own-read budget, and flow-conservation certificates.
- **Compacted gene-centered GFA v1.0 graphs (`PREFIX.gfa`):** lossless blunt-end (`0M`) DAGs split at exact certainty boundaries with Bandage colors (`gene`, `certain`, `end_only`, `candidate`) and local per-segment depths (`DP:f:`).
- **Certified two-sided `[left + gene + right]` pairing (`PREFIX.pairs.tsv`):** four strictly ordered evidence tiers (`span_pair`, `snv_phase`, `unbranched` / Kirchhoff flow partitioning, and post-fork `depth_match`) that abstain cleanly rather than emit chimeras.
- **Multi-sample cohort co-abundance deconvolution (`argjacent-cohort` / `argjacent cohort`):** resolves tied-depth and shared-arm two-sided flank pairs across $S \ge 2$ samples, and deconvolves post-repeat exit arms across long shared repeats (e.g., `IS26`, `Tn3`, $\ge 400$ bp) across $S \ge 4$ samples directly from `flanks.json` (`depth_bins` + `gene_depth`) or optionally `aw.reads.tsv`.

---

## Install

Requirements: Python 3.8+, `numpy >= 1.20`, a C compiler (`gcc` or `clang`), and `zlib` headers (`xcode-select --install` on macOS; `apt install gcc zlib1g-dev` on Debian/Ubuntu).

### With `pip` (recommended)

```bash
git clone https://github.com/amulyabaral/argjacent.git
cd argjacent
python -m pip install -U pip setuptools    # setuptools >= 61 required
python -m pip install .                    # compiles aw, awx, awv into the package
argjacent-selftest                         # ~15 s self-test; ends with "All checks passed."
```

### Without installing outside the repository folder

```bash
sh install.sh
sh demo/run_demo.sh
./bin/argjacent --help
./bin/argjacent-cohort --help
```

### Building a local Conda package

A Conda recipe (which builds the C binaries and runs `argjacent-selftest` plus the full 68-test `pytest` suite) is provided in `recipe/`:

```bash
conda build recipe -c conda-forge --python 3.12
```

---

## Usage

### Single-sample flank calling (`argjacent`)

```bash
argjacent -g genes.fa -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o results/sample
```

- **`genes.fa`**: FASTA of target genes to extend — for example, the reference alleles of resistance genes detected in the sample by a read mapper (`kma -ipe R1.fq.gz R2.fq.gz -o sample -t_db resfinder -1t1 -cge -ef -nf`). The sample's copy does not need to match the database sequence identically: the last 300 bp of each gene end are inferred from the sample's reads and used if $\ge 90\%$ identical to the input sequence. Header names (first whitespace-delimited token) must be unique printable ASCII.
- **Reads**: Paired 4-line FASTQ (plain or `.gz`), Illumina 2×150 bp or longer (adapter-trimmed, e.g. with `fastp`). Because reads are streamed twice, inputs must be seekable files rather than pipes.

| Option | Default | Description |
|---|---:|---|
| `-D`, `--flank` | `10000` | Target flank length per side (bp) |
| `-t`, `--threads` | `min(8, CPUs)` | Worker threads (output is deterministic and thread-count invariant) |
| `--max-mem GB` | `60%` RAM | Memory budget (also respects Linux cgroup / Slurm limits). If the sample anchor map exceeds the budget, it is partitioned transparently across extra read passes with identical output |
| `--max-paths N` | `16` | Maximum candidate branches followed per gene side |
| `--min-share F` | `0.002` | Minimum read share retained during the anchor walk (`0` follows all branches) |
| `--save-map`, `--map FILE` | — | Save the Pass-1 sample anchor map and reuse it for additional gene sets (1 read pass instead of 2) |
| `-w`, `-k`, `-b` | `20`, `31`, auto | Minimizer window (`-w 40` reduces map memory by ~40%), anchor $k$-mer size (`11–32`), and $\log_2$ Bloom filter bits |
| `--reads-out` | off | Also write recruited read pairs near target genes (`PREFIX.near_1.fq`, `PREFIX.near_2.fq`) |
| `--pair-bp BP` | `500` | Certified flank radius $W$ used for two-sided flank pairing |
| `--no-gfa`, `--no-pair` | off | Disable writing `PREFIX.gfa` or `PREFIX.pairs.tsv` |
| `--keep` | off | Retain intermediate files (`aw.reads.tsv`, candidates, read-pass stats) in `PREFIX.tmp/` |

### Multi-sample cohort deconvolution (`argjacent-cohort` / `argjacent cohort`)

When a gene occurs across multiple metagenomes ($S \ge 2$), cross-sample abundance variation resolves ambiguities that single-sample short reads cannot separate (equal-depth plasmid ties, shared flanking cassettes, and shared mobile elements such as `IS26`):

```bash
argjacent-cohort \
  -g genes.fa \
  -i results/s1.flanks.json results/s2.flanks.json results/s3.flanks.json results/s4.flanks.json \
  -o results/cohort \
  --gfa-dir results/cohort_gfa
```

Also callable as `argjacent cohort ...`. For cohorts with sample-specific gene FASTA files, pass `--sample-sheet samples.tsv` with tab-separated columns `sample`, `flanks_json` (and optional `genes_fa`, `reads_tsv`).

- **`PREFIX.cohort_pairs.tsv`**: two-sided `[left + gene + right]` pairs resolved across $S \ge 2$ samples via bidirectional weighted Lawson–Hanson NNLS, single-sample `1×1` anchor transfer, and `2×2` complement subtraction gated by leave-one-sample-out $\Delta\chi^2$.
- **`PREFIX.cohort_repeats.tsv`**: post-repeat exit arms resolved past $\ge 400$ bp shared repeats across $S \ge 4$ samples by regressing `gene_depth` on post-repeat `depth_bins` (runs directly on `flanks.json` with no `--keep` required, and automatically upgrades to `aw.reads.tsv` k-mer measurement if `--keep` directories exist).
- **`PREFIX.cohort.json`** & **`--gfa-dir`**: complete JSON summary plus updated per-sample GFA v1.0 graphs containing cohort-resolved `pair` and `cohort_ext` `P` paths.

---

## Output Files

| File | Description |
|---|---|
| `PREFIX.summary.tsv` | One row per `(gene, side)` (`left` = upstream, `right` = downstream): `status` (`complete`, `partial`, `ended`, `not_found`, `failed`), `contexts`, `certain_bp`, `longest_bp`, `contexts_certain_1kb`, `read_depth`, `gene_end` (`sample`/`database`), `certain_end_bp` |
| `PREFIX.contexts.tsv` | One row per candidate context (`c1`, `c2`, … ordered by certainty): `length_bp`, `certain_bp`, `read_depth`, `stop`, `certain_stop`, `other_genes` (co-localized genes from `genes.fa`), `gene_end`, `main` (`yes` for the best-supported branch at every fork), `certain_end_bp`, `shared_end_bp` |
| `PREFIX.flanks.fa` | FASTA of all candidate flanks (left flanks oriented upstream 5'→3', so their certain bases are the trailing bases; right flanks start at the gene's 3' end) |
| `PREFIX.flanks.json` | Complete machine-readable record including `joins`, `open_joins`, `bridged`, strict `gene_depth`, 100-bp binned `depth_bins`, and per-gene certified `pairs` |
| `PREFIX.pairs.tsv` | Certified single-sample two-sided `[left + gene + right]` pairings with evidence tier (`span_pair`, `snv_phase`, `unbranched`, `depth_match`) and score |
| `PREFIX.gfa` | Lossless compacted GFA v1.0 graph (`[left suffix tree] -> [GENE] -> [right prefix tree]`) with `TP:Z:`, `DP:f:`, and Bandage `CL:Z:` colors (`#d95f02` gene, `#1b9e77` certain, `#66a61e` end-only, `#b3b3b3` candidate) |
| `PREFIX.log` | Detailed timing, memory, and status log |

**Two certainty claims (`certain_bp` vs `certain_end_bp`):**
- `certain_bp` certifies that molecules carrying the gene's full last 300 bp carry this flank that far.
- When another locus joins inside the gene's last 100–300 bp (e.g., mosaic tetracycline genes or truncated gene fragments), `shared_end_bp` reports how many terminal bases are shared and `certain_end_bp` reports how far the flank is certain for *a gene ending like this one* (its last 100 bp, within $\le 10$ differences). Whenever `shared_end_bp == 0`, both claims are equal.

---

## Method Overview

1. **Streaming minimizer anchor walk (`src/aw.c`):** Pass 1 extracts canonical minimizer anchors ($k=31, w=20$), discards singletons via an auto-sized Bloom filter, and records anchor adjacency links in a compact graph (switching nodes with degree $>16$ to an open-addressed hash table to eliminate hub-anchor blowup). A bounded outward walk (`--min-share`) marks all anchors within `-D` bp of any target gene, and Pass 2 recruits read pairs containing any reached anchor.
2. **Variant-aware multi-branch proposer (`src/awx.c`):** Extends base-by-base from each gene end (using both the database end and the sample's inferred 300-bp gene end), collapsing small strain patches (`-P 1 -Q 4`) while forking cleanly at structural junctions up to `--max-paths` branches per side.
3. **Dual strict + tolerant read verification (`src/awv.c`, `argjacent/verify.py`):** Lays all recruited reads onto every candidate in a single fused strict + tolerant pass (`-Z -q`). A flank position is called **certain** only when molecules originating in the gene (or prior to the last incoming foreign locus join) demonstrably reach it via spanning read pairs, an own-read budget accounting for ±20% coverage unevenness (`UNEV = 0.20`), or read depth exceeding the upper bound of all incoming foreign loci (with up to 100-bp sustained bridging). The final 30 bp of verified sequence are withheld as a safety margin (`certain_cut_bp = 30`).
4. **Graph compaction, two-sided pairing, and cohort deconvolution (`gfa.py`, `pair.py`, `cohort.py`):** Compacts per-gene candidates into a round-trip verified GFA v1.0 DAG with local `depth_bins` coverage, links upstream and downstream trunks across four strict evidence tiers, and optionally pools cohorts via weighted Lawson–Hanson NNLS.

---

## Validation Summary

Every version of the verification rules was frozen prior to evaluation on unseen real metagenomes with matched single-molecule **PacBio HiFi** long reads from the same DNA samples as ground truth (see [docs/heldout-validation.md](docs/heldout-validation.md)).

### 1. Held-out real metagenomes (unseen Block F, 12 library-concordant samples across 3 projects)
Eleven Malawian children's gut metagenomes (~2×114–131 bp TELL-seq vs PacBio Revio HiFi) plus two 2×151 bp NovaSeq vs Sequel II/IIe HiFi digester/bioreactor samples:
- **Native sample DNA primary claims ($\ge 1$ kb):** **0 / 1,388 failed (0.00%, exact 95% CI 0.00–0.27%)** — including 0 / 480 at $\ge 2$ kb and 0 / 167 at $\ge 5$ kb.
- **Native end-anchored claims ($\ge 1$ kb):** **0 / 1,415 failed (0.00%, 95% CI 0.00–0.26%)**.
- *(Raw mechanical check before separating two sequence-verified library artifacts: 73 / 1,498 [4.87%], consisting of 34 claims on one mislabeled BioSample `SAMN42799965` refuted by both HiFi and ONT, and 39 claims reconstructing a 4,724-bp `pMB1`/`ColE1` *bla*TEM-116 reagent expression plasmid present at 175–1,279× in 11/11 TELL-seq libraries and 0/2 NovaSeq libraries.)*

### 2. Head-to-head comparison against assembly and targeted flank tools
Three real metagenomes with PacBio HiFi truth (`H4` unseen gut, `H1` gut, `B-A2` sludge):

| Tool | Failed claims $\ge 1$ kb | Total wall time (3 samples) | Peak RAM |
|---|---:|---:|---:|
| **ARGjacent (`v0.4` / `v0.5` rules)** | **0 / 834 (0.0%)** | **9.1 min** | **1.1–3.7 GB** |
| ARGextender | 14% | hours | high |
| MEGAHIT (whole-metagenome assembly) | 22% | 3 h 43 min | high |
| ARGContextProfiler / Sarand | 50–52% | — | — |

### 3. Two-sided pairing, GFA compaction, and cohort deconvolution (`v0.5.0`)
- **Lossless GFA v1.0 round-trip:** 19,075 / 19,075 candidate paths across 32 real metagenomes reconstructed byte-for-byte with **1.56× compaction** (120.25 → 76.93 Mbp).
- **Real PacBio HiFi two-sided pairing (22 metagenomes: 10 dev + 12 unseen held-out):** **0 / 91 chimeras (0.00%, 95% CI 0.00–3.97%)** at ±500 bp and **0 / 58 (0.00%, 95% CI 0.00–6.16%)** at ±1,000 bp, vs ~20% for Cartesian pairing, 14–45% for Main-to-Main, and 8–10% for unguarded greedy depth matching.
- **Controlled suites (±18% coverage waves):** 0 chimeras across `span_pair`, `snv_phase`, shared-promoter `depth_bins` post-fork pairing, 3-plasmid cohort pairing with shared arms and 45×/45× ties (100% recall vs 33% single-sample), and 900-bp `IS26` post-repeat exit deconvolution (100% certified $\ge 2.5$ kb past `IS26`, 0 false exits, and 0 claims emitted on an adversarial 1:1 collinear decoy trap).

### 4. Typical runtime and memory (Apple M4, `-t 8`, gzip FASTQ)

| Sample | Read pairs | Genes | Wall time | Peak RAM |
|---|---:|---:|---:|---:|
| *E. coli* isolate (30×, 2×150) | 0.49 M | 100 | 11 s | 0.5 GB |
| Dried sludge metagenome | 35.2 M | 18 | ~65 s | 4.5 GB |
| Dried sludge with `--max-mem 3` (3 map parts) | 35.2 M | 18 | ~169 s | 2.9 GB |
| Human gut metagenome (`--max-mem 7`) | 54.1 M | 45 | ~5.5 min | 5.4 GB |

---

## Repository Layout

- `argjacent/` — Python package (`cli.py`, `pipeline.py`, `verify.py`, `gfa.py`, `pair.py`, `cohort.py`, `selftest.py`).
- `src/` — C engines (`aw.c`, `awx.c`, `awv.c`).
- `bin/` — Source-folder CLI wrappers (`argjacent`, `argjacent-cohort`).
- `tests/` — 68-test automated suite (`pytest`).
- `demo/` — Quick self-test runner (`demo/run_demo.sh`).
- `recipe/` — Conda recipe (`meta.yaml`, `conda_build_config.yaml`).
- `docs/` — Validation documentation (`heldout-validation.md`).

---

## License & Citation

Released under the [MIT License](LICENSE). See [CITATION.cff](CITATION.cff) for citation metadata.
