# ARGjacent

**Find trustworthy DNA neighborhoods around antibiotic resistance genes (or any target genes) straight from raw short-read metagenomes, without assembling the whole sample.**

Knowing which resistance genes are in a sample is only half the story. To tell whether a gene sits on a mobile plasmid, inside a transposon, or safely on a chromosome, you need the flanking DNA on both sides of it. Full metagenome assemblers are slow, memory-heavy, and often break or stitch false chimeras whenever the same gene appears in multiple places or next to repeats.

ARGjacent skips whole-sample assembly. It pulls out only the reads near your target genes, reconstructs every genomic neighborhood separately, and tells you base by base how far the reads actually prove each flank (`certain_bp`). When it hits an unresolvable repeat, it stops cleanly instead of guessing.

![How ARGjacent works](docs/graphical-abstract.png)

---

## How It Works

1. **Start from raw reads and target genes.** Give ARGjacent paired short reads (trimmed `.fastq.gz`) and a FASTA of genes you care about (for example, resistance genes found in your sample).
2. **Pull only nearby reads.** In two fast streaming passes, ARGjacent builds a lightweight map of short DNA anchors around your target genes and keeps only read pairs within your chosen distance (default 10 kb on each side). Most unrelated reads are ignored, so memory stays low.
3. **Extend outward and prove every base.** Starting from each end of a gene, ARGjacent extends base by base and splits cleanly whenever two different contexts diverge (such as a plasmid copy and a chromosomal copy). It then lays the reads back onto every branch and checks read pairs and read depth to mark the exact distance that is genuinely supported by molecules carrying the gene (`certain_bp`).
4. **Link both sides and pool cohorts.** For each gene, ARGjacent pairs left and right flanks only when the evidence rules out chimeras, and packs all branches into a compact Bandage-ready graph (`.gfa`). If you have multiple samples, `argjacent-cohort` uses natural abundance differences across samples to pair tied copies without re-reading the FASTQ files.

*(For the full algorithms, statistical proofs, and benchmark tables, see the [validation report](benchmarks/v05_validation/RESULTS.md) and the forthcoming paper.)*

---

## Big-Picture Performance

Tested with frozen code against matched **PacBio HiFi** single-molecule long reads from the same DNA samples and official **CAMI II/III** benchmarks:

- **High accuracy where it counts:** **0.33% failed flanks** (`9 / 2,757`) at 1 kb or longer across 21 unseen real metagenomes, compared to **22.5%** for whole-metagenome assembly (MEGAHIT) and **50% or more** for assembly-graph extractors.
- **No chimeric two-sided pairs on real data:** **0 / 151 chimeric pairs** across 31 real long-read metagenomes (compared to ~17% for naive left-by-right pairing).
- **Fast and light on RAM:** ~2 to 6 minutes and ~4 to 6 GB RAM for a typical human gut metagenome on a laptop (~20x faster than full assembly).
- **Scales cleanly to clusters:** Includes ready-to-run Slurm array and cohort scripts in [`hpc/`](hpc/README.md).

---

## Install

Requires Python 3.8+, `numpy`, a C compiler (`gcc` or `clang`), and `zlib` headers (`xcode-select --install` on macOS; `apt install gcc zlib1g-dev` on Debian/Ubuntu).

**With `pip` (recommended):**

```bash
git clone https://github.com/amulyabaral/argjacent.git
cd argjacent
python3 -m pip install -U pip setuptools
python3 -m pip install .
argjacent-selftest
```

**Or keep everything inside this folder:**

```bash
sh install.sh
sh demo/run_demo.sh
```

---

## Quick Start

### Single sample

```bash
argjacent -g genes.fa -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o results/sample
```

Useful everyday flags: `-D 10000` (flank size per side), `-t 8` (CPU threads), `--max-mem 16` (RAM limit in GB), `-w 40` (uses ~40% less RAM on very deep samples). Run `argjacent -h` to see all options.

### Multiple samples (cohort pairing)

```bash
argjacent-cohort \
  -g genes.fa \
  -i results/s1.flanks.json results/s2.flanks.json results/s3.flanks.json \
  -o results/cohort
```

For cohorts where each sample has its own gene list, pass `--sample-sheet samples.tsv`.

### Large runs on a Slurm HPC

See [`hpc/README.md`](hpc/README.md) for a ready-to-run Slurm workflow (`make_samplesheet.py`, `submit_pipeline.sh`, local scratch staging, automatic cohort run, and merged QC summaries).

---

## Main Outputs

| File | What is in it |
|---|---|
| `PREFIX.summary.tsv` | One line per gene side: status, number of neighborhoods found, and longest proven flank (`certain_bp`). |
| `PREFIX.contexts.tsv` | One line per branch (`c1`, `c2`, ...): total length, proven length, read depth, why it stopped, and nearby genes. |
| `PREFIX.flanks.fa` | FASTA sequences of all reconstructed flanks. |
| `PREFIX.pairs.tsv` | Certified two-sided (`left + gene + right`) pairs. |
| `PREFIX.gfa` | Compact, lossless gene-centered graph ready to open in Bandage. |
| `PREFIX.flanks.json` | Complete machine-readable record (also used as input to `argjacent-cohort`). |

---

## Want the Details?

- **Expanded v0.5.0 validation (25 CAMI II/III + real HiFi samples):** [benchmarks/v05_validation/RESULTS.md](benchmarks/v05_validation/RESULTS.md) and [PROTOCOL.md](benchmarks/v05_validation/PROTOCOL.md)
- **Earlier held-out validation history:** [docs/heldout-validation.md](docs/heldout-validation.md)
- **Slurm HPC guide:** [hpc/README.md](hpc/README.md)
- **Release history:** [CHANGELOG.md](CHANGELOG.md)

Released under the [MIT License](LICENSE). See [CITATION.cff](CITATION.cff) for citation metadata.
