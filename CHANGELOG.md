# Changelog

## 0.6.0 (2026-11-27)

Low-coverage sensitivity upgrades (`2.5x-10x`) across the minimizer scout (`aw`), exact read extender (`awx`), and flank verifier (`verify.py`), while preserving `< 0.5%` primary certified chimera rates and multi-copy repeat-wall abstention.

New
- **Two-stage saturating Bloom pre-pass (`src/aw.c`):** `pass0` records `>= 2x` minimizer sightings before freeing `bloom1` prior to `pass1`, eliminating first-read link loss on `2x` anchors without increasing peak graph RAM; gene-seed minimizers are pre-inserted into `bloom` before `pass0`.
- **Dead-end paired-end `FR` mate bridges & adapter-hub guard (`src/aw.c`):** records `FR` mate-pair links (`-i 380` default) flagged with bit 14 (`Edge.mate`) and traverses them strictly at within-read dead ends (`tot == 0`) so exact within-read Dijkstra coordinates remain untouched; `is_chimeric_hub` (`deg_sv >= 24 && max_count <= 3` on entry side `sv`) blocks random-offset adapter fan-ins without clipping high-depth plasmid or multi-copy IS backbones.
- **Stall-only `FR` mate & `K2 = 19` read-tip rescue (`src/awx.c`):** when exact extension stalls (`tot3 < 1`, `-R 1` default), rescues concordant `FR` mates (`ov = 14..150 bp`, mismatch-tiered with exact terminal 8 bp) and `19..43 bp` exact overlaps from a compact `K2 = 19` read-tip index (`<= 32` bucket hits, `>= 80%` forward 10-bp consensus), plus isolated low-Q (`Q < 15`) pass-through and homopolymer/dinucleotide/Illumina adapter tail filtering (`has_junk_tail`).
- **Rule U low-coverage mate-bridged certification (`argjacent/verify.py`):** certifies low-depth (`d0_u <= 20x`) single-context positions `x >= g` backed by `>= 2` concordant spanning read pairs (`UMARG = 25`, bridging `<= 90 bp` 1-pair Poisson necks when `pspan >= 2` sustains for `>= 180 bp` beyond) while enforcing **zero** non-allele entering reads, **zero** non-patch leaving reads, `FaF == 0`, `Lseen == 0`, no multi-copy depth jump, and no ragged-end pileup.

Fixed
- **Terminal gene-end paralog hijack guard (`argjacent/pipeline.py`, `argjacent/cli.py`):** rejects `redo = 300` sample-end replacements when mismatches crowd near the gene boundary (`> 15` in the terminal 100 bp or `> 9` in the terminal 30 bp), preventing C/N-terminal paralog or co-occurring allele hijacks while retaining uniform `~90%`-identity orthologs.

## 0.5.1 (2026-11-25)

Cohort scaling optimizations (PR #1) and targeted bug fixes from the 25-sample expanded validation and second-round code audit.

Fixed
- **Enforced short-gene ceiling on `span_pair` (`argjacent/pair.py`):** strictly gates both accumulation and emission of direct read-pair bridging (`span_pair`) to `len(gene) <= 650 bp`. Previously, multi-trunk genes `> 650 bp` entering `index_sample_reads` for `snv_phase` could still emit `span_pair` across a co-existing cassette-deletion variant (the sole 2-sided chimera observed on CAMI II Strain Madness `s30 aadA2`).
- **Preserved read-backed single-sample pairs in `argjacent-cohort` (`argjacent/cohort.py`):** `cohort_pair_family` now restores single-sample `span_pair` and `snv_phase` pairs stored in `flanks.json["genes"][g]["pairs"]` whenever `parameters.pair_bp >= W`, keeping them in `cohort_pairs.tsv`, `cohort.json`, and `--gfa-dir` outputs and allowing them to act as `cohort_anchor` seeds across samples.
- **Safe cross-database gene-family grouping (`argjacent/cohort.py`):** replaced unconditional `g.split("_")[0]` grouping with `--family-mode {auto,exact,resfinder}` (default `auto`), which collapses ResFinder allele copies (`<allele>_<copy>_<accession>`) while keeping RefSeq (`NC_...`, `NZ_...`, `WP_...`), CARD, Prokka, and custom FASTA names distinct.
- **Outlier-resistant Part B repeat fork selection (`argjacent/cohort.py`):** replaced `min(fork_votes)` with per-sample `±50 bp` modal cluster selection so a single outlier sample in a large cohort cannot drag `fork_bp` inside a shared repeat for all other samples.
- **Collision-free GFA segment IDs and accurate read-length `RC:i:` (`argjacent/gfa.py`, `argjacent/cli.py`):** disambiguates distinct gene names that sanitize to the same GFA ID (such as `bla-1` and `bla_1`) and scales `RC:i:` by actual sample `read_len` instead of hardcoded `150`.
- **Aligned C binary version strings (`src/{aw,awx,awv}.c`)** to `0.5.1` and removed hardcoded local machine paths in `benchmarks/v05_validation/run_*.py`.

Improved
- **Exact multi-sample cohort scaling (`argjacent/cohort.py`, `argjacent/pair.py`):** `co_present` gating before NNLS, sorted leave-one-sample-out early-break, sparse column/row pruning with Gram-Cholesky + SVD fallback in `_nnls_weighted`, NumPy `uint8` vectorized `seq_sim_fast` / `divergence_bp` / `_ExactSeqClusters`, cached FASTA loading in `_load_cohort_inputs`, slimmed per-sample GFA worker payloads, and `-t/--threads` `spawn` multiprocessing (`5.5-21x` faster with exact equivalence; `--cohort-cpus` added to `hpc/submit_pipeline.sh`).
- **CI:** added GitHub Actions workflow (`.github/workflows/ci.yml`) testing Linux (GCC) and macOS (Clang) across Python 3.10 and 3.12 (`pip install .`, `argjacent-selftest`, `pytest`, and `hpc/test_hpc_local.sh`).

## 0.5.0 (2026-11-22) — Renamed to ARGjacent

Renamed package and CLI from `bigargs-anchorwalk` (`aw-flank`, `aw-cohort`, `aw-flank-selftest`) to **`argjacent`**
(`argjacent`, `argjacent-cohort` / `argjacent cohort`, `argjacent-selftest`) under the MIT License
(`https://github.com/amulyabaral/argjacent`).

Adds lossless gene-centered GFA v1.0 graph export (`PREFIX.gfa`), 100-bp binned exact depth profiles (`depth_bins` and
`gene_depth`), certified single-sample two-sided `[left + gene + right]` flank pairing (`PREFIX.pairs.tsv`), and
multi-sample cohort co-abundance deconvolution (`argjacent-cohort` / `argjacent cohort`) for both left–right pairing and
post-repeat exit resolution across shared mobile elements (`IS26`, `Tn3`, etc.).

**100 % backwards compatible with 0.4.0 / 0.4.1:** the three C binaries (`aw`, `awx`, `awv`), all verification rules and
thresholds in `verify.py`, all existing CLI flags, and all columns/contents of `PREFIX.summary.tsv`,
`PREFIX.contexts.tsv`, and `PREFIX.flanks.fa` are unchanged.

New
- **Per-context local depth profile (`verify.py`, `flanks.json`):** every verified candidate now records `gene_depth`
  (median strict exact-match depth over the last 200 bp of the gene) and `depth_bins` (100-bp binned strict exact-match
  median depth along the candidate flank, always preserved from the strict pass even when tolerant verification extends
  `verified`). Enables local per-segment GFA depths, post-fork exclusive depths for shared-promoter pairing, and
  JSON-only post-repeat cohort deconvolution without `--keep`.
- **Compacted gene-centered GFA v1.0 export (`argjacent/gfa.py`, `PREFIX.gfa`, opt-out `--no-gfa`):** lossless blunt-end
  (`0M`) DAG `[left suffix tree] -> [GENE] -> [right prefix tree]` split at exact `certain_bp` and `certain_end_bp`
  boundaries with Bandage `CL:Z:` colors (`gene` orange, `certain` teal, `end_only` green, `candidate` grey), `TP:Z:`
  tier tags, local `DP:f:` depths from `depth_bins`, and round-trip verified `P` records for every left/right candidate
  and certified two-sided pair. Verified 100 % lossless on 19,075 / 19,075 paths across 32 real metagenomes (120.25 →
  76.93 Mbp, 1.56× compaction).
- **Certified single-sample two-sided flank pairing (`argjacent/pair.py`, `PREFIX.pairs.tsv`, `flanks.json["genes"][g]["pairs"]`,
  opt-out `--no-pair`, `--pair-bp 500`):** links upstream and downstream flanks across a gene through four strictly
  ordered evidence tiers: (1) `span_pair` direct read-pair bridging across short genes (`<= 650 bp`); (2) `snv_phase`
  internal gene-SNV haplotype chaining; (3) `unbranched` (`1 × 1` and Kirchhoff flow-conserved `1 × k` / `k × 1` with
  pruned-fork and distal-depth safeguards); and (4) `depth_match` (symmetric log-normal + Poisson delta-method flow
  bijection with `UNEV = 0.18`, `z <= 1.8`, and permutation `Delta_chi2 >= 9.0`, using `depth_bins` post-fork exclusive
  depths whenever branches share a proximal cassette/promoter).
- **Multi-sample cohort deconvolution (`argjacent/cohort.py`, `argjacent-cohort` command and `argjacent cohort` subcommand):**
  - *Part A (`PREFIX.cohort_pairs.tsv`):* resolves tied-depth and shared-arm two-sided flank pairs across `S >= 2`
    samples via bidirectional weighted Lawson–Hanson NNLS, single-sample `1 × 1` cohort-anchor transfer, and `2 × 2`
    complement subtraction gated by leave-one-sample-out `Delta_chi2`.
  - *Part B (`PREFIX.cohort_repeats.tsv`):* resolves post-repeat exit arms beyond shared repeats (`>= 400 bp`, e.g.
    `IS26`) across `S >= 4` samples by fitting `g_trunk ~ Y_exit @ w` (`w >= 0`) with completeness (`chi2_red <= 2.2`),
    stoichiometry (`0.65 <= w <= 1.40`), `Delta_chi2 >= 12.0`, and leave-one-sample-out `Delta_chi2 >= 6.0` guards.
    Works directly from 0.5.0 `flanks.json` (`depth_bins` + `gene_depth`, tolerating proposer dropouts cleanly) and
    automatically upgrades to two-pass `aw.reads.tsv` k-mer measurement when `--keep` temp directories are present.
  - Optional `--gfa-dir DIR` writes updated per-sample `.gfa` graphs incorporating cohort-resolved `pair` and
    `cohort_ext` `P` paths.

Validation (`benchmarks/dev-tools/experiments/gpc/RESULTS.md`)
- **Real PacBio HiFi two-sided pairing across 22 real metagenomes** (10 Beijing human gut dev + 12 frozen unseen
  held-out metagenomes): **0 / 91 = 0.00 % chimeras (95 % CI 0.00–3.97 %)** at ±500 bp (`0 / 79` single-sample, `0 / 12`
  cohort-added) and **0 / 58 = 0.00 % (95 % CI 0.00–6.16 %)** at ±1,000 bp, vs ~20 % chimeras for Cartesian pairing,
  14–45 % for Main-to-Main, and 8–10 % for unguarded greedy depth matching.
- **Controlled benchmarks (with ±18 % coverage waves):** 0 chimeras across `span_pair`, `snv_phase`, shared-promoter
  `depth_bins` post-fork pairing, 3-plasmid cohort pairing with shared arms and equal-depth ties (100 % recall vs 33 %
  single-sample), and 900-bp `IS26` post-repeat exit deconvolution (100 % certified ≥ 2.5 kb past `IS26`, 0 false exits,
  and 0 claims emitted on an adversarial 1:1 collinear decoy trap).

Tests & packaging
- 68 tests (`pytest tests`, adding `tests/test_gfa_pair_cohort.py` covering GFA round-trip and local `DP:f:`,
  `depth_bins` shared-promoter unlocking, pruned-fork chimera guard, `span_pair`, equal-depth abstention, 3-plasmid
  cohort pairing, JSON-only `IS26` repeat deconvolution without `--keep`, and 1:1 collinear trap refusal) and 6
  self-test cases (`argjacent-selftest`). Added `./bin/argjacent`, `./bin/argjacent-cohort`, `LICENSE` (MIT), and
  `CITATION.cff`. Public repository initialized at `https://github.com/amulyabaral/argjacent`.

## 0.4.1 (2026-11-22)

Pure speed and pathology fixes across the three C programs (`aw`, `awx`, `awv`) and the Python verifier (`verify.py`,
`pipeline.py`, `cli.py`). **No rule, threshold, or tie-breaker changed:** `.summary.tsv`, `.contexts.tsv`, `.flanks.fa`,
`.map`, and `.flanks.json["genes"]` are byte-for-byte identical to 0.4.0 across the test suite, *E. coli* 30× (919
candidates), human gut metagenome E-G8 (2.06 M kept pairs, 1,041 candidates), and hub-anchor stress tests.

Speed and pathology fixes
- **`aw` (`src/aw.c`): open-addressed per-node hash table (`hub_ht`) for high-degree hub anchors (`deg > 16`).** On
  Block F's NovaSeq bioreactor sample `F-J2` (`SRR16526566`), chimeric poly-C + adapter 31-mers created anchor nodes
  with tens of thousands of distinct neighbour edges, spending 1,601.5 s (99 % of the run) scanning singly-linked lists
  in `edge_half()`. Nodes with at most 16 half-edges keep the compact linked list (0 extra bytes); once a node exceeds
  16 half-edges, lookups switch via bit 31 of `Node.head` to an open-addressed hash table (`HubEnt`, 8 bytes/slot,
  counted in `--max-mem` and freed before Pass 2). On a 30,000-degree hub-anchor benchmark (175k pairs), the read pass
  drops from 4.8 s to 0.1 s (48× faster; 14× end-to-end) and writes the identical `.map` byte for byte.
- **`awv` (`src/awv.c`) + `pipeline.py` + `cli.py`: fused single-pass strict + tolerant record generation (`-Z`), `-q`
  fast array mode, packed 21-mer read pre-filter, and multi-allele candidate deduplication.** 0.4.0 invoked `awv` twice
  over `aw.reads.tsv` (once strict, once tolerant), unpacked every read's bases/qualities/reverse-complements before
  k-mer lookup, computed `low3` quality scans and nine prototype `.i32` arrays that `verify.py` does not read, and
  re-evaluated identical candidate sets shared by close alleles of one gene family. Now one `awv` invocation loads
  `aw.reads.tsv` once and writes both strict (`vs_`) and tolerant (`vt_`) record files in one parallel worker pass,
  checks candidate 21-mers directly on packed 2-bit `RS` words before unpacking non-matching reads, skips unused arrays
  under `-q` (writing zero-length `.idx` spans so the 20-slot binary record layout stays compatible with `verify_ref.py`),
  and deduplicates identical `(glen, cands)` verification jobs across multi-allele genes.
- **`verify.py`: single-load binary record reader (`_load_records`) and `O(L)` window filters.** Binary record files
  (`.idx`, `.i32`, `.ent`, `.lev`, `.tpl`) are parsed and stably sorted once per gene instead of twice when a side has
  an inside-gene join (`shared_end_bp > 0`); the 41-bp rolling median uses `int32` `np.partition` instead of `float64`
  `np.median`, and the 67-bp rolling maximum uses 7-step 1-D doubling `np.maximum` instead of a `67 × L` temporary
  matrix.
- **`awx` (`src/awx.c`): rolling 31-mer index scan, multi-threaded bucket sort, and duplicate 400-bp gene-end task
  deduplication.** Index construction scans 31-mers in `O(L)` per read and sorts the `2^24` prefix buckets across `-t`
  threads; identical 400-bp gene ends shared by multiple ResFinder alleles are extended once and copied to each allele.

Measured wall-clock speedup (Apple M4, `-t 8`, identical Python 3.12 + NumPy 2.5 environment, 0 output differences)
- *E. coli* 536 at 30× (0.49 M pairs, 100 genes, 919 candidates): verify step **13.2 s -> 7.9 s** (1.67× faster), total
  run **16.0 s -> 10.9 s** (1.47× faster).
- Human gut metagenome E-G8 kept reads (2.06 M pairs, 39 genes, 1,041 candidates): verify step **49.6 s -> 30.0 s**
  (1.65× faster), total run **63.1 s -> 42.9 s** (1.47× faster).
- Chimeric hub-anchor stress case (175k pairs, degree 30,000): read pass **4.8 s -> 0.1 s** (48× faster), total run
  **5.07 s -> 0.36 s** (14.0× faster), `.map` byte-identical.

Tests
- 62 tests (`pytest tests`, adding `test_hub_anchors_and_duplicate_gene_ends_give_identical_results`) and 6 self-test
  cases (`aw-flank-selftest`). Conda packages built and tested for `osx-arm64` across Python 3.10–3.13 (`dist/conda/`).

## 0.4.0 (2026-11-22)

What 0.3.0's held-out test, adversarial made-up cases (even and uneven depth), and sixteen known-answer genome cases
showed, what was changed for each, and the held-out result of the frozen 0.4.0 package on block F. The method is the
same: propose, then verify.

**Held out as a package on block F** (`benchmarks/heldout/PROTOCOL.md`, addendum 8 and its result;
`benchmarks/dev-tools/experiments/pkg04/NOTE`): thirteen completely unseen samples from three new projects (eleven
Malawian children's gut metagenomes with variable-length ~2×114–131 bp TELL-seq short reads and Revio HiFi, plus an
acid-phase digester sludge and a stillage bioreactor with 2×151 bp NovaSeq and Sequel II/IIe HiFi), run once alongside
paired frozen 0.3.0. On native sample DNA across all twelve library-concordant samples, **0 of 1,388 primary claims
certain for 1 kb or more failed (0.00 %, exact 95 % interval 0.00 to 0.27 %)** and **0 of 1,415 end-anchored claims
failed (0.00 to 0.26 %)**, with 0 of 167 at 5 kb and 139 / 62 / 22 gene-family sides confirmed at 1 / 2 / 5 kb (paired
0.3.0: 3 of 1,702 native claims ≥ 1 kb failed at two *erm*(F) soft-join forks on `F-M7` and `F-M9`, where 0.4.0 stopped
cleanly before the fork; shorter native claims 100–999 bp: 34 / 2,800 = 1.2 % vs 90 / 2,881 = 3.1 %). Before separating
two sequence-verified library-level artifacts (`F-M3`'s mislabeled `SRR29981012` short-read deposit, refuted by both
Revio HiFi and PromethION ONT of the same BioSample, and a 4,724-bp circular `pMB1`/`ColE1` *bla*TEM-116 reagent vector
present in 11 / 11 TELL-seq libraries and 0 / 2 NovaSeq libraries), the raw mechanical count across all thirteen samples
is 73 / 1,498 (4.87 %; 0.3.0: 92 / 1,852, 4.97 %).

Found wrong in 0.3.0 (and in the 19 Nov development rules)
- **The flow rule trusted the count of reads seen coming on** (block E: all 22 failed claims of 3,016, three gene
  sides). A locus whose reads differ from the path in single bases before the place is under-counted there. Now:
  what came on is at least what the depth has risen; the count is also taken from the reads that reach far before the
  place; and the flow rule can find no more of the gene's own reads than the budget allows.
- **Reads that leave the path in a single base were not charged to the budget.** Where the gene's molecules go off
  along a relative 7 to 12 % away while another locus holds the depth, most of their reads leave like that (made-up
  cases: 33 wrong claims in 37,440, most at 240× and more). Now they are charged where reads are also seen leaving
  nearby and another locus is on the path (`SOFTW = 300`, `SOFTF = 0.25`); 1 of 37,440 is left. Closing it altogether
  (a check of the budget by the depth) took a fifth to a half of the claims on real samples, every one of them borne
  out by long reads, and was not done.
- **Read pairs "from the gene" included another locus's.** Where a locus shares the gene's last bases, its reads
  begin in the gene too and their mates lie along its own flank. Now such pairs must begin before any place in the
  gene where another locus was seen coming on, with a 60-bp margin (`GINSIDE`, `GMARGIN = 60`).
- **Joins inside the gene were looked for in its last 100 bases only**, while the path is read off the sample over
  the last 300. A candidate that went over to another locus further in was judged as if that locus's reads were the
  gene's. Now the whole 300 bases are looked at (`INSIDE = 300`), and the record program marks what each read is like
  before the place where it comes on, so that the sample's other allele of the gene (within 15 %, also shifted) is
  told from another locus (foreign there).
- **The budget's bar** did not allow for the chance in what had been charged to it (`BSIG = 3.0`), nor for **uneven
  read depth along a molecule**: comparing a depth crest at the gene with a trough at a fork produced 35 wrong claims
  in 25,740 made-up cases at 40 % depth swings (6 at 20 %). Now the budget must stand clear of 20 % of what was
  charged (`UNEV = 0.20`, combined block-wise over 300-bp blocks `UBLOCK = 300`), and under the tolerant measure the
  budget is first taken down by `2 * 0.20` of all that other loci brought on (`UPOW = 2.0`). An extra cap allowance on
  the tolerant flow rule (`FCAPU`) was tested and left switched off (`False`): it closed 2 wrong claims in 25,740
  adversarial cases (20 % swing with a 4× relative at the gene end), but cut the rarer context in the 80×/50× self-test
  from 10 kb to 5.9 kb and cost 1.6 % of certified bases (4 % of 2-kb molecules) on two gut samples.
- **Charges before the gene end (19 Nov fault, fixed 20 Nov; `AFTERG = True`).** Counting reads that leave in the
  gene's last 20 bases was harmless in 0.3.0, but inflated `BSIG` and `UNEV`, crushed `FCAP`, and disabled the
  tolerant depth charge whenever relatives of the gene left right before the gene end. Only reads leaving after the
  gene end are now charged.

New
- **A second, tolerant way of counting reads** (`awv -T 50 -B 1 -I 12 -g 30 -w 9`), beside the strict one. A sample
  often holds several versions of a stretch, strains that differ in a base here and a few bases there; counted
  strictly, the other versions' reads count for nothing and are seen leaving and coming on at every difference, and
  certainty ended after a few hundred bases. Counted tolerantly they stay on the path. The tolerant count has
  safeguards of its own (its budget is checked against the depth; among a crowd of other loci it does not vouch; a
  locus that came on stays supply for the rest of the candidate; template "from before" uses a full read-length
  margin), because under it a relative's reads can come and go unseen. A position is certain if the evidence holds
  under either count.
- **The proposer keeps reads through places where versions differ** (`awx -P 1 -Q 4`). 0.3.0 dropped a read at its
  first disagreement with the path, so a fork that only the less common version takes was seen by few of its reads
  or none; and a gene side's sixteen candidates could be versions of one context.
- **Bridging short un-vouched stretches** (`BRIDGE = 100`, `SUSTAIN = 300`, `BRQ = 0.2`): a stretch of up to 100 bp
  that nothing vouches for (for example, a short repeat's reads lying on the path briefly) does not end certainty if
  reads carry the path across it, at most 20 % of the gene's own reads were charged in it, and the evidence holds again
  beyond it for 300 bp without a break.
- **End-anchored second claim (`certain_end_bp`, `shared_end_bp`).** Where another sequence comes onto the gene's last
  100 to 300 bases, `certain_bp` stops (it cannot tell which sequence owns the flank), while `certain_end_bp` gives
  how far the flank is certain as a context of *a gene that ends like this one* (its last 100 bases, provided those
  100 bases differ from the given gene in at most 10 places, `NEARDIFF = 10`). Never less than `certain_bp`.
- **Known-answer cases** as a standing check, besides held-out real samples: made-up cases built to provoke wrong
  claims with even and uneven depth (`synth_sound5.py`), and reads simulated from real genomes with the answer read
  off the genome (`refcase.py`, `refcase_w.py`, `battery_t.sh`, `battery_w.sh`): single *E. coli*, two-strain
  mixtures, *Klebsiella* with six resistance plasmids, *Shigella* with its hundreds of IS copies, and four uneven-depth
  variants (20 % and 40 % swings). The self-test has six cases (adding `strains`, `alleles`, and `shared`).

Development measurements (these samples and cases shaped the rules; NOT held out)
- **Known-answer battery, sixteen cases** (twelve even depth, `battery_t6.log`; four uneven depth at 20 % and 40 %
  swings, `battery_w2.log`): **0 wrong claims** in either the primary claim or the end-anchored claim (0.3.0 had 3
  wrong on *E. coli* 536 + CFT073). Gene sides of 200 with a context certain for 1 / 5 / 10 kb: one *E. coli* genome
  192 / 167 / 142 (30×) and 195 / 182 / 166 (60×); synthetic two-strain mixtures 167–188 / 90–128 / 37–69 (0.3.0:
  98–156 / 15–28 / 1–3); *E. coli* K-12 + CFT073 144 / 32 / 5 (0.3.0: 63 / 6 / 2); 536 + CFT073 163 / 98 / 46 (0.3.0:
  115 / 26 / 7); K-12 + *Shigella* 166 / 54 / 15 (0.3.0: 76 / 12 / 6); *Shigella* alone 170 / 126 / 94 (0.3.0:
  175 / 128 / 93); *Klebsiella* + 6 plasmids 23 / 14 / 6 of 24 (`kp1`) and 30 / 20 / 9 of 30 (`kp2`).
- **Seventeen real metagenomes with kept reads, tolerant proposer** (what the installed package runs; `devq17.log`,
  pooled via `toleval.py`):
  - Primary claim (`uF`): **1,882 claims ≥ 1 kb, 0 failed (0.00 %)**, 762 ≥ 2 kb, 99 ≥ 5 kb, **4.22 Mb**; gene-family
    sides borne out at 1 / 2 / 5 kb: 163 / 89 / 32; long reads certain at 1 / 2 / 5 / 8 kb: 42 / 26 / 7 / 4 % (874
    molecules at 8 kb); shorter claims 100–999 bp: 134 of 4,367 checkable failed (3.1 %).
  - End-anchored claim (`uFE`): **1,935 claims ≥ 1 kb, 0 failed (0.00 %)**, 804 ≥ 2 kb, 103 ≥ 5 kb, **4.42 Mb**;
    gene-family sides borne out at 1 / 2 / 5 kb: 169 / 94 / 34; long reads certain at 1 / 2 / 5 / 8 kb: 42 / 27 / 7 /
    4 % (896 at 8 kb).
  - Version 0.3.0 on the same seventeen samples: 1,953 claims ≥ 1 kb, **22 failed (1.13 %)**, 4.82 Mb, 166 / 102 / 40
    gene-family sides; a candidate of the tolerant proposer matches 77 / 64 / 39 / 22 % of long reads at 1 / 2 / 5 / 8 kb
    (0.3.0's proposer: 73 / 61 / 36 / 20 %).
- **Fifteen real metagenomes with the standard proposer** (isolating the rules from the proposer; `devq18.log` +
  `devq18b.log`, all fifteen including E-G1): **0 failed** in every new variant — final form (`uOz20`) 1,849 claims
  ≥ 1 kb (90 ≥ 5 kb, 4.12 Mb), end-anchored (`uEOz20`) 1,890 claims ≥ 1 kb (85 ≥ 5 kb, 4.25 Mb).
- **Adversarial made-up cases** (25,740 per setting; `sndo*`, `sndg*`, `sndk*`, `sndj*`): 0 wrong at even depth (with
  or without a 4× relative at the gene end), 1 at a 20 % swing, 2 at a 20 % swing with a 4× relative at the gene end,
  12 at a 40 % swing.
- **Installed package against development pipeline** (`regress5.py` on E-G8, E-G7, E-G2; 2,263 candidates): 0 of 2,263
  differ on the primary claim in strict or tolerant (4,526 comparisons), and 0 differ on `verified_end = max(claim,
  end)`.

Changed for users
- `certain_bp` keeps its meaning, now strictly anchored to the gene's last 300 bases. A context's certain part may run
  over places where the sample's strains differ (the sequence shown is one version) and over bridged stretches up to
  100 bp.
- **New output fields for the end-anchored claim and bridging:**
  - `summary.tsv`: new column `certain_end_bp` (never less than `certain_bp`).
  - `contexts.tsv`: new columns `certain_end_bp` and `shared_end_bp`.
  - `flanks.fa`: adds `certain_end=N` in the header where the end-anchored claim is longer than `certain_bp`.
  - `flanks.json`: adds `end_len`, `shared_end`, `end_stop`, and `bridged` per context.
- Under 0.3.0 the copy of a gene counted as the sample's own if its last 300 bases were within 10 % of the sequence
  given. That still holds for where a candidate starts. Reads of another copy up to about 15 % away before an inside
  join now count as the gene's other allele, so such a copy can lend its flank to the gene.
- A run takes somewhat longer because reads are laid on the candidates twice (strict and tolerant): an *E. coli*
  isolate at 30× takes 25 s and 0.6 GB peak (0.3.0: 9 s, 1.1 GB).

Tests
- 61 tests (`pytest tests`: adds tests for the tolerant measure, uneven-depth guard, bridging, charges-after-gene-end,
  and end-anchored claim) and 6 self-test cases (`aw-flank-selftest`: `clean`, `repeat`, `mixed`, `strains`, `alleles`,
  `shared`). All pass on the installed build (`rebuild_t8.log`).

## 0.3.0 (2026-11-19)

A new method for the part after the reads have been chosen: **propose, then verify**. 0.2.1 walked a map of anchors
and vouched for its own paths; held out on 34 real metagenomes with long-read truth, 3.2 % of its claims "certain for
1 kb or more" failed (`docs/heldout-validation.md`). Four prototypes of the new method were frozen and held out in
turn; the fourth failed in 5 of 1,585 claims (0.3 %) on fifteen fresh samples, with a confirmed claim on a quarter
more gene sides at 1 kb and nearly twice as many at 2 and 5 kb. 0.3.0 is that prototype with a repair for each of the
two places where it failed, as a package.

**Held out as a package** on sixteen further samples the same day (protocol, addendum 7 and its result): 22 of 3,016
claims certain for 1 kb or more failed, 0.73 % (95 % interval 0.46 to 1.10 %); none of 220 claims of 5 kb or more;
a confirmed claim of 1 / 2 / 5 kb on 359 / 227 / 87 gene sides. Under the bar of 1 %, not settled under it. The 22
are three gene sides on two samples and one fault of the flow rule (a locus that comes on softly is under-counted);
a repair is in development and not in this version. The five samples of earlier blocks on which 0.2.1's read pass
never fit in memory ran in 51 to 297 s and 4.0 to 6.8 GB and gave 389 claims, none failed.

How it works now
- **Read pass** (`aw`, as before): a thin map of the sample, a walk outward from every gene, and the read pairs near
  the genes. Its result is the same as that of 0.2.1, read for read.
- **Propose** (`awx`, new): from each gene end the path is extended base by base with the kept reads that match its
  end without a difference; where they disagree for more than a base or a small insertion, it forks, and up to 16
  arms per gene side are followed, the best-supported first. This is done twice: starting from the gene as given, and
  with the gene's last 300 bases taken from the sample's own reads. Candidates of the second kind are used if those
  300 bases are at least 90 % identical to the gene as given; else those of the first kind.
- **Verify** (`awv` and `anchorwalk/verify.py`, new): every kept read is laid on every candidate of a gene, and for
  each position it is recorded how many reads match, which reads come onto the path or leave it there, and what
  their mates do. A position is certain only if molecules that carry the gene demonstrably reach it: templates from
  the gene or from before the last place where another locus came on, a budget of the gene's own reads that have not
  been seen leaving, or more reads than everything that came on can supply. Certainty also ends at a place that no
  read crosses while reads pile up at its edges. The last 30 bases of what is vouched for are not called certain.
- The caller of 0.2 (`anchorwalk/contexts.py`) is gone, and with it what it alone reported.

Changed for users
- `certain_bp` means what it meant; how it is arrived at is new, and on real data it is right far more often and
  reaches further (numbers above). A whole flank of 10,000 bp that is certain is given as 9,970.
- **No `exact_bp`.** 0.2.1 marked a part of each context as exact, base for base. The new candidates are built from
  reads that match without a difference, but no claim about single bases was held out for them, so none is made.
- `summary.tsv`: columns `exact_bp`, `partial_gene_copies`, `reads_per_gene_anchor` and `branches_not_followed` are
  gone; `read_depth` and `gene_end` (`sample` or `database`: where the candidates of that side start from) are new.
  Status `ambiguous` is now `partial` (certain up to a point, candidates beyond).
- `contexts.tsv`: `exact_bp` and `events` are gone; `certain_stop` says why the certain part ends
  (`end_of_candidate`, `another_locus_on_the_path`, `own_reads_spent`, `pairs_do_not_cross`, `reads_do_not_cross`,
  `too_few_reads`), `gene_end` as above. `stop` is why the candidate itself ends (`reached_target`,
  `ran_out_of_reads`, `loop`, `work_limit`).
- `flanks.json`: per context `id`, `flank`, `flank_len`, `context_len`, `depth`, `stop`, `certain_stop`, `main`,
  `gene_end`, `gene_end_diffs`, `joins`, `open_joins`, `genes`.
- `--strict` and `--min-link` belonged to the old caller; they are accepted and do nothing. `--max-paths` is the
  number of candidates followed per gene side (default 16, as before).
- Gene names may now hold any printable ASCII character, `|` included (names of CARD); inside the run a gene is
  known by its number.
- Reads shorter than 140 bp: the rules' read-length terms are scaled (not held out); under 70 bp nothing is called
  certain, and the run says so.

Memory and time
- **`--max-mem` now bounds the read pass.** The map of the sample was the part that did not fit: five of 41 real
  samples could not be run in 9.5 GB. If the map does not fit in the memory given (default: 60 % of the machine's),
  it is built in parts, each with one more pass over the reads, and written to a temporary file that the walk reads;
  the kept reads are the same. On a sludge sample of 35 million pairs: 4.5 GB and 65 s without a limit, 2.9 GB and
  169 s with `--max-mem 3` (0.2.1: 6.6 GB, 75 s).
- The map itself takes about a quarter less (links and anchors in fixed pieces that are never copied, a denser table,
  no per-anchor walk tables), and nothing of it is held during the second read pass.
- The saved map (`--save-map`) has a new format; maps saved by 0.2 are refused with a message.
- After the read pass a sample takes seconds to about a minute: the proposer and the record program are C, threaded
  over gene ends and genes, and hold the kept reads in memory once (about 400 bytes per kept read at the peak).

Tests
- 57 tests: the three C programs alone (threads, gzip, cut-off files, a map in parts and under a memory limit, saved maps),
  the command (outputs, errors, repeatability), the merge of the two proposer runs, and the rules against the
  prototype's rule file (`benchmarks/dev-tools/rules-reference/`). The tests of the old caller went with it.

## 0.2.1 (2026-11-14)

The first real metagenome: a human gut sample, with long reads of the same DNA as the check. Of the 116 contexts that
0.2.0 called certain for 1 kb or more, the long reads contradict 3, and one shorter claim: two places, two faults of
the method. A third fault had hidden two wrong claims of the second kind on the simulated data. With 0.2.1 none of
the 496 claims on that sample (certain for 100 bp or more) is contradicted. The same day the test machine
ran out of memory, which led to the memory changes below. Numbers, method and the comparison with ARGextender:
`docs/benchmarks.md`.

Fixed
- **Read depth as a witness at forks was measured from the wrong level.** A branch that carries more reads than
  everything that can be on the path from other loci must hold reads of the gene's own loci. 0.2.0 took "what can
  be on the path" as the rise of the depth above the depth next to the gene. A path's own depth drops along the way
  (it leaves forks; where strains differ it carries only its own strain's reads), and another locus fills the gap
  without the depth rising: on the gut sample both ways out of a repeat passed, one of them another locus's
  (erm(B)). Now the rise is measured from the lowest level the path has had since its last fork, up to the highest
  level that lasted for a read length (or the level at the fork itself), and what was on board at the last fork is
  carried over.
- **A branch with few reads was vouched for whenever no read had been seen coming on.** A locus with up to about 9
  reads per link can come onto a path without one read showing it. Such a branch is now vouched for only by reads
  that reach back to the gene itself (tet(40) on the gut sample; two contexts of CAMI III sample 15).
- **A candidate with a short or empty flank cut the certain and exact parts of its siblings**, down to nothing for
  an empty one. On sample 15 this hid the two wrong claims above, and cost 31 exact contexts.
- **Another locus that comes onto the path right next to the gene hid itself.** Where most reads on a path are
  another locus's, nothing is called exact ("crowded"). That was judged against the median depth of the first eight
  links, and noticed at the fourteenth link at the earliest: a locus coming on within those links raised the median,
  and the bases up to where it was noticed were called exact. Found on CAMI III sample 8, which no version had seen,
  run once with the code otherwise final: five contexts of one gene side with a base of the other locus 88 bp from
  the gene. Now the gene's own level is not taken higher than its first four links or its own anchors show when the
  eight carry three times as much, and exactness ends where the crowd began, not where it was noticed.
- `flanks.json` lists the genes in the order of the input file (it was the order in which the calls finished).

Changed
- 16 candidate paths per gene side instead of 8 (`--max-paths`), and they are spent differently: when a path has
  ended, the waiting branch whose fork lies nearest to the gene is walked next (it was the last one pushed, so the
  far forks of the strongest branch used the paths up). On the gut sample 158 instead of 115 of the 294 contexts
  that the long reads show are found, with no claim more that they contradict; the list of contexts is 1.6 to 1.75
  times as long and the calling takes longer (CAMI III sample 15: 83 s in all instead of about 60; on the gut sample the packed read tables below make up for it). Where a path parts into branches of which none has
  15 % of the reads, the part they share is listed as a context of its own (stop reason `forks`); without that,
  more paths lost exact contexts.
- `--strict` now means: vouch for fork branches on reads alone, never on read depth. (It used to add a rule that the
  new handling of thin branches makes the default.)
- `contexts.tsv` has a last column `main`, and `flanks.fa` headers `main=yes`, for the one context per gene side that
  took the best-supported branch at every fork. Up to `certain_bp` nothing changes; beyond it this is the candidate to
  take if one must be taken.

Memory
- **The reads of a gene are held packed.** Calling a gene held every read's anchors as Python objects, about 2.4 kB
  per read, and the columns of a whole pile-up until its end: up to 1.9 GB for one gene of the gut sample, and `-t 8`
  calls eight genes at a time (the eight largest together: 13 GB). Now the anchors are in arrays and a pile-up is
  added up every 4,096 reads: 0.67 GB for that gene, 4.9 GB for the whole run. The results are the same, gene
  by gene, on every data set of `docs/benchmarks.md` (451 genes).
- **`--max-mem GB`**: what the genes called at the same time may use together. Default: half of the machine's
  memory (on Linux: of what a container or batch job is allowed, if that is less). Genes are started largest first,
  and only while the estimated memory of those being called fits; a gene alone always runs. The first pass over the
  reads is not covered (it takes what the map of the sample takes; `-w 40` makes that smaller).

Benchmarks
- `benchmarks/real_metagenome/`: a real sample checked against long reads, without storing the long-read run.
- `benchmarks/argextender/`: ARGextender (the flank step of ARGprofiler) on the same reads.
- All eleven CAMI III samples run again; a twelfth (sample 8), never seen by any version, run once with the final
  code. It showed the last fault of the list above; with that fixed, the twelve simulated metagenomes have
  4 flagged contexts among 4,486, at two places, the same as with 0.2.0.
- The scorer of the simulated samples padded true left contexts the wrong way round (no number of sample 15 changes).
- The download of real reads is verified (`benchmarks/real/run.sh`): three of three read files fetched from the
  archive on the day arrived damaged, and fastp passes a damaged file on, cut short, with exit status 0.

Tests: the faults as small simulated cases that fail with the old rules (`AW_FLOWREF=0 AW_THINBR=0`); the depth rule
on the real series that made its first form too cautious; the memory limit, the packed tables against the reads file,
and pile-ups added up in blocks (`tests/test_memory.py`).

## 0.2.0 (2026-11-14)

First version tested on metagenomes with known source genomes (CAMI II and CAMI III toy data sets) and on real
reads (two E. coli K-12 isolates). 0.1 was correct on the simulated metagenomes but nearly useless: half of the gene
sides had under 30 trusted bases, and the walk reached millions of anchors (6.7 GB and about half an hour for one
sample). Numbers and method: `docs/benchmarks.md`.

Results
- One context per place of the gene. A gene in three genomes gives three contexts, each with its own lengths,
  instead of one "trusted" length that ends at the first difference between them.
- Two claims per context: `certain_bp` (a real context of the gene) and `exact_bp` (base for base).
- New output `PREFIX.contexts.tsv`; `summary.tsv` has new columns; `flanks.fa` headers changed
  (`len= certain= exact= depth= stop=`). The JSON keeps `trusted_len` (= exact) and adds `context_len` (= certain),
  events with reasons, the other input genes on each context, and partial copies of the gene.
- Statuses are now `complete`, `ambiguous`, `ended`, `not_found`, `failed`.

Method
- Bounded walk: a branch is followed only if it carries at least 0.2 % of the gene's reads still on their way
  (`--min-share`). The walk no longer floods repeats.
- Joins (reads of another locus coming onto the path) are found from the reads; a fork is vouched for by reads that
  reach back before the join, or by read counts that the other loci cannot supply.
- Rules for low coverage, found on samples held back from development: a single entering read bounds what its locus
  can supply; a dead end with real reads next to a join is ambiguous; a gene or branch with few reads is not followed
  through a repeat longer than a read; exact bases need three agreeing reads; a branch with under 15 % of the reads
  has no exact bases before its fork.
- Reads on a link that skips an anchor are counted with the way they take.
- From the last held-back round (four samples, run after everything else was frozen; two places in one of them
  were wrong): a fork of a gene or branch with few reads is a context fork only if reads of the branch reach back to
  the gene, or the branch carries more reads than an unseen locus could, or the branches together carry no more
  reads than the gene has (at about 10x a second copy of a repeat comes onto the path without any read showing it
  about once in ten times). And where the read depth is well above the depth next to the gene, a second base is
  looked for with every anchored base of a read, also under its anchors on a variant arm the path did not take
  (there the other loci can be the majority, and their reads hardly counted at the variant itself).

Real reads (first tests: E. coli K-12, Genome Analyzer IIx 2×100 and MiSeq 2×300)
- k-mers that are a run of one base or of a 2- or 3-base unit are no anchors (junk reads made paths loop for hours).
- The walk writes every link of every reached anchor; before, a path that passed a repeat the other way round than
  the walk had could see one exit where the reads show several, and call it certain.
- Reads are laid onto a path by the path's own anchors only (adapter read-through no longer ends up in a flank); a
  branch that ends within 200 bp of its fork is not vouched for; a path stops at the third visit of an anchor.
- Inside a repeat longer than a read, a branch that looks like another locus's makes the path ambiguous instead of
  being dropped silently.
- The read length is taken from the reads (was fixed at 150).
- Genes share one file of reads through an index each, instead of a copy each (disk use about halved, more for isolates).

Speed
- Threads in the read passes (`-t`), with output identical for any number of threads: both passes over a
  16.6 M pair metagenome take about 30 s (100 s with one thread).
- In total such a sample takes about one minute and 2.7 GB; 0.1 needed 27 minutes for calling the flanks alone
  and 6.7 GB (CAMI III sample 0).
- `--save-map` / `--map`: keep the map of a sample and skip the first read pass next time.
- Contexts are called one gene per process; memory per gene instead of for all genes at once.

Packaging
- Python package `anchorwalk` with the commands `aw-flank` and `aw-flank-selftest`; the C program is built into the
  package. Conda recipe in `recipe/`, test suite in `tests/`. Built conda packages for macOS arm64, Linux x86-64 and
  Linux aarch64 (Python 3.10 to 3.13) in `dist/conda/`; same output on all three, byte for byte. Each package build
  runs the whole test suite against the package it has just made.
- A gzip file that was cut off is reported with every zlib version. With zlib 1.3.2 (the one conda-forge ships) the
  warning was missing: that version reports the unexpected end only with the read that meets it. Two files cut at
  the same pair would have been used without a word.
- Clear messages (and no left-over temporary folder) on a full disk, an unwritable output folder, Ctrl-C, a worker
  that fails or dies. Any number of genes (0.2 development versions failed above 126 on macOS). The same output in
  every run and for every Python version.
- Input checks: blanks and Windows line ends in the gene file and the reads; gene names must be plain ASCII; reads
  given as pipes are refused (they are read twice); a saved map is checked against the reads and for damage.

Changed defaults and names
- `--keep` no longer implies writing the reads near the genes; use `--reads-out` (`PREFIX.near_1.fq`, `PREFIX.near_2.fq`).
- The source layout changed: `src/aw_flank.py` and `src/aw_paths.py` became `anchorwalk/cli.py` and `anchorwalk/contexts.py`.

## 0.1 (2026-11)

Prototype, tested on simulated reads from one genome. Archived in `archive/`.
