# Pre-Registered Validation Protocol for Frozen ARGjacent v0.5.0 (CAMI II, CAMI III Untouched, and Fresh Real HiFi)

**Written:** 23 November 2026, before downloading or running any sample named below.
**Purpose:** Evaluate the complete frozen `argjacent v0.5.0` release — one-sided certified flanks, lossless GFA v1.0 compaction, single-sample two-sided `[left + gene + right]` pairing with `depth_bins`, and multi-sample `argjacent-cohort` (Part A pairing + Part B post-repeat exit resolution in JSON-only mode) — across proper official CAMI II challenge datasets, the untouched CAMI III Toy Human Gut held-out reserve, and fresh unseen real paired Illumina + PacBio HiFi metagenomes, under strict zero-tarball streaming storage management.

## 1. Frozen Software Checksums (SHA-256)

`argjacent v0.5.0` and `argjacent-cohort v0.5.0` are frozen for the entire validation. No rule, threshold, or code path in `argjacent` is modified during this benchmark; every failure or chimera is reported verbatim alongside its triage/hand-check.

- `argjacent/bin/aw` : `8a2e418c2740836b49de396b9ad06f06469894a81327520fe04a0ae3ccd58d77`
- `argjacent/bin/awx` : `ce7dad8f58192483f5ed952b95a042d487ae0b73ee90cbf0a16bb2aa13205707`
- `argjacent/bin/awv` : `90ea98473c41527ece8c140b4d07ed023d0796cbd9f394cc7e780021f39b828a`
- `argjacent/cli.py` : `2131b5d191ad0494d79fa5df9c64a0ad3ffda94273ba7c795681d75ccacec22b`
- `argjacent/verify.py` : `0c5f28c022251f88d4899904a9727a14e906f7624ea15f6a7acec650a113c72c`
- `argjacent/pipeline.py` : `c0cfe0df3562ba411296c786a45bfc76aaa3f3b8edfb96291dfff7f9095973b7`
- `argjacent/gfa.py` : `eaa23575a51cbf7f1818a3a6e8803199baade49e607230455e646443bcb23c8d`
- `argjacent/pair.py` : `313d2960eb1cd545624afdd8bf607e93a50c10fd5fda1e025965a905084c42c9`
- `argjacent/cohort.py` : `20259a2b7526fc8b176d9b097df499fb7d9f3d78b5334feeeac56c3f26832745`
- `benchmarks/v05_validation/score_v05_cami.py` : `691699ffb1bef4c01adfe2dd34932702ebc239df293824b996d5b8612b7a91ef`
- `benchmarks/v05_validation/score_v05_cohort.py` : `14b7a02049e490163d4423e5f05c18ed7c5234b9525b70b8f8f98eab5fe6e664`
- `benchmarks/v05_validation/run_cami_suite.py` : `c8c7e20207fceef4fd345b591ba1ca61e539c804d67be9a2efad0eca8d119e91`
- `benchmarks/v05_validation/run_real_suite.py` : `61dedeeaa11bfb2601672f18a8c87f2e52854f656164566ba6b08f1340a25fdf`

Command run on every sample:
- Seed identification: `kma -ipe R1.fq.gz R2.fq.gz -o kma -t_db kmadb -1t1 -cge -ef -nf -t 8` against ResFinder (`all.fsa`, 3,212 sequences), keeping consensus templates with `>= 90%` identity and `>= 60%` coverage (`benchmarks/cami/kma_genes.py`).
- ARGjacent: `bin/argjacent -g genes.fa -1 R1.fq.gz -2 R2.fq.gz -o run -D 10000 -t 8 --max-mem 7` (with `--keep` retained only for the ~20 s needed to score `W = 1,000 bp` pairing, then deleted immediately).
- Multi-sample cohort: `bin/argjacent-cohort --sample-sheet sheet.tsv -o cohort --pair-bp {500,1000}` in pure JSON-only mode (`depth_bins` + `gene_depth` from `run.flanks.json`; zero read files).

## 2. Pre-Registered Datasets and Sample Selection Rules

All sample IDs below were fixed **before** downloading any sample reads or running ARGjacent:

1. **Suite 1 — CAMI II Strain Madness (`strmgCAMI2`, Charité clinical MDRO multi-strain benchmark; 408 closely related clinical strains, 251 carrying ResFinder ARGs):**
   - **Rule:** 8 evenly spaced samples across the 100-sample series: **samples `0, 10, 20, 30, 40, 50, 60, 70`** (`2.07 GB` reads each).
   - **Cohort test:** All 8 samples jointly via `argjacent-cohort` (Part A pairing + Part B post-repeat exits).
   - **Why:** Official CAMI II benchmark specifically designed to break strain-level assembly and flanking-context callers, with 405–408 closely related MRSA / ESBL / *Enterobacteriaceae* clinical strains co-existing in every sample (20–29 abundant strains `>= 5x` plus ~380 lower-abundance strains).

2. **Suite 2 — CAMI III Toy Human Gut, Reserved Untouched Held-Out Samples (`cami3_toydata/human-gut-toy`):**
   - **Rule:** All 4 samples explicitly held in reserve in `benchmarks/heldout/PROTOCOL.md` section 8: **samples `14, 16, 18, 19`** (`4.87 GB` reads each). Never downloaded or inspected by any version of ARGjacent (`v0.1`–`v0.5`).
   - **Cohort test:** All 4 samples jointly via `argjacent-cohort` (Part A pairing + Part B post-repeat exits).

3. **Suite 3 — CAMI II Marine (`marmgCAMI2`, environmental metagenome; 977 marine genomes + plasmids, 72 carrying ResFinder ARGs across 160 families):**
   - **Rule:** First 4 samples of the 10-sample series: **samples `0, 1, 2, 3`** (`5.17 GB` reads each).
   - **Cohort test:** All 4 samples jointly via `argjacent-cohort` (Part A pairing + Part B post-repeat exits).

4. **Suite 4 — Fresh Unseen Real Paired Illumina + PacBio HiFi Metagenomes (`run_real_suite.py`, 9 brand-new unseen samples across 4 projects, never used in any prior block A–F):**
   - **Rule:** Pre-registered before downloading any read file; every sample has direct ENA HTTPS paired-end Illumina FASTQs (`<= 4.7 GB` total) and PacBio Sequel II / Revio HiFi reads (`>= 6.7 Gbp`):
     - **Longitudinal infant gut cohort (`PRJNA1365828`, 6 untouched timepoints reserved in `PROTOCOL.md` section 8):**
       - `INF134_8G1` (`SRR36082382` 2.02 GB / `SRR36082377` 11.5 Gbp Revio)
       - `INF133_3G1` (`SRR36082375` 2.17 GB / `SRR36082374` 7.8 Gbp Revio)
       - `INF133_8G1` (`SRR36082386` 2.03 GB / `SRR36082371` 7.7 Gbp Revio)
       - `INF134_2G1` (`SRR36082367` 2.07 GB / `SRR36082381` 7.2 Gbp Revio)
       - `INF134_3G1` (`SRR36082366` 1.80 GB / `SRR36082384` 7.1 Gbp Revio)
       - `INF134_4G1` (`SRR36082365` 1.78 GB / `SRR36082380` 6.8 Gbp Revio)
       - *Cohort test:* all 6 samples jointly via `argjacent-cohort` (JSON-only mode), validated directly against real single-molecule PacBio Revio HiFi reads.
     - **US adult gut (`PRJNA754443`, untouched `2_treatment`):** `SRR15489030` (3.90 GB) / `SRR15489017` (7.4 Gbp Sequel II).
     - **UK adult gut (`PRJEB90666`, untouched `21BR240s1`):** `ERR15463326` (4.67 GB) / `ERR15479330` (6.7 Gbp Revio).
     - **ZymoBIOMICS Gut Microbiome Standard D6331 (`PRJNA1446191` / `PRJNA680590`, 21 staggered gut strains including 5 co-existing *E. coli* strains):** `SRR37882064` (3.20 GB) / `SRR13128014` (18.0 Gbp Sequel II).

## 3. Endpoints and Scoring Rules (Frozen Before Running)

Evaluated by `benchmarks/v05_validation/score_v05_cami.py` and `score_v05_cohort.py`:

1. **One-sided certified flank accuracy (`context_len` = `certain_bp`, and separately `end_len` = `certain_end_bp`):**
   - Strict alignment verification matching `claims_lr3.py` v3.1 (`<= 3%` differing bases per 3-kb window, `<= 25` differing bases per 300-bp window, no indel `>= 50 bp`, from `-100 bp` inside the seed gene to `certain_bp`) against all source genomes present in that sample (`cov > 0`).
   - Reported at `>= 100 bp`, `>= 1 kb`, `>= 2 kb`, and `>= 5 kb` — both **per claim** and **per independent 31-mer sequence group**, with exact **95% Clopper–Pearson confidence intervals**.
   - Every single failing claim is automatically triaged (`gene_end_diff`, `past_contig_end`, `elsewhere_in_sample_99pct`, or `unexplained`) and inspected by hand; primary raw numbers are always reported alongside triaged numbers.
   - Structural recall on all true contexts with `cov >= 5x` at `>= 1 kb`, `>= 2 kb`, `>= 5 kb`, and full `10 kb` (`min(L, achievable)`), reported both raw and with CAMI III circular plasmid rotations collapsed.

2. **Lossless GFA v1.0 compaction:**
   - Independent parse of `run.gfa`, checking `100%` byte-for-byte path reconstruction of every `P` line against `genes.fa` + `run.flanks.json`, plus compaction ratio (`raw_bp / compacted_bp`).

3. **Two-sided `[left + gene + right]` pairing (`W = 500 bp` and `W = 1,000 bp`):**
   - Direct evaluation of `run.pairs.tsv` (`v0.5.0` with `depth_bins`) and `W = 1,000 bp` against exact source genomes (`cov > 0`), compared head-to-head against Cartesian, Main-to-Main, and Unguarded Greedy Depth baselines on the identical `run.flanks.json`.
   - Reports chimera rate (`false pairs / checkable`, 95% CI, plus distinct family-pairs), overall recall, multi-context (`>= 2` true contexts) recall, tier breakdown, and number of pairs unlocked by `v0.5.0` `post_fork_bins`.

4. **Multi-sample `argjacent-cohort` (JSON-only `depth_bins` mode):**
   - **Part A (`cohort_pairs.tsv`):** chimera rate and multi-context recall of cohort-added pairs (`cohort_flow`, `cohort_anchor`, `cohort_complement`) vs single-sample pairing.
   - **Part B (`cohort_repeats.tsv`):** strict alignment verification (`claims_lr3.py` v3.1 rules) of every post-repeat certified exit arm (`cohort_certified_bp`) and mean certified extension gain beyond single-sample `certain_bp`.

## 4. Storage Management Policy

- Zero multi-GB `.tar.gz` archives are ever stored on disk (HTTP stream -> on-the-fly tar extraction -> deinterleave -> `R1.fq.gz` / `R2.fq.gz`).
- Disk guard enforces `>= 6.5 GiB` free before every sample starts.
- `R1.fq.gz` and `R2.fq.gz` are deleted immediately after `argjacent` completes; `run.tmp/` is deleted immediately after `score_v05_cami.py` completes; per-suite temporary genomes + BLAST DBs are deleted immediately after each suite's `argjacent-cohort` scoring finishes.
- Only compact JSON/TSV/GFA/log outputs (~5–15 MB/sample) are retained in `benchmarks/v05_validation/results/`.
