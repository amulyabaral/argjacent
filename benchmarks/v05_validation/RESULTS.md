# ARGjacent v0.5.0 Expanded Validation Results (25 Unseen Samples Across CAMI II, CAMI III, and Fresh Real PacBio HiFi Metagenomes)

**Date:** 23 November 2026
**Protocol:** [PROTOCOL.md](PROTOCOL.md) (frozen `ARGjacent v0.5.0` and `argjacent-cohort v0.5.0` checksums recorded before running any sample; zero code changes made during validation).
**Storage management:** Zero multi-GB read archives or full HiFi files were stored on disk (streaming HTTP/S3 extraction and immediate FASTQ deletion kept ≥ 12 GiB free throughout; final archived outputs total 352 MB across all 25 samples).

---

## 1. Executive Summary Across All 4 New Suites (25 Unseen Samples)

| Suite (all untouched by `v0.5.0`) | Samples | 1-sided ≥ 1 kb (raw strict / unexplained faults) | 1-sided ≥ 5 kb | GFA v1.0 lossless (paths, ratio, mismatches) | Single-sample 2-sided `±500 bp` (ARGj vs Cartesian) | Single-sample 2-sided `±1 kb` (ARGj vs Cartesian) | Cohort Part A added pairs (`±500 / ±1 kb`) |
|---|---:|---|---|---|---|---|---|
| **Suite 1: CAMI II Strain Madness** (408 clinical MDRO strains) | 8 | 329 / 4,822 raw (`12` faults = `1 / 264` groups) | **0 / 292** (`0.00%`) | 10,455 paths, **1.65×**, **0** | **1 / 120** (`0.83%`) vs 16 / 249 (`6.43%`) | **1 / 77** (`1.30%`*) vs 18 / 173 (`10.40%`) | **0 / 26** (`0.00%`) / **0 / 22** (`0.00%`) |
| **Suite 2: CAMI III Human Gut reserve** (`14, 16, 18, 19`) | 4 | 47 / 1,083 raw (**0 faults**, `0 / 108` groups) | 9 / 188 raw (**0 faults**) | 2,998 paths, **1.50×**, **0** | **0 / 54** (`0.00%`) vs 10 / 89 (`11.24%`) | **0 / 42** (`0.00%`) vs 7 / 62 (`11.29%`) | **0 / 1** / 0 / 0 |
| **Suite 3: CAMI II Marine** (908 marine genomes + plasmids) | 4 | 9 / 19 raw† (**0 faults**, `0 / 12` groups) | 7 / 15 raw† (**0 faults**) | 193 paths, **1.20×**, **0** | **0 / 7** (`0.00%`) vs 0 / 7 | **0 / 6** (`0.00%`) vs 0 / 6 | 0 / 0 / 0 / 0 |
| **Suite 4: Fresh real paired Illumina + PacBio HiFi** (6 infant gut, 2 adult gut, Zymo D6331) | 9 | **9 / 1,369 (`0.66%`**, 95% CI `0.30–1.24%`; `3 / 132` groups; **1 / 1,369** structural overshoot) | **2 / 268** (`0.75%`) | 3,098 paths, **1.73×**, **0** | **0 / 58** (`0.00%`, CI `0–6.16%`) vs 9 / 101 (`8.91%`) | **0 / 45** (`0.00%`, CI `0–7.87%`) vs 9 / 70 (`12.86%`) | **0 / 2** / **0 / 2** (6-sample infant cohort) |
| **Pooled all 25 new samples** | **25** | **16 CAMI:** `1 / 384` unexplained group (`0.26%`)<br>**9 real HiFi:** `9 / 1,369` (`0.66%`) | **16 CAMI:** `0` faults<br>**9 real HiFi:** `2 / 268` | **16,744 paths, 1.63×, 0 mismatches** | **1 / 239 (`0.42%`, CI `0.01–2.31%`)** vs **35 / 446 (`7.85%`)** | **1 / 170 (`0.59%`, CI `0.01–3.23%`)** vs **34 / 311 (`10.93%`)** | **0 / 29 (`0.00%`)** / **0 / 24 (`0.00%`)** |

\* *In Strain Madness `s40`, the single `±1 kb` flag (`erm(C)_13_M13761`) was a real 2,400-bp circular plasmid (`0` mismatches over its 1,665-bp circular wraparound) cut short by linear SPAdes draft contigs; after hand-checking, true `±1 kb` chimeras are `0 / 76` on Strain Madness and `0 / 169` across all 25 new samples.*
† *In CAMI II Marine, all 9 linear-alignment flags (`dfrD`, `blaZ`, `str`) were small circular plasmids reconstructed at 100.00% identity past the linear FASTA cutpoint.*

### Combined Unseen Real PacBio HiFi Held-Out Totals (Block F + Suite 4 = 21 Unseen Real Metagenomes)
- **One-sided certified flanks ≥ 1 kb (primary, native DNA):** **9 / 2,757 = 0.33%** (exact 95% CI **0.15–0.62%**) across **3 / 259** independent 31-mer sequence groups (**1.16%**, 95% CI **0.24–3.35%**).
- **One-sided certified flanks ≥ 2 kb:** **2 / 1,364 = 0.15%**. At **≥ 5 kb:** **2 / 435 = 0.46%**.
- **Single-sample two-sided pairing (unseen real HiFi):** **0 / 69 = 0.00%** (95% CI **0.00–5.21%**) at `±500 bp`, and **0 / 53 = 0.00%** (95% CI **0.00–6.72%**) at `±1,000 bp` (closing the small `0 / 11` Block F sample-size gap).
- **Lossless GFA v1.0 across all 51 tested metagenomes (26 prior archive + 25 new):** **29,720 paths, 0 byte-for-byte round-trip mismatches (`100.000%` lossless)**.

---

## 2. Detailed Suite Breakdowns and Forensic Hand-Checks

### Suite 1 — CAMI II Strain Madness (8 samples: `0, 10, 20, 30, 40, 50, 60, 70`)
- **Dataset:** 408 closely related Charité Berlin clinical multidrug-resistant bacterial isolates (MRSA, ESBL *E. coli*, *Enterobacteriaceae*; 251 genomes carrying ResFinder ARGs across 995 families). Every sample contains 405–408 co-existing strains (20–29 at ≥ 5× coverage up to 71.8×, plus ~380 low-coverage background strains).
- **Runtime & memory:** 31–35 s per sample, ~850 MB peak RAM.
- **One-sided certified flanks:**
  - `≥ 5 kb`: **0 / 292** raw failures (`0 / 51` groups).
  - `≥ 1 kb`: 4,493 / 4,822 (`93.2%`) passed strict alignment directly (4,355 molecule-exact). Triage of the 329 raw flags across 24 / 264 groups: `189` `past_contig_end` (source SPAdes draft contigs ending inside shared regions), `91` `gene_end_diff` (flank strict from gene boundary, copy differs inside the gene), `21` `elsewhere_in_sample_99pct` (≥ 99% over ≥ 98% in unseeded locus), `16` `lenient_strain_diff` (closely related strain SNPs), and **12 `unexplained` claims in 1 single group (`1 / 264 = 0.38%`)**: `s70 tet(D)_1_AF467077` right (`1,919 bp`), which matches `G0281` (49.0×) at 100.00% identity over 92.1% (`1,767 bp`) up to the end of its 3,906-bp draft contig (`1,894 bp` available right flank).
  - **Recall (≥ 5× true contexts):** `1,908 / 2,428` (**78.6%**) found; `1,397 / 2,428` (**57.5%**) certified ≥ 1 kb; `393 / 2,428` (**16.2%**) certified ≥ 5 kb.
- **Two-sided pairing & forensic diagnosis of the single `±500 bp` chimera:**
  - Single-sample `±500 bp`: **1 / 120** (`0.83%`) vs Cartesian **16 / 249** (`6.43%`), Main-to-Main **6 / 203** (`2.96%`), Greedy **3 / 210** (`1.43%`).
  - **Forensic mechanism of the 1 false pair (`s30 aadA2_2_JQ364967`, `tier = span_pair`, `pairs = 10`):** `aadA2` is 792 bp long (too long for a 2×150 short-read pair to physically span). `L1` (`10×`, carrying `catB3–blaOXA-1–aac(6')-Ib-cr`) and `R1` (`67×`, carrying the Class 1 integron `qacEdelta1–sul1` 3'-CS) occur directly adjacent to each other *without* `aadA2` on a separate co-existing integron cassette array in the same sample; 10 read pairs from that `aadA2`-lacking integron bridged `L1` to `R1`, fooling `span_pair` because `pair.py` does not disable `span_pair` when `len(gene) > 650 bp` during `snv_phase` scanning.
  - **8-sample `argjacent-cohort` Part A (JSON-only `depth_bins` mode):** **0 / 143** (`0.00%`, CI `0.00–2.55%`) at `±500 bp`, including **0 / 26** newly resolved cohort pairs (`cohort_flow` 0/20, `cohort_anchor` 0/6), nearly **doubling multi-context recall** from `7.6%` (`14/184`) to `14.1%` (`26/184`). At `±1,000 bp`: **0 / 22** newly resolved cohort pairs (`cohort_flow` 0/12, `cohort_anchor` 0/10), doubling multi-context recall from `5.5%` (`10/182`) to `11.0%` (`20/182`).

### Suite 2 — CAMI III Toy Human Gut Untouched Reserve (`14, 16, 18, 19`)
- **One-sided certified flanks:** **0 unexplained faults** across `1,083` certified claims ≥ 1 kb (`108` groups; `1,030 / 1,083 = 95.1%` molecule-exact) and `2,627` claims ≥ 100 bp (`205` groups). Raw flags ≥ 1 kb (`47`): `28` `gene_end_diff`, `17` `elsewhere_in_sample_99pct`, `2` `past_contig_end`.
- **Recall (≥ 5×, circular plasmid rotations collapsed):** `386 / 478` (**80.8%**) found; `196 / 478` (**41.0%**) certified ≥ 1 kb; `69 / 478` (**14.4%**) certified ≥ 5 kb.
- **Two-sided pairing:** ARGjacent **0 / 54** (`0.00%`) at `±500 bp` and **0 / 42** (`0.00%`) at `±1,000 bp`, vs Cartesian **10 / 89** (`11.24%`) / **7 / 62** (`11.29%`), Greedy **4 / 75** (`5.33%`) / **4 / 54** (`7.41%`), Main-to-Main **2 / 73** (`2.74%`) / **3 / 54** (`5.56%`). Cohort Part A: **0 / 55** (`±500 bp`) and **0 / 42** (`±1,000 bp`).

### Suite 3 — CAMI II Marine (`0, 1, 2, 3`; 908 active marine genomes + plasmids)
- **Runtime & memory:** 26–30 s per 16.6 M read-pair sample, 2.15–2.22 GB RAM.
- **One-sided & two-sided:** **0 faults** across all 19 certified claims ≥ 1 kb and all 50 claims ≥ 100 bp (all 9 linear flags were small circular plasmids carrying `dfrD`, `blaZ`, or `str` reconstructed at 100.00% identity around the circular plasmid wraparound). Two-sided pairing: **0 / 7** (`±500 bp`) and **0 / 6** (`±1,000 bp`). GFA: `193` paths, `1.20×`, **0** mismatches.

### Suite 4 — 9 Fresh Unseen Real Paired Illumina + PacBio HiFi Metagenomes
Per-sample summary (exact frozen `claims_lr3.py` v3.1 + `truth_2side.py` single-molecule HiFi scoring):

| Sample | Project / Biome | KMA genes | 1-sided ≥ 1 kb (primary / groups) | 1-sided ≥ 5 kb | GFA ratio (mismatches) | 2-sided `±500 bp` (ARGj vs Cart) | 2-sided `±1 kb` (ARGj vs Cart) |
|---|---|---:|---|---|---|---|---|
| `INF134_8G1` | `PRJNA1365828` infant gut | 29 | 2 / 153 (`1 / 21`) | 2 / 60 | 1.70× (0) | **0 / 9** vs 2 / 14 | **0 / 6** vs 1 / 9 |
| `INF133_3G1` | `PRJNA1365828` infant gut | 23 | 6 / 247 (`1 / 25`) | 0 / 52 | 1.97× (0) | **0 / 14** vs 0 / 17 | **0 / 11** vs 1 / 14 |
| `INF133_8G1` | `PRJNA1365828` infant gut | 21 | **0 / 195** (`0 / 18`) | **0 / 31** | 1.83× (0) | **0 / 6** vs 2 / 10 | **0 / 4** vs 1 / 6 |
| `INF134_2G1` | `PRJNA1365828` infant gut | 21 | **0 / 235** (`0 / 19`) | **0 / 42** | 1.81× (0) | **0 / 9** vs 2 / 17 | **0 / 8** vs 3 / 15 |
| `INF134_3G1` | `PRJNA1365828` infant gut | 16 | **0 / 196** (`0 / 17`) | **0 / 15** | 1.79× (0) | **0 / 9** vs 2 / 13 | **0 / 8** vs 2 / 12 |
| `INF134_4G1` | `PRJNA1365828` infant gut | 21 | **0 / 107** (`0 / 11`) | **0 / 10** | 1.81× (0) | **0 / 6** vs 2 / 10 | **0 / 4** vs 0 / 6 |
| `US_Adult2T` | `PRJNA754443` US adult gut | 19 | **0 / 85** (`0 / 5`) | **0 / 26** | 1.53× (0) | 0 / 0 vs 0 / 9 | 0 / 0 vs 0 / 1 |
| `UK_Adult` | `PRJEB90666` UK adult gut | 19 | 1 / 67 (`1 / 8`) | **0 / 6** | 1.50× (0) | **0 / 1** vs 1 / 7 | **0 / 1** vs 1 / 2 |
| `Zymo_D6331` | Zymo Gut Standard (5 *E. coli*) | 4 | **0 / 84** (`0 / 8`) | **0 / 26** | 2.16× (0) | **0 / 4** vs 0 / 4 | **0 / 3** vs 1 / 5 |
| **Total (9 real)** | **4 projects** | **173** | **9 / 1,369 (`0.66%`, CI `0.30–1.24%`; `3 / 132` grp)** | **2 / 268 (`0.75%`)** | **1.73× (0 / 3,098)** | **0 / 58 (`0.00%`)** vs **9 / 101 (`8.91%`)** | **0 / 45 (`0.00%`)** vs **9 / 70 (`12.86%`)** |

- **Base-level forensic inspection of all 9 failed `≥ 1 kb` real HiFi claims (3 distinct loci):**
  1. `INF134_8G1` `erm(B)` left `c3` (2 allele claims, `5,403 bp`, verdict `detail`): 65 single Revio HiFi reads span the full 5,403 bp with only **1 single SNP** plus one **72-bp short tandem repeat indel** at `5,105 bp` (`99.98%` identical outside that 72-bp unit).
  2. `INF133_3G1` `mecA` right `c6–c11` (6 duplicate candidates sharing one sequence, `1,862–1,874 bp`, verdict `detail`): 284 single Revio HiFi reads span the full 1,874 bp with **0 SNPs** plus one **80-bp direct-repeat unit indel** at `1,708 bp` (the well-known 40-bp tandem repeat *S. aureus* SCC*mec* `dru` locus).
  3. `UK_Adult` `lnu(C)` left `c2` (1 claim, `1,228 bp`, verdict `contradicted`): 61 HiFi reads match with 0 mismatches to `1,154 bp` and then depart — **1 genuine structural overshoot of 74 bp (`1 / 1,369 = 0.07%`, `1 / 132` groups = `0.76%`)**.
- **6-sample real longitudinal infant gut cohort (`PRJNA1365828`, JSON-only `argjacent-cohort`):**
  - Part A pairing: **0 / 55** chimeras at `±500 bp` (including `0 / 2` newly resolved cohort pairs: `cohort_anchor` 0/1, `cohort_flow` 0/1) and **0 / 43** chimeras at `±1,000 bp` (`0 / 2` newly resolved).

---

## 3. Honest Evaluation of `argjacent-cohort` Part B (Post-Repeat Exit Deconvolution)

Testing `argjacent-cohort` Part B beyond isolated single-repeat IS26 benchmarks revealed an important, mechanistically clear limitation when real genomes contain **chained repeats**:
- **CAMI II Strain Madness (8 samples):** `113 / 157` (`72.0%`) passed strict; `9` had pre-existing single-sample flags (`gain = 0` or draft contig end); **`35 / 148` (`23.6%`)** newly extended exits failed strict verification.
- **CAMI III Human Gut (4 samples):** `6 / 7` passed strict (extending `tet(Q)` in `s19` cleanly from `545 bp` to `4,970 bp`), while **`1 / 5`** newly extended exits (`mef(A)` in `s19`) failed strict.
- **Real infant gut HiFi cohort (6 samples):** `6 / 8` passed strict (all 6 had `gain = 0 bp`, matching single-sample `certain_bp`), while both newly extended exits (**`2 / 2`**, `fosB4` right in `INF134_3G1` `+204 bp` and `INF134_4G1` `+2,154 bp`) failed strict HiFi verification.
- **Root cause:** In `argjacent/cohort.py` (`_measure_repeat_from_bins`), after NNLS validates co-abundance across the first repeat (`fork_bp`), `cohort_certified_bp` extends up to 5 kb using only coarse 100-bp depth stability (`0.30×–1.85×`) without re-checking downstream forks or joins. Consequently, while Part B works accurately across isolated single repeats (0 / 108 errors on synthetic IS26), in complex clinical metagenomes with chained mobile elements it should be treated as **exploratory / hypothesis-generating** until downstream join verification is added in a future release, whereas one-sided certification (`0.33%` across 21 unseen real metagenomes) and Part A two-sided pairing (`0 / 69` on unseen real HiFi; `0 / 29` cohort-added pairs across 18 cohort samples) meet the `< 1%` reliability bar.
