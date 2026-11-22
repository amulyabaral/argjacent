> *Historical note: Versions 0.1–0.4.1 were developed and held-out under the name `aw-flank` (`bigargs-anchorwalk`) before renaming to **ARGjacent** in v0.5.0.*

# Held-out validation of aw-flank 0.2.1

How often are aw-flank's two claims wrong on data that no version of the tool has seen, with the code frozen? This is
the report of that test, carried out on 14-17 November 2026 by the protocol in `benchmarks/heldout/PROTOCOL.md`
(written before any of its data was fetched; its three addenda were written as the work went on, each before the
results it could affect were scored). Per-sample tables, scripts and the hand-check log are in `benchmarks/heldout/`.

## 1. The claims

For every context of a gene the tool reports two lengths, counted outward from the gene's end:

- **certain** (`certain_bp`): the order of sequence is that of a real DNA molecule next to this gene; bases may
  differ where strains differ;
- **exact** (`exact_bp`): base for base the sequence of one real molecule.

Everything beyond the certain part is a candidate, not a claim.

## 2. What was frozen

The built conda package (`bigargs-anchorwalk-0.2.1-py312h8a0c2ca_0.conda`, sha256 `5a39d33c...af212c6`), default
options, every sample run once with it (and once more with `--strict`). The tool was not changed during the test.

## 3. Data

Fourteen real samples with short reads and PacBio HiFi reads of the same sample, and four simulated ones. None had
been seen by any version of the tool. Every sample of the protocol's tables was run; the one that gave no result is
listed.

| Sample | What | Read pairs | Resistance genes found | 
|---|---|---:|---:|
| H1, H2 | adult gut, Japan (the cohort of the development sample) | 32.5 M, 55.0 M | 47, 36 |
| H3, H4 | infant gut | 45.4 M, 16.7 M | 17, 15 |
| H5 | adult gut, USA | 12.4 M | 7 |
| H6 | adult gut, Korea | 54.3 M | 39 |
| H7 | adult gut, UK (reported apart, as the protocol says) | 50.4 M | 16 |
| B-A1 | municipal waste water, influent | 32.3 M | 102 |
| B-A2 | pig faeces | 30.9 M | 78 |
| B-A4 | chicken caecum | 44.8 M | 30 |
| B-A5, B-A6 | two anaerobic digesters (B-A6: 2×300 reads) | 38.1 M, 23.5 M | 21, 20 |
| B-A10 | cheese rind (2×251 reads) | 29.5 M | 4 |
| mock | 71 strains mixed from purified DNA (reported apart) | 20.6 M | 1 |
| B-A3 | activated sludge | 50.3 M | **no result**: the tool's read pass needs more than this laptop's 16 GB |
| CAMI III 9, 10, 11, 13 | simulated human gut (toy set), 0.3 % read errors | 16.6 M each | |

The "main set" of the protocol is H1 to H6 and the six samples B-A1 to B-A10 that ran: twelve samples. All reads were
trimmed with fastp; genes were found with KMA against ResFinder, as a user would.

## 4. How a claim was checked

**Truth.** For real samples, PacBio HiFi reads of the same sample: a long read that carries the end of the gene is a
real molecule next to it. For simulated samples, the source genomes (and, for a strict reading, their neighbourhood of
every gene copy taken as "long reads"). For the mock community, both.

**The scorer** (`benchmarks/heldout/scorer/claims_lr3.py`, version 3.1). The long reads are aligned to each claim
(the gene's last 300 bp + the certain part). A claim is

- *confirmed* if one read holds all of it: at most 3 % different over any 3 kb, at most 25 differing bases in any
  300 bp (insertions and deletions counted by their length), no insertion or deletion of 50 bp or more;
- *tiled* if every 3-kb window of it is held so by some read;
- *differs* if the reads follow it only leniently (somewhere a stretch differs by 50 bp to 3 kb, or by more than
  8 % locally);
- *contradicted* if support stops and three or more of the reads that held it until then go on differently;
- *not checkable* otherwise (the reads end, or none follows it).

*Failed* = differs + contradicted, out of the checkable claims. An exact part is *refuted* if a 31-mer of it occurs
in no long read that carries this gene end, or if no long read matches it with differences of its own only (a
difference counts as a read's own error only if no other read shares it and at least two other reads show the
claimed base there).

**How the scorer itself was tested.** Before any held-out sample was scored:

- About 1,400 errors were planted into the development sample's claims: wrong joins, inversions, stretches of 300
  to 2,500 bp replaced or left out, single wrong or missing bases.
- Two reviews were made: an adversarial audit of the code, with its own labelled test sets (218 mutants of real
  claims, synthetic samples), and an independent second checker written without the scorer's code.
- The project's earlier scorer turned out to pass 90 to 99 % of claims with a stretch replaced or left out. The
  version used here passes none of the planted structural errors that are errors (one planted join that passes is
  carried whole by 325 long reads, so it exists in the sample), refutes every planted wrong base in the part it can
  check, and still catches the four claims of 0.2.0 that the long reads contradict.
- The same planted errors were put into every held-out sample's claims as well (section 6).
- The independent checker was run on every held-out sample; every disagreement was looked at.

**Hand checks.** Every failed claim and every refuted exact part was looked at, first by an automated sorting
(`real_triage.py`, `sim_triage.py`), then by eye (`benchmarks/heldout/results/HANDCHECKS.md`). The sorting asks:

- Does the claim occur exactly in a source genome?
- Does the sample's copy of the gene differ from the seed inside the gene, while the flank holds from the gene's
  boundary?
- Is the order of sequence that of a real molecule with its bases differing (strains, a paralogous copy)?
- Is it a rare context that a 3-kb window cannot tile?
- Or is it a fault?

The primary numbers keep the scorer's verdict; the hand-checked numbers are given next to them.

**Statistics.** Claims of one place stand or fall together, so besides the exact 95 % interval on claims, an
interval with **groups** as the unit is given: claims whose certain parts share a 31-mer are one group. With no
failure in m groups the exact one-sided upper limit is about 3/m.

**Known blind spots of the check** (measured on labelled data):

- a wrong exit from a repeat longer than 3 kb;
- a wrong join that one chimeric long read happens to carry;
- a wrong base that one read's error shares;
- a true rare context next to an abundant one, which is called contradicted when no read holds the gene and 3 kb of
  it (the hand check is the safeguard).

## 5. Results

### Real samples: the primary endpoints

Scorer verdicts as they stand, pooled over the main set (`benchmarks/heldout/gather_final.py`, `scorer/pool3.py`):

| Endpoint | Failed / checked | Rate | Exact 95 % interval | By groups (places) |
|---|---:|---:|---:|---:|
| **P1: certain claims of 1 kb or more** | **34 / 1,127** | **3.0 %** | 2.1-4.2 % | 5 of 132 groups (3.8 %, 1.2-8.6 %) |
| certain claims of 100 bp or more | 82 / 3,529 | 2.3 % | 1.9-2.9 % | 12 of 351 |
| certain claims of 5 kb or more | 0 / 79 | 0 % | 0-4.6 % | 0 of 26 |
| **P3: exact parts** | **11 / 1,182** | **0.9 %** | 0.5-1.7 % | |

Six claims of 1 kb or more of 1,133 were not checkable.

**By the rule set down beforehand, the README's statement does not hold.** It said to treat "certain" as right in
roughly 99 of 100 contexts. The protocol asked for at most 1 % failed claims of 1 kb or more and an upper limit
below 3 % with places as the unit; measured are 3.0 % and an upper limit of 8.6 %. The statement is replaced by the
measured rate: **about 97 of 100 certain claims of 1 kb or more are borne out by long reads (95 % interval 96 to 98),
and the failures sit at about 1 place in 26.**

Per sample:

| Sample | Claims ≥ 1 kb: failed / checkable | Places failed / all | Claims ≥ 100 bp | Exact parts refuted / checked |
|---|---:|---:|---:|---:|
| H1 | 5 / 276 | 1 / 20 | 23 / 761 | 2 / 338 |
| H2 | 6 / 155 | 1 / 11 | 6 / 576 | 2 / 145 |
| H3 | 0 / 35 | 0 / 14 | 0 / 45 | 0 / 24 |
| H4 | 4 / 114 | 1 / 11 | 4 / 128 | 2 / 84 |
| H5 | 0 / 39 | 0 / 4 | 0 / 77 | 0 / 36 |
| H6 | 0 / 95 | 0 / 18 | 5 / 318 | 0 / 106 |
| B-A1 | 3 / 162 | 1 / 22 | 8 / 559 | 3 / 143 |
| B-A2 | 0 / 49 | 0 / 10 | 1 / 434 | 1 / 132 |
| B-A4 | 16 / 104 | 1 / 10 | 35 / 306 | 0 / 89 |
| B-A5 | 0 / 31 | 0 / 4 | 0 / 105 | 0 / 18 |
| B-A6 | 0 / 51 | 0 / 6 | 0 / 199 | 1 / 46 |
| B-A10 | 0 / 16 | 0 / 2 | 0 / 21 | 0 / 21 |
| *H7 (apart)* | 0 / 38 | 0 / 5 | 8 / 141 | 0 / 39 |
| *mock (apart)* | 0 / 7 | 0 / 2 | 0 / 7 | 0 / 1 |

By kind of sample, claims of 1 kb or more: adult gut 11 / 565 (1.9 %), infant gut 4 / 149 (2.7 %), waste water
3 / 162 (1.9 %), animal gut 16 / 153 (10.5 %, one place on one sample), digester 0 / 82, food 0 / 16. With
`--strict`: 34 / 964 (3.5 %; the same five places, fewer claims).

### What the failed claims are (hand check)

Every failed claim was read against the long reads (`benchmarks/heldout/results/HANDCHECKS.md`). The 34 failed
claims of 1 kb or more are five places:

| Place | Claims | What the long reads show |
|---|---:|---|
| B-A4, ant(6)-Ia right | 16 | **Fault**: the certain part runs 10 bp past the place where the gene's reads leave the path (a step back from a fork that was 10 bp too short) |
| H2, catP right | 6 | **Fault**: for its last 1.5 kb the certain part follows a related copy at another locus (12-35 differing bases per 300 bp) |
| B-A1, aph(3'')-Ib / aph(6)-Id left | 3 | **Fault**: right for 3,725 bp, then about 200 bp to 1 kb that no long read with the gene carries |
| H1, sul3 right | 5 | Order of sequence right; in the last 450 bp every long read differs by 16-21 bases per 300 bp (the bases of a related copy) |
| H4, fosA left | 4 | **Limit of the check**: the sample's copy of the gene differs from the database sequence inside the gene; the flank holds from the gene's boundary |

So by hand **25 of 1,127 claims of 1 kb or more are faults (2.2 %), at 3 places**; 5 more have the right order and
wrong bases at their end; 4 are not faults. None of the faults is a wrong context from the start: each claim is
right for most of its length and wrong at its far end (10 bp, about 1 kb, 1.5 kb). All three places are where reads
of another locus had come onto the path. Under 1 kb there are 48 more failed claims at 7 places; they are listed in
the hand-check log (faults at 5 places; a limit of the check at 1; right order with differing bases at 1).

Of the 11 refuted exact parts, 9 are faults at 5 places: a consensus of related contexts (3), one or two bases at
the gene's boundary where the sample's gene begins or ends differently from the database sequence (3), a context
that the long reads do not carry (2), the bases of another version of the gene (1). Two are limits of the check.

### How much is claimed

Over the main set (511 gene-family sides):

- A context that is certain and borne out for 1 kb or more: 139 sides (27 %); for 5 kb or more: 27; for 8 kb: 12.
- Of the long reads that carry a gene with 1 kb or more of flank (122,496 molecules), a candidate follows 49 % for
  1 kb, and a certain part covers 18 %. At 5 kb: 17 % and 1.7 %.
- The `main` candidate taken whole is right for 1 kb on 321 of 363 sides where it can be checked, for 5 kb on 106
  of 173, for 8 kb on 73 of 128: a good guess, not a claim.
- Time and memory on a 16-GB laptop: 49 seconds to 13 minutes and 2.0 to 6.6 GB per sample.

### Simulated samples

CAMI III samples 9, 10, 11 and 13, scored against the source genomes. By the lenient reading that the README's
simulated numbers use (a claim fails if no source genome supports it), 5 of 2,021 contexts are flagged (0.25 %):
four for a wrong base in the exact part, at two places, and one claim under 1 kb that no genome holds. By the strict reading used for the real samples (the genomes' neighbourhoods
taken as long reads), 32 of 665 claims of 1 kb or more fail (4.8 %): 0 of 214, 0 of 183 and 0 of 133 on samples 9,
10 and 13 after the hand check (their 5 failures are claims that occur exactly in a source genome given as ten
rotations of a plasmid), and 27 of 135 on sample 11, of which 9 are faults of order at two ends of draft contigs,
14 are stretches of another genome's copy of the sequence, and 4 are bases of mixed strains just over the limit.
Sample 11 holds many strains that share erm(B), erm(G) and cfxA.

## 6. Do the checks catch errors on these samples?

Every claim of 1 kb or more of every real sample was damaged in ten ways and scored again, and one base was
changed and one deleted in every exact part (fourteen samples):

| Planted error | Planted | Passed | Called differs | Called contradicted | Not checkable |
|---|---:|---:|---:|---:|---:|
| outer part replaced by another gene's context | 1,155 | 0 | 0 | 1,111 | 44 |
| outer part replaced by random sequence | 1,178 | 0 | 0 | 1,134 | 44 |
| outer part reverse-complemented | 1,178 | 0 | 1 | 1,132 | 45 |
| 300 / 1,000 / 2,500 bp replaced | 1,155 / 773 / 193 | 0 / 0 / 0 | 4 / 0 / 0 | 1,101 / 727 / 175 | 50 / 46 / 18 |
| 300 / 1,000 / 2,500 bp left out | 1,178 / 784 / 200 | 0 / 0 / 0 | 320 / 12 / 2 | 834 / 719 / 180 | 24 / 53 / 18 |
| join to another context through shared sequence | 50 | 6 | 2 | 38 | 4 |
| **all structural errors** | **7,844** | **6** | 341 | 7,151 | 346 |
| one wrong base in the exact part | 795 (760 could be checked) | 0 | | 760 refuted | |
| one deleted base in the exact part | 795 (758 could be checked) | 0 | | 758 refuted | |

The six that pass are all joins through shared sequence, on three samples (B-A2, B-A4, B-A5). Five of them are
spanned whole by 6 to 145 long reads, gene end included, so the "planted error" is a real context of that sample.
The sixth (B-A5, lnu(B) left, candidate 2) is held only piece by piece and was not resolved; it may be the blind
spot named in section 4 (a wrong exit from a repeat). Every other planted error is caught in the part the check can
look at; 4 % of the damaged claims could not be checked.

## 7. Other tools on the same data

Full tables: `benchmarks/heldout/results/COMPARISON.md`. On three real samples (H1, B-A2, H4), claims or
flanks of 1 kb or more that the long reads do not bear out:

| Tool | Failed / checkable | |
|---|---:|---|
| aw-flank, certain parts | 9 / 439 (2.1 %) | none a fault of order |
| MEGAHIT contigs (flanks cut from them) | 31 / 138 (22 %) | |
| ARGextender | 11 / 79 (14 %) | |
| Sarand (k141 graph) | 106 / 204 (52 %) | |
| ARGContextProfiler (k141 graph) | 105 / 208 (50 %) | |

aw-flank's *candidates* hold more true molecules than any other tool's output on H1 and B-A2 and are level with the
best on H4; its *certain parts* cover far fewer, and almost nothing at 5 kb (on H4: 0.2 % of the molecules, where
MEGAHIT's and ARGextender's flanks cover 60 % and 76 %). The four simulated samples are in the same file.
The graph tools need a MEGAHIT assembly first (14 minutes to 1 hour 45 minutes and up to 8 GB here); MetaCherchant
ran out of its heap on one sample.

## 8. What this does and does not show

**What it shows.**

- On twelve fresh real samples of six kinds, 97 % of the claims "certain for 1 kb or more" are borne out by long
  reads, not 99 %. By hand, 2.2 % are faults, at 3 places of 132, and each of those claims is wrong only at its far
  end. Claims of 5 kb or more: none of 79 failed.
- The exact parts hold at the stated level by count (0.9 % refuted), with faults of two kinds: a consensus of
  related contexts, and one or two bases at a gene boundary where the sample's gene differs from the database's.
- Every fault of a certain part sits where reads of another locus had come onto the path. That is the one place to
  repair.
- The check is able to see errors: of 7,844 planted structural errors it passed 6, of which 5 exist in the sample,
  and it refuted every planted wrong or missing base that it could look at.
- Against four other tools on the same data the certain parts are by far the most reliable output (2 % failed
  against 14 to 52 %), and the shortest.

**What it does not show.**

- It is one laptop's worth of samples: 132 places. The upper limit for the share of places with a failed claim is
  8.6 %.
- Two samples could not be run in 16 GB: B-A3 (activated sludge, 50 M pairs) and, among the seven later ones, a
  cattle-faeces sample of 34 M pairs. The read pass grows with the diversity of the sample, not only with its size.
- The long reads are the truth only where they are deep enough: a rare true context next to an abundant one can
  be called contradicted, and a wrong claim inside a repeat longer than 3 kb can pass (section 4).
- Low recall is not an error, but it is the tool's main weakness on real data: it vouches for 1 kb on a quarter of
  the gene sides and for 5 kb on one in twenty.

**Seven more samples (17 November).** For the prototypes below, nine further pairs were taken from the same
table by a rule written beforehand (`PROTOCOL.md`, addendum 5); seven gave a figure (mouse, red panda and giant panda
faeces, a second cheese, digester, waste-water and chicken sample). On two, 0.2.1 gave nothing: a soil sample with no
resistance gene that passes the gene step, and a cattle sample on which the read pass reached 9.6 GB and was stopped,
as on B-A3. For 0.2.1 these are held-out samples like the others: 30 of its 744 claims of 1 kb or more fail (4.0 %),
all at one place on the waste-water sample (the place of B-A1's fault, in a sample of the same plant), and none of
829 exact parts is refuted. **Over all nineteen samples: 64 of 1,871 claims of 1 kb or more fail (3.4 %, 95 % interval
2.6 to 4.3 %), at 6 of 198 places; 11 of 2,011 exact parts are refuted (0.5 %).**

**What follows.** The statement in the README was changed to the measured rate. Two prototypes of a next
generation (a base-level proposer and a separate verifier) were built and held out in turn, each frozen and its
measure declared before its samples were scored (`benchmarks/dev-tools/experiments/frozen3/NOTE`, `frozen4/NOTE`,
`proposer/RESULTS.md`):

| On its own held-out samples | Claims ≥ 1 kb: failed / checkable | Places failed | Gene sides with a confirmed claim ≥ 1 / 2 / 5 kb | 0.2.1 on the same samples |
|---|---:|---:|---:|---|
| second generation (six samples) | 20 / 886 = 2.3 % | 3 | 80 / 47 / 27 | 0 / 238; 46 / 23 / 15 |
| third generation (seven fresh samples) | 20 / 1,552 = 1.3 % (0.8-2.0 %) | 1 | 131 / 69 / 23 | 30 / 744 = 4.0 %; 91 / 35 / 3 |
| fourth round (fifteen fresh samples, 18 November) | **5 / 1,585 = 0.32 % (0.10-0.73 %)** | 2 spots (5 gene sides) | 232 / 160 / 71 | 32 / 1,088 = 2.9 %; 182 / 97 / 38 |

The second and the third missed the bar of 1 %. The third generation, which repairs the three faults of the second,
confirms a claim of 2 kb on twice as many gene sides as 0.2.1 and of 5 kb on eight times as many, and fails less
often than 0.2.1 on the same samples, but its one failed place counts twenty claims. That fault was read to the base
(a count of leaving reads that misses a read whose first difference lies 31 bp before the place).

**The fourth round (18 November) meets the bar.** The repair turned out to be wider than that one place: the third
generation's record program missed about a sixth of all reads that leave the path or come onto it, and its long
claims had leaned on the undercount. The fourth-round prototype records every such read and lets the read's mate say
whether the molecule was on the path (`benchmarks/dev-tools/experiments/f5/NOTE` for the path to it, `frozen5/NOTE`
for what was frozen). Before it was frozen it was put through a stress test on the thirteen samples that still had
their reads: four times the usual number of candidates (33,862), every claim of 1 kb or more judged by the long
reads; that test found a fault of the round's own first variant before any fresh sample could, and none in the
version that was frozen (4,412 claims).
It was then run once on eighteen fresh samples chosen by a rule written beforehand (`PROTOCOL.md`, addendum 6): eight
human-gut samples and a laboratory reactor from two projects not used before, and nine samples of projects used
before (waste water, chicken, a digester). Fifteen gave a figure; on three waste-water samples aw-flank 0.2.1's read
pass, which the prototype needs, went over the memory allowed.

- **Claims of 1 kb or more: 5 failed of 1,585 (0.32 %, exact 95 % interval 0.10 to 0.73 %).** The third generation,
  run unchanged on the same samples: 10 of 1,694 (0.59 %); 0.2.1: 32 of 1,088 (2.9 %). On the nine samples of the two
  new projects alone: 4 of 1,305 (0.31 %, 0.08 to 0.78 %).
- It confirms a claim of 1 / 2 / 5 kb on 232 / 160 / 71 gene sides where 0.2.1 does on 182 / 97 / 38, and 54 %, 40 % and
  18 % of the true molecules lie inside a certain claim of that length, against 37 %, 21 % and 9 %.
- The five failed claims are two spots, both read to the base. At the first (one claim, right of *erm*(F) on the
  chicken sample) the path takes 88 bases that a fifth of the molecules carry and then runs on along the sequence of
  the other four fifths; the rule that vouched there took a small difference of two large depths for evidence. At
  the second (four claims that cross one spot 428 bp right of *tet*(A) on a gut sample) the short reads do not cross
  a 44-bp tandem copy at all: no read in the archived data does, the long reads do, and the path joins the two sides
  without the copy on the strength of 2 reads. (This report first called the two spots one kind and named a rule on
  the share of reads as the repair. That reading of the second spot was wrong; the correction, with the counts, is in
  `benchmarks/dev-tools/experiments/frozen5/NOTE`.)
- Two things the pooled figure hides. Long claims cross more places: of the claims of 5 kb or more, 4 of 248 fail
  (1.6 %; all four at one spot). And claims shorter than 1 kb, which the declared measure leaves out, fail in 79 of
  3,670 (2.2 %); 36 of those are the gene-copy class (the sample's gene differs from the database sequence inside
  the gene, where the check's first window begins).

**What this adds up to.** For 0.2.1, the released tool, block D gives the same picture as before: over all 34 fresh
real samples with a figure, 96 of 2,959 claims of 1 kb or more fail (3.2 %, 2.6 to 4.0 %). The propose-and-verify
design has now passed its bar once, on samples it had not seen, with a quarter more gene sides than 0.2.1 at 1 kb
and 1.6 to 1.9 times as many at 2 and 5 kb. It is still a prototype outside the package (two C programs, a rule file and
shell scripts that reuse 0.2.1's read pass, and therefore its memory limit).

**After the fourth round: two repairs, in development only.** Both spots have a repair that costs nothing on the
samples at hand (`benchmarks/dev-tools/experiments/f6/NOTE`). For the first, the flow rule's bar now counts that
the supply of a locus that came on is known to about a tenth of itself. For the second, a claim stops at a place
that hardly any read spans, where no reads are seen leaving or coming on either, and where the reads that stop or
start there with a few disagreeing bases at their very end outnumber the reads that cross (about 50 against 2 at
that spot; 4 against 15 at a like hole on another sample, where the long reads bear the path out). On the eleven
seen samples that still have their reads this gives 841 claims of 1 kb or more with none failed and the same gene
sides as before (109 / 79 / 36; the fourth-round prototype: 844 claims, 5 failed). Under the stress test on the seven
block-D samples that still have reads (12,644 candidates, four times the usual number), the fourth-round prototype
fails at the same two spots and nowhere else (5 of 2,909 claims), and the repaired rules at none (0 of 2,899, the
same gene sides). These are development figures on samples the repairs were shaped on, not a held-out result; a
version with them needs fresh samples of its own before anything is claimed for it.

**The package, aw-flank 0.3.0, held out on block E (19 November).** 0.3.0 is the fourth-round prototype with the two
repairs above, built as one package whose read pass no longer needs the whole map of a sample in memory. It was
frozen and then run once on sixteen further samples chosen by a rule written beforehand (`PROTOCOL.md`, addendum 7):
the eleven unused samples of the Beijing gut cohort with 3 Gbp or more of long reads, and the five samples of
earlier blocks on which 0.2.1's read pass had never fit (cattle faeces, three waste-water samples, an activated
sludge). All sixteen gave a figure.

- **Claims of 1 kb or more: 22 failed of 3,016 (0.73 %, exact 95 % interval 0.46 to 1.10 %).** Of 2 kb or more: 20 of
  1,291; of 5 kb or more: 0 of 220. On the same eleven gut samples 0.2.1 fails in 7 of 1,254 (0.56 %), and the
  fourth-round prototype, unchanged, in 25 of 1,760 on the nine where it fits in memory (1.42 %; three of those 25,
  at two spots, do not fail in the package).
- It confirms a claim of 1 / 2 / 5 kb on 359 / 227 / 87 gene sides. On the eleven gut samples a certain part covers
  53 %, 35 % and 8 % of the true molecules with that much flank, against 28 %, 11 % and 2 % for 0.2.1; a candidate
  matches 83 %, 72 % and 44 % of them.
- The five large samples gave 389 claims and none failed. They ran in 51 to 297 s with a peak of 4.0 to 6.8 GB under
  `--max-mem 7`; on the activated sludge (49.8 million read pairs) the map was built in two parts.
- **The 22 failed claims are three gene sides on two samples, and one fault** (read to the reads,
  `benchmarks/dev-tools/experiments/pkg03/NOTE`). Right of *erm*(F) on one gut sample, 14 candidates take the same
  wrong arm at a fork 1,956 bp out: the flow rule ("more reads than everything that came on can supply") put the
  supply of a locus that had come on at 2,135 reads where the depth shows about 3,500, and vouched for an arm of
  2,900. Left of *dfrA7* and *sul1* on another (one spot in a cassette array, seen from two genes): supply put at
  52 reads where about 190 were on the path. In both the other locus comes on softly: its reads differ from the
  path in single bases before, few of them are seen coming on, and the rule trusts that count. This is the class of
  the fourth round's failure on the chicken sample, which the fifth round's change of one constant was meant to
  repair. It did not.
- Claims shorter than 1 kb, outside the declared measure: 229 of 5,250 fail the check (4.4 %). For 162 of them long
  reads hold the claimed flank itself and the verdict comes from inside the gene, where the sample's copy differs
  from the database sequence; 65, at 11 places, are flanks that the long reads do not show (1.2 %).

**What block E adds up to.** Against the bar written down beforehand the package is under 1 % with an interval that
reaches over it: met, not settled. The block is sixteen samples the package had not seen, but no new kind of sample:
every project in it had been used before, and eleven of the sixteen are one cohort. The fault it found is not new
either; it is the one fault class that has now survived two rounds. What the package does reliably: claims of 5 kb
or more (0 of 220 here; the fourth round's 4 of 248 were one spot, repaired), the large samples that 0.2.1 could
not hold, and 1.6 to 1.8 times the gene sides of 0.2.1 on the gut samples (325 / 209 / 84 against 208 / 115 / 49 at
1 / 2 / 5 kb).

**After block E: what limits the length of certain parts, and a next version in development** (19 November, evening;
`benchmarks/dev-tools/experiments/tol/NOTE`; development figures on samples that were read while developing, not
held-out results).
- A candidate matches 83 % of the gene-carrying long reads for 1 kb and 44 % for 5 kb on the gut samples, but only
  53 % and 8 % lie inside a certain part. On two thirds of the gene sides some candidate is borne out a kilobase or
  more beyond everything that is called certain.
- Read at the reads, the largest single cause is this: a gut sample holds several strains' versions of one stretch
  that differ in a base here and a few bases there. 0.3.0 sees a read leaving or coming on by a structural measure
  but counts depth and templates only for reads without any difference, so the other versions' reads count for
  nothing, are charged to the gene's own reads at every difference, and count as a foreign locus at every return. In
  a made-up case (one context carried by two strains at 60x and 40x that differ in a few bases every 150 bp) 0.3.0
  calls 260 and 448 bp certain where the whole 10 kb is one context.
- A second, tolerant measure of the same records was built (small differences count as the same path, only
  structural divergence as leaving), and probed with made-up cases designed to provoke wrong claims. As first built
  it failed 19 of 180 of them (reads that drift off the path at 7 to 12 % difference are never seen leaving); with a
  check of the budget by the depth and a guard for crowded paths it fails none of 270, where the strict measure
  (0.3.0's counting, with the repaired flow rule) fails 2. On the three block-E samples that still have their reads, with both measures and the repaired flow rule: 781
  claims of 1 kb or more and none failed, against 648 with 22 failed for 0.3.0; 47 % and 27 % of the true molecules
  with 1 and 2 kb of flank inside a certain part, against 44 % and 23 %.
- What that will not change: where the path is a minority branch behind a stretch shared with another locus that is
  longer than a read pair reaches (300 bp in this cohort's libraries), the short reads do not say which molecules
  take it. The best-supported branch at every fork (the context marked `main`) is right at 1 / 2 / 5 kb on 93 / 85 /
  68 % of gene sides; what is called certain is right in over 99 % of claims. That difference is what the certain
  part is for.

**The package, aw-flank 0.4.0, held out on block F (21–22 November).** `0.4.0` adds the repaired flow rule (`ZR = 4.0`,
far-read counts, budget cap), single-base soft-leaver charges, inside-gene allele vs foreign classification over 300 bp,
post-gene leaver charging, a 20 % uneven-depth budget allowance, tolerant proposing and dual strict + tolerant counting,
100-bp bridging with 300-bp sustain, and the end-anchored second claim (`certain_end_bp`). Frozen on 21 November
(`benchmarks/dev-tools/experiments/pkg04/FROZEN`) alongside paired frozen `0.3.0`, it was run once on thirteen completely
unseen samples from three projects never used before (`PROTOCOL.md`, addendum 8; `results/blockF/`): eleven Malawian
children's gut samples (`PRJNA1139951`, Illumina TELL-seq barcode-stripped short reads of ~2×114–131 bp vs PacBio Revio
HiFi) and two 2×151 bp NovaSeq anaerobic digester / stillage bioreactor samples (`PRJNA768577`, `PRJNA708972` vs Sequel
II/IIe HiFi). All thirteen ran to completion (`27` to `411 s` and `1.0` to `7.0 GB` on twelve samples; `1,629 s` and
`1.2 GB` on `F-J2`).

- **Unfiltered mechanical scorer across all thirteen samples:** `73 / 1,498` (`4.87 %`, `103` uncheckable) for `0.4.0`
  primary (`p4`), `73 / 1,553` (`4.70 %`, `118` uncheckable) for `0.4.0` end-anchored (`p4e`), and `92 / 1,852` (`4.97 %`,
  `135` uncheckable) for paired `0.3.0` (`p3`). Every one of the `73` raw contradictions in `0.4.0` comes from two
  cross-library artifacts outside the algorithm, verified base by base (`benchmarks/dev-tools/experiments/pkg04/NOTE`):
  1. **SRA short-read sample mismatch on `F-M3` (`SAMN42799965`, `Hsap.meta.a4`; `34 / 71` contradicted and `88 / 159`
     uncheckable in `0.4.0`; `33 / 76` contradicted and `119 / 195` uncheckable in `0.3.0`):** in `SRR29981012` (Illumina
     TELL-seq), `816 / 841` single reads spanning the `erm(F)` 5' boundary carry `AATTCAGCCCTGTAAGAAGTTACTA` (`0` carry
     `AAATTTTAACGAACCATTGTTACTA`), `310 / 344` 3' reads carry `...GCATAACGGTTTGGCGC`, `121 / 122` `tet(Q)` left reads
     match `aw-flank`'s junction across `-1,114 bp`, and `dfrA14` is at `283×` depth; whereas both PacBio Revio HiFi
     `SRR29980953` (`2,235 / 2,276` 5' reads) **and** the third independent run of `SAMN42799965`, Oxford Nanopore
     PromethION `SRR29980994` (`20 / 20` sampled 5' reads), carry `AAATTTTAACGAACCATTGTTACTA` (`0` carry the `SRR29981012`
     sequence) and `0` reads of `dfrA14`.
  2. **TELL-seq recombinant reagent plasmid across the remaining samples (`39 / 39` non-`F-M3` contradictions in `0.4.0`;
     `56 / 59` in `0.3.0`):** all `39` claims are `blaTEM-116` / `blaTEM-1B` reconstructing the exact same `4,724-bp`
     circular `pMB1`/`ColE1` recombinant expression vector (`861-bp` `blaTEM-116` CDS + `3,863-bp` `pMB1 ori` backbone,
     stopping as a loop at `4,354 bp`), present at `175–1,279×` depth in **11 of 11** Malawian TELL-seq short-read
     libraries (from *E. coli*-expressed TELL-seq transposase/enzyme reagents) and **0 of 2** non-TELL-seq NovaSeq
     libraries (`F-J1`, `F-J2`), and absent from the unamplified stool HiFi libraries (whose native `Tn2`/`Tn3`
     `blaTEM-1B` copies diverge at `-201 bp` left and `+149 bp` right).
- **Accuracy on native sample DNA across all twelve concordant unseen samples:** **`0 / 1,388` primary claims ≥ 1 kb
  failed (`0.00 %`, exact 95 % interval `0.00 to 0.27 %`)** and **`0 / 1,415` end-anchored claims ≥ 1 kb failed (`0.00 %`,
  `0.00 to 0.26 %`)**, with `0 / 480` ≥ 2 kb, `0 / 167` ≥ 5 kb, and `139 / 62 / 22` (`141 / 63 / 23` end-anchored)
  gene-family sides confirmed at `1 / 2 / 5 kb`. That includes **`0 / 1,285`** on the ten variable-length ~2×114–131 bp
  TELL-seq gut samples (first held-out test of `0.4.0`'s scaled read-length constants) and **`0 / 103`** on the two
  2×151 bp JGI digester/bioreactor samples (`11` vs `7` confirmed gene-family sides ≥ 1 kb, `48.9 %` vs `44.5 %` 1-kb
  claimed molecule recall).
- **Against paired `0.3.0` on the same twelve samples:** `0.3.0` had **`3 / 1,702` native claims ≥ 1 kb failed (`0.18 %`)**
  — `F-M7` `erm(F)` left `c1` (`4,394 bp`, contradicted by `305 / 324` HiFi reads) and `F-M9` `erm(F)` left `c1, c2`
  (`2,824` and `2,819 bp`, contradicted by `487 / 497` HiFi reads), both instances of Block E's soft-join flow-rule fault
  class where `0.4.0` stopped cleanly before the fork (`1,593 bp` and `1,201 bp`) — and `90 / 2,881` (`3.12 %`) vs
  `34 / 2,800` (`1.21 %`) failed on shorter native claims (`100–999 bp`; `3 / 2,244 = 0.13 %` vs `36 / 2,251 = 1.60 %` on
  the ten gut samples). Candidate coverage of true long-read molecules rose from `78.4 / 67.0 / 42.2 %` (`0.3.0`) to
  `80.1 / 69.5 / 43.5 %` (`0.4.0`) at `1 / 2 / 5 kb`.
