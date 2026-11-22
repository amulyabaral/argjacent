#!/usr/bin/env python3
"""Self-test of ARGjacent on small simulated data with a known answer (about 15 s).

Six cases, each a 1 kb gene with 10 kb flanks, 2x150 bp read pairs (fragment ~400 bp, 0.3% substitution errors):

  clean   120 kb genome, no repeats near the gene, 40x      -> one context per side, certain over the whole 10 kb
  repeat  the same genome with 5 copies of a 1.3 kb repeat, one of them 3 kb downstream of the gene
                                                             -> right side certain up to the end of the repeat, then candidates
  mixed   the gene in two genomes, at 80x and at 50x         -> both contexts on both sides, each certain for itself
  strains the gene in one context that two strains carry (60x and 40x); the second strain's version differs in a
          cluster of 3 to 5 bases about every 150 bp and in single bases between
                                                             -> the context certain over the whole 10 kb on both sides
  alleles the gene in one context, as given (50x) and as another allele (30x: 2 % of its bases differ, a cluster of
          differing bases and 9 bases fewer near either end)   -> the context certain over the whole 10 kb on both sides
  shared  the clean genome (40x), and a second genome (60x) without the gene in which another locus holds a copy of
          the gene's last 160 bases                           -> the gene's own context certain over the whole 10 kb; the
                                                                flank of the other locus is proposed but not certain

For every context reported, the part given as certain must follow a true context base for base; in the case "strains",
where a path may carry one strain's bases here and the other's there, it must follow a true context but for small
differences (no 300 bases with more than 25 differing). (The last 30 bases of what the reads vouch for are never called
certain, so a full flank of 10,000 bp is certain for 9,970.)

Usage: argjacent-selftest [out_dir]      (default: a temporary folder, removed afterwards)
"""
import json
import os
import random
import shutil
import sys
import tempfile

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def simulate(out, seed=42):
    """write gene.fa, truth.json and the read pairs of the six cases into `out`"""
    R = random.Random(seed)

    def rand_seq(n):
        return "".join(R.choice("ACGT") for _ in range(n))

    def mutate(s, rate):
        o = list(s)
        for i in range(len(o)):
            if R.random() < rate:
                o[i] = R.choice([b for b in "ACGT" if b != o[i]])
        return "".join(o)

    def sim_pairs(genomes_cov, prefix, rlen=150, err=0.003):
        q = "I" * rlen
        with open(prefix + "_1.fq", "w") as f1, open(prefix + "_2.fq", "w") as f2:
            for gi, (g, cov) in enumerate(genomes_cov):
                L, gg = len(g), g + g[:2000]
                for k in range(int(cov * L / (2 * rlen))):
                    ins = max(rlen + 20, int(R.gauss(400, 40)))
                    start = R.randrange(0, L)
                    frag = gg[start:start + ins]
                    if R.random() < 0.5:
                        frag = rc(frag)
                    f1.write("@g%d_%d/1\n%s\n+\n%s\n" % (gi, k, mutate(frag[:rlen], err), q))
                    f2.write("@g%d_%d/2\n%s\n+\n%s\n" % (gi, k, mutate(rc(frag)[:rlen], err), q))

    os.makedirs(out, exist_ok=True)
    gene, IS = rand_seq(1000), rand_seq(1300)
    A = rand_seq(50000) + gene + rand_seq(69000)
    gp = 50000
    Arep = A
    for pos, ori in [(105000, -1), (85000, 1), (gp + 1000 + 3000, 1), (25000, -1), (10000, 1)]:
        Arep = Arep[:pos] + (IS if ori == 1 else rc(IS)) + Arep[pos:]
    gpr = Arep.find(gene)
    B = rand_seq(60000) + gene + rand_seq(59000)
    gpb = 60000
    # another strain's version of genome A: the same sequence, with a cluster of 3 to 5 differing bases within 12 bp about every
    # 150 bp and single differing bases between (substitutions only, so that positions stay comparable); the gene itself is the same
    R2 = random.Random(seed + 1)
    a2 = list(A)
    pos = 60
    while pos < len(a2) - 20:
        if not (gp - 5 <= pos <= gp + 1005):
            if R2.random() < 0.75:
                for q in R2.sample(range(pos, pos + 12), R2.randint(3, 5)):
                    a2[q] = R2.choice([b for b in "ACGT" if b != a2[q]])
            else:
                a2[pos] = R2.choice([b for b in "ACGT" if b != a2[pos]])
        pos += R2.randint(110, 190)
    A2 = "".join(a2)
    # another allele of the gene in the same context: 2 % single differences, a cluster of 4 within 12 bp 250 bp from either end,
    # and 9 bases fewer 150 bp from either end
    R3 = random.Random(seed + 2)
    g2 = list(gene)
    for i in range(len(g2)):
        if R3.random() < 0.02:
            g2[i] = R3.choice([b for b in "ACGT" if b != g2[i]])
    for c0_ in (250, 740):
        for q in R3.sample(range(c0_, c0_ + 12), 4):
            g2[q] = R3.choice([b for b in "ACGT" if b != g2[q]])
    gene2 = "".join(g2[:150] + g2[159:841] + g2[850:])
    A3 = A[:gp] + gene2 + A[gp + 1000:]
    # a genome without the gene in which another locus holds a copy of the gene's last 160 bases
    Cg = rand_seq(40000) + gene[840:] + rand_seq(40000)
    truth = {
        "alleles": {"A": {"left": A[gp - 10000:gp], "right": A[gp + 1000:gp + 11000]}},
        "shared": {"A": {"left": A[gp - 10000:gp], "right": A[gp + 1000:gp + 11000]},
                   "other": {"left": "", "right": Cg[40160:50160]}},
        "strains": {"A": {"left": A[gp - 10000:gp], "right": A[gp + 1000:gp + 11000]},
                    "A2": {"left": A2[gp - 10000:gp], "right": A2[gp + 1000:gp + 11000]}},
        "clean": {"A": {"left": A[gp - 10000:gp], "right": A[gp + 1000:gp + 11000]}},
        "repeat": {"A": {"left": Arep[gpr - 10000:gpr], "right": Arep[gpr + 1000:gpr + 11000]}},
        "mixed": {"A": {"left": A[gp - 10000:gp], "right": A[gp + 1000:gp + 11000]},
                  "B": {"left": B[gpb - 10000:gpb], "right": B[gpb + 1000:gpb + 11000]}},
    }
    with open(os.path.join(out, "gene.fa"), "w") as f:
        f.write(">demo_gene\n")
        for i in range(0, len(gene), 80):
            f.write(gene[i:i + 80] + "\n")
    with open(os.path.join(out, "truth.json"), "w") as f:
        json.dump(truth, f)
    sim_pairs([(A, 40)], os.path.join(out, "clean"))
    sim_pairs([(Arep, 40)], os.path.join(out, "repeat"))
    sim_pairs([(A, 80), (B, 50)], os.path.join(out, "mixed"))
    sim_pairs([(A, 60), (A2, 40)], os.path.join(out, "strains"))
    sim_pairs([(A, 50), (A3, 30)], os.path.join(out, "alleles"))
    sim_pairs([(A, 40), (Cg, 60)], os.path.join(out, "shared"))


def _lcp(a, b):
    n = 0
    while n < min(len(a), len(b)) and a[n] == b[n]:
        n += 1
    return n


def _follows(c, truths, n):
    """does the first n bases of c follow one of the truths but for small differences: no 300 bases with more than 25 differing
    from EVERY truth at once (a path may change from one strain's bases to the other's)"""
    bad = [0] * (n + 1)
    for i in range(n):
        bad[i + 1] = bad[i] + (0 if any(i < len(t) and c[i] == t[i] for t in truths) else 1)
    return all(bad[min(n, i + 300)] - bad[i] <= 25 for i in range(0, max(1, n - 299)))


EXPECT = {"strains": "one context of two strains: certain for the whole flank on both sides, in at least one candidate",
          "clean": "one context per side, certain for the whole flank",
          "repeat": "left as clean; right certain up to the end of the repeat, the true path among the candidates",
          "mixed": "both genomes' contexts on both sides, each certain for itself",
          "alleles": "two alleles of the gene in one context: certain for the whole flank on both sides",
          "shared": "the gene's context certain for the whole flank; the other locus's flank proposed, not certain as a context of the gene, certain as a context of a gene that ends like it"}
WHOLE = 9900            # a flank of 10,000 bp that is certain throughout is given as certain for 9,970
CASES = ("clean", "repeat", "mixed", "strains", "alleles", "shared")


def check(out, echo=print):
    """compare the results in `out` with the known answer; True if all is as expected"""
    truth = json.load(open(os.path.join(out, "truth.json")))
    ok_all = True
    for case in CASES:
        res = json.load(open(os.path.join(out, case + ".flanks.json")))["genes"]["demo_gene"]
        echo("%-6s expected: %s" % (case, EXPECT[case]))
        for side in ("left", "right"):
            paths = [p for p in res[side]["paths"] if p["flank"]]
            truths = {k: (v[side][::-1] if side == "left" else v[side]) for k, v in truth[case].items()}
            cands = [((p["flank"][::-1] if side == "left" else p["flank"]), p) for p in paths]
            certain_ok = all(any(_lcp(c, t) >= p["context_len"] for t in truths.values()) for c, p in cands)
            best = {k: max((min(p["context_len"], _lcp(c, t)) for c, p in cands), default=0) for k, t in truths.items()}
            full = {k: any(c == t for c, _ in cands) for k, t in truths.items()}
            if case == "strains":
                # every certain part follows the two strains' context but for small differences, and one is certain throughout
                certain_ok = all(_follows(c, list(truths.values()), p["context_len"]) for c, p in cands)
                longest = max((p["context_len"] for c, p in cands), default=0)
                good = certain_ok and longest >= WHOLE
                ok_all &= good
                echo("   %-5s %-4s  contexts %d | longest certain part %d bp%s" % (side, "OK" if good else "FAIL", len(cands), longest,
                                                                              "" if certain_ok else " (a certain part does not follow the true context)"))
                continue
            if case == "alleles":
                good = certain_ok and best["A"] >= WHOLE
                ok_all &= good
                echo("   %-5s %-4s  contexts %d | true context certain %d bp" % (side, "OK" if good else "FAIL", len(cands), best["A"]))
                continue
            if case == "shared":
                # a context that is not the gene's own (it follows the other locus, or nothing) must not be certain at all
                mine = [(c, p) for c, p in cands if _lcp(c, truths["A"]) >= min(200, len(c))]
                others = [p["context_len"] for c, p in cands if _lcp(c, truths["A"]) < min(200, len(c))]
                proposed = side == "left" or any(_lcp(c, truths["other"]) >= 200 for c, p in cands)
                good = best["A"] >= WHOLE and all(_lcp(c, truths["A"]) >= p["context_len"] for c, p in mine) and max(others, default=0) == 0 and proposed
                # ... but it is a context of a gene that ENDS like this one: the end claim of the other locus's flank holds as far
                # as that flank goes, and says how many of the gene's last bases the two loci share (160)
                endnote = ""
                if side == "right":
                    oth = [(c, p) for c, p in cands if _lcp(c, truths["other"]) >= 200]
                    e_ok = bool(oth) and all(_lcp(c, truths["other"]) >= p["end_len"] for c, p in oth) \
                        and max(p["end_len"] for c, p in oth) >= WHOLE and all(140 <= p["shared_end"] <= 180 for c, p in oth if p["end_len"])
                    good = good and e_ok
                    endnote = "; as a context of a gene that ends like it: %d bp (last %s bases shared)" % (
                        max((p["end_len"] for c, p in oth), default=0), "/".join(sorted({str(p["shared_end"]) for c, p in oth})) or "-")
                ok_all &= good
                echo("   %-5s %-4s  contexts %d | true context certain %d bp; contexts that are not the gene's: %d, certain for at most %d bp%s" % (
                    side, "OK" if good else "FAIL", len(cands), best["A"], len(others), max(others, default=0), ("" if proposed else " (the other locus's flank was not proposed)") + endnote))
                continue
            if case == "clean" or (case == "repeat" and side == "left"):
                good = certain_ok and len(cands) == 1 and best["A"] >= WHOLE
            elif case == "repeat":
                good = certain_ok and 2900 <= best["A"] <= 4300 and full["A"]
            else:
                good = certain_ok and best["A"] >= WHOLE and best["B"] >= WHOLE and len(cands) == 2
            ok_all &= good
            note = ", ".join("%s certain %d bp%s" % (k, best[k], " (whole 10 kb among the candidates)" if full[k] else "")
                             for k in truths)
            echo("   %-5s %-4s  contexts %d | true context: %s" % (side, "OK" if good else "FAIL", len(cands), note))
    return ok_all


def main(argv=None):
    from . import cli
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] in ("-h", "--help"):
        print(__doc__)
        return 0
    out = argv[0] if argv else tempfile.mkdtemp(prefix="argjacent-selftest-")
    print("== simulating six small test cases in %s" % out)
    simulate(out)
    for c in CASES:
        print("\n== case: %s" % c, flush=True)
        rc_ = cli.main(["-g", os.path.join(out, "gene.fa"), "-1", os.path.join(out, c + "_1.fq"), "-2", os.path.join(out, c + "_2.fq"),
                        "-o", os.path.join(out, c), "-t", "2"])
        if rc_:
            return rc_
    print("\n== checking against the known answer")
    ok = check(out)
    print("All checks passed." if ok else "SOME CHECKS FAILED - please report the output above.")
    if not argv:
        shutil.rmtree(out, ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
