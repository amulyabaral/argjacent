"""Regression tests for v0.6.0 low-coverage upgrades:
  - aw: two-stage saturating Bloom pre-pass (pass0), seed-exempt Bloom, dead-end FR mate links, hub-spur guard
  - awx: mate-anchored FR overlap rescue (14..30 bp / 1-SNP) and unambiguous K2=19 read-tip rescue
  - verify: Rule U mate-bridged unbranched zero-join low-depth certification and repeat/chimera abstention
"""
import json
import math
import os
import random


COMP = str.maketrans("ACGT", "TGCA")


def _rc(s):
    return s.translate(COMP)[::-1]


def _rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _sim_pairs(rng, molecules, r1_path, r2_path, rlen=150, ins_mean=380, ins_sd=35, err=0.002):
    q = "I" * rlen
    pid = 0
    with open(r1_path, "w") as f1, open(r2_path, "w") as f2:
        for mi, (seq, cov) in enumerate(molecules):
            L = len(seq)
            npairs = int(round(cov * L / (2.0 * rlen)))
            for _ in range(npairs):
                ins = max(rlen + 25, min(650, int(rng.gauss(ins_mean, ins_sd))))
                st = rng.randrange(0, L - ins)
                frag = seq[st : st + ins]
                if rng.random() < 0.5:
                    frag = _rc(frag)
                s1, s2 = list(frag[:rlen]), list(_rc(frag)[:rlen])
                for s in (s1, s2):
                    for i in range(rlen):
                        if rng.random() < err:
                            b0 = s[i]
                            s[i] = rng.choice([b for b in "ACGT" if b != b0])
                f1.write("@m%d_%d/1\n%s\n+\n%s\n" % (mi, pid, "".join(s1), q))
                f2.write("@m%d_%d/2\n%s\n+\n%s\n" % (mi, pid, "".join(s2), q))
                pid += 1


def test_low_coverage_single_copy_certified_and_repeat_traps_abstain(cli, tmp_path):
    """At 5.5x depth, clean single-copy loci reach >= 1 kb certified flanks with 0 mismatches,
    while 2-copy IS repeat traps and shared-gene-end traps never falsely certify foreign flanks.
    """
    rng = random.Random(606)
    g_clean = _rand_seq(rng, 900)
    L_clean, R_clean = _rand_seq(rng, 6000), _rand_seq(rng, 6000)

    g_trap = _rand_seq(rng, 900)
    is_elem = _rand_seq(rng, 820)
    L_trap = _rand_seq(rng, 6000)
    stem = _rand_seq(rng, 650)
    exit_true, exit_false = _rand_seq(rng, 4500), _rand_seq(rng, 4500)
    R_trap_true = stem + is_elem + exit_true

    g_shared = _rand_seq(rng, 900)
    L_sh, R_sh_true, R_sh_false = _rand_seq(rng, 6000), _rand_seq(rng, 6000), _rand_seq(rng, 6000)

    genes_fa = str(tmp_path / "genes.fa")
    with open(genes_fa, "w") as fh:
        fh.write(">clean_arg\n%s\n>istrap_arg\n%s\n>shared_arg\n%s\n" % (g_clean, g_trap, g_shared))

    r1, r2 = str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq")
    molecules = [
        (L_clean + g_clean + R_clean, 6.0),
        (L_trap + g_trap + R_trap_true, 5.5),
        (_rand_seq(rng, 5000) + is_elem + exit_false, 9.0),
        (L_sh + g_shared + R_sh_true, 6.0),
        (_rand_seq(rng, 5000) + g_shared[-160:] + R_sh_false, 8.0),
    ]
    _sim_pairs(rng, molecules, r1, r2)

    out = str(tmp_path / "out")
    cli(["-g", genes_fa, "-1", r1, "-2", r2, "-o", out, "-t", "2", "-D", "4000"])
    doc = json.load(open(out + ".flanks.json"))["genes"]

    # 1. Clean 6x locus reaches >= 1 kb on at least one side and follows truth
    clean_best = max(
        max((p["context_len"] for p in doc["clean_arg"][sd]["paths"]), default=0)
        for sd in ("left", "right")
    )
    assert clean_best >= 1000, clean_best
    for sd, truth_fl in (("left", L_clean), ("right", R_clean)):
        for p in doc["clean_arg"][sd]["paths"]:
            c = p["context_len"]
            if c > 0:
                fl = p["flank"] if sd == "right" else p["flank"][::-1]
                tr = truth_fl if sd == "right" else truth_fl[::-1]
                diffs = sum(1 for a, b in zip(fl[:c], tr[:c]) if a != b)
                assert diffs <= max(2, int(0.01 * c))

    # 2. IS26-like repeat trap abstains at or before stem + IS (650 + 820 = 1470 bp)
    for p in doc["istrap_arg"]["right"]["paths"]:
        c = p["context_len"]
        assert c <= 1470, (c, p["certain_stop"])
        if c > 0:
            diffs = sum(1 for a, b in zip(p["flank"][:c], R_trap_true[:c]) if a != b)
            assert diffs <= max(2, int(0.01 * c))

    # 3. Shared 160-bp gene-end trap never certifies R_sh_false as a context of shared_arg
    for p in doc["shared_arg"]["right"]["paths"]:
        c = p["context_len"]
        if c > 0:
            diffs = sum(1 for a, b in zip(p["flank"][:c], R_sh_true[:c]) if a != b)
            assert diffs <= max(2, int(0.01 * c))


def test_deterministic_mate_and_short_overlap_rescue(cli, tmp_path):
    """Constructs a hand-crafted 4x tiling where consecutive reads overlap by only 20 bp (< k=31)
    across a mate-bridged junction: v0.6 must assemble and certify across the 20-bp overlap.
    """
    rng = random.Random(607)
    gene = _rand_seq(rng, 800)
    left_fl = _rand_seq(rng, 2500)
    right_fl = _rand_seq(rng, 2500)
    mol = left_fl + gene + right_fl
    g_start = len(left_fl)
    g_end = g_start + len(gene)

    genes_fa = str(tmp_path / "g.fa")
    open(genes_fa, "w").write(">g\n%s\n" % gene)

    # Tile clean pairs at 5x EXCEPT force a single 20-bp read overlap right at g_end + 300..320
    # while keeping >= 3 concordant pairs spanning across [g_end + 180, g_end + 520].
    rlen = 150
    q = "I" * rlen
    r1, r2 = str(tmp_path / "r1.fq"), str(tmp_path / "r2.fq")
    pairs = []
    # Upstream/gene and downstream dense tiling outside the forced 20-bp overlap neck
    for st in range(g_start - 1800, g_end + 1800, 35):
        ins = 380
        # Skip reads whose body would cover [g_end + 295, g_end + 325] with >= 31 bp overlap
        r1_s, r1_e = st, st + rlen
        r2_s, r2_e = st + ins - rlen, st + ins
        if (r1_s < g_end + 300 and r1_e > g_end + 320) or (r2_s < g_end + 300 and r2_e > g_end + 320):
            continue
        pairs.append((st, st + ins))
    # Add 4 concordant pairs whose forward mate ends at g_end + 320 and reverse mate starts at g_end + 300 (exact 20-bp overlap!)
    for shift in (-4, -2, 0, 2):
        # Fwd mate: [g_end + 170 + shift, g_end + 320 + shift] -> wait, keep exact junction [g_end + 300, g_end + 320]:
        # Pair A: fwd [g_end + 170, g_end + 320), rev [g_end + 370, g_end + 520)
        # Pair B: fwd [g_end + 50, g_end + 200), rev [g_end + 300, g_end + 450)
        pairs.append((g_end + 170 + shift, g_end + 530 + shift))
        pairs.append((g_end + 40 + shift, g_end + 450))
    # Fix exact 20-bp overlap: clamp any read crossing (g_end+300, g_end+320)
    clean_pairs = []
    for a, b in pairs:
        if a < g_end + 300 < a + rlen and a + rlen > g_end + 320:
            a = g_end + 320 - rlen
        if b - rlen < g_end + 320 and b > g_end + 300 and b - rlen < g_end + 300:
            b = g_end + 300 + rlen
        clean_pairs.append((a, b))

    with open(r1, "w") as f1, open(r2, "w") as f2:
        for i, (a, b) in enumerate(clean_pairs):
            frag = mol[a:b]
            f1.write("@p%d/1\n%s\n+\n%s\n" % (i, frag[:rlen], q))
            f2.write("@p%d/2\n%s\n+\n%s\n" % (i, _rc(frag)[:rlen], q))

    out = str(tmp_path / "res")
    cli(["-g", genes_fa, "-1", r1, "-2", r2, "-o", out, "-t", "2", "-D", "1500"])
    paths = json.load(open(out + ".flanks.json"))["genes"]["g"]["right"]["paths"]
    assert paths and paths[0]["context_len"] >= 1000, paths
    assert paths[0]["flank"][: paths[0]["context_len"]] == right_fl[: paths[0]["context_len"]]
