"""The tolerant measure of the record program (awv -T 25 -B 1 -I 12) and what it is for: several versions of one stretch."""
import json
import os
import random
import subprocess

import numpy as np

from argjacent import cli as awcli
from argjacent import pipeline, verify

COMP = str.maketrans("ACGT", "TGCA")
G = 400                                                      # the gene part of the made-up candidate


def _molecules(tmp_path):
    """one candidate of 2,400 bp and read pairs of five kinds of molecule: the path itself (40 deep); a version with 5
    differing bases within 20 bp at flank 400 (40); a version that lacks 3 bases at flank 800 (30); a molecule that follows
    the path to flank 1200 and then goes elsewhere (30); a foreign molecule that comes on at flank 1500 (30)"""
    R = random.Random(7)

    def rnd(n):
        return "".join(R.choice("ACGT") for _ in range(n))

    def mut(s, pos):
        s = list(s)
        for p in pos:
            s[p] = R.choice([b for b in "ACGT" if b != s[p]])
        return "".join(s)
    path, pre = rnd(2400), rnd(600)
    mols = [(pre + path + rnd(300), 40),
            (pre + mut(path, [G + 400, G + 404, G + 409, G + 413, G + 419]) + rnd(300), 40),
            (pre + path[:G + 800] + path[G + 803:] + rnd(300), 30),
            (pre + path[:G + 1200] + rnd(1500), 30),
            (rnd(600 + G + 1500) + path[G + 1500:] + rnd(300), 30)]
    reads, cands = str(tmp_path / "reads.tsv"), str(tmp_path / "cands.fa")
    with open(reads, "w") as f:
        pid = 0
        for m, depth in mols:
            for _ in range(int(depth * len(m) / 300)):
                ins = max(160, int(R.gauss(320, 40)))
                a = R.randrange(0, len(m) - ins)
                f.write("%d\t1\t%s\t\t%s\n%d\t2\t%s\t\t%s\n" % (pid, m[a:a + 150], "I" * 150, pid, m[a + ins - 150:a + ins].translate(COMP)[::-1], "I" * 150))
                pid += 1
    open(cands, "w").write(">0|right|1|%d|0|0\n%s\n" % (G, path))
    return reads, cands, len(path)


def _records(tmp_path, reads, cands, n, tag, options):
    jobs, pre = str(tmp_path / ("jobs." + tag)), str(tmp_path / tag)
    open(jobs, "w").write("%s\t%s\n" % (cands, pre))
    subprocess.run([awcli.find_bin("awv")] + options + [jobs, reads], check=True, capture_output=True)
    A = np.fromfile(pre + ".i32", dtype=np.int32).reshape(verify.NARR, n + 1)
    ent = np.fromfile(pre + ".ent", dtype=np.int32).reshape(-1, 8)
    lev = np.fromfile(pre + ".lev", dtype=np.int32).reshape(-1, 8)
    return pre, A, ent, lev


def _near(rows, col, lo, hi):
    r = rows[(rows[:, col] >= G + lo) & (rows[:, col] <= G + hi)]
    return [int((r[:, 4] % 10000 // 1000 == k).sum()) for k in (0, 1, 2)]      # by the event's class (tolerant measure; strict: all 0)


def _marks(rows, lo, hi):
    """reads that come on between lo and hi, by what they are like before the place: not told, like the path, foreign"""
    r = rows[(rows[:, 1] >= G + lo) & (rows[:, 1] <= G + hi)]
    return [int((r[:, 4] // 10000 == k).sum()) for k in (0, 1, 2)]


def test_strict_and_tolerant_records_of_five_kinds_of_molecule(tmp_path):
    reads, cands, n = _molecules(tmp_path)
    _, S, sent, slev = _records(tmp_path, reads, cands, n, "strict", [])
    pre, T, tent, tlev = _records(tmp_path, reads, cands, n, "tolerant", pipeline.TOLERANT)
    d = verify.A_EXACT
    # strict: inside the patch only the path's own version counts; the other version's reads leave before it and come on after it
    assert S[d][G + 410] < S[d][G + 200] - 20 and _near(slev, 3, 390, 425)[0] >= 15 and _near(sent, 1, 400, 440)[0] >= 15
    # tolerant: they count through the patch, and none of them is an event of class 0; those that pass are counted (class 1)
    assert T[d][G + 410] > S[d][G + 410] + 20
    assert _near(tlev, 3, 390, 425)[0] == 0 and _near(tent, 1, 400, 440)[0] == 0
    assert _near(tlev, 3, 390, 425)[1] >= 15 and _near(tent, 1, 400, 440)[1] >= 15
    # the version that lacks three bases counts across the place
    assert T[d][G + 801] > S[d][G + 801] + 15
    # the fork and the foreign molecule are seen as before
    assert _near(tlev, 3, 1190, 1215)[0] >= 15 and _near(slev, 3, 1190, 1215)[0] == _near(tlev, 3, 1190, 1215)[0]
    assert _near(tent, 1, 1490, 1515)[0] >= 10
    # what a read is like before the place where it comes on: the other version's reads after the patch are like the path there,
    # the foreign molecule's are not
    m_patch, m_foreign = _marks(sent, 400, 440), _marks(sent, 1490, 1515)
    assert m_patch[1] >= 10 and m_patch[2] == 0 and m_foreign[2] >= 8 and m_foreign[1] == 0
    # the windows of 31 and 61 bases do not run on over the patch: the path's own version must be carried by reads
    assert T[verify.A_S61][G + 430] == S[verify.A_S61][G + 430]
    # and the rules read both
    assert verify.verify_gene(str(tmp_path / "strict"))[0]["verified"] > 1000
    assert verify.verify_gene(pre, patches=True)[0]["verified"] > 1000


def test_one_context_in_two_strains_is_certain_throughout(data, cli, tmp_path):
    """the self-test's fourth case: under the strict measure alone the certain part ends within a few hundred bases (the other
    strain's reads leave and come on at every cluster of differing bases); with both measures the whole flank is certain"""
    out = str(tmp_path / "s")
    r = subprocess.run([os.sys.executable, "-m", "argjacent", "-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "strains_1.fq"),
                        "-2", os.path.join(data, "strains_2.fq"), "-o", out, "-t", "2", "--keep"],
                       capture_output=True, text=True, env=dict(os.environ, AW_KEEP_RECORDS="1"))
    assert r.returncode == 0, r.stdout + r.stderr
    res = json.load(open(out + ".flanks.json"))["genes"]["demo_gene"]
    for side in ("left", "right"):
        assert max(p["context_len"] for p in res[side]["paths"]) >= 9900
    strict = verify.verify_gene(out + ".tmp/v/0")
    tolerant = verify.verify_gene(out + ".tmp/vt/0", patches=True)
    assert max(m["verified"] for m in strict) < 1000 < 9900 < max(m["verified"] for m in tolerant)
