"""The rules that say how far a candidate is certain."""
import importlib.util
import json
import os

import pytest

from argjacent import verify


def test_limits_of_a_count():
    assert verify.ucb(0) == 3.0 and verify.lcb(0) == 0.0
    for n in range(0, 60):
        assert verify.lcb(n) <= n <= verify.ucb(n)
        assert verify.ucb(n + 1) > verify.ucb(n) and verify.lcb(n + 1) >= verify.lcb(n)
    assert abs(verify.ucbf(4.0) - verify.ucb(4)) < 1e-9 and verify.ucb(4) < verify.ucbf(4.5) < verify.ucb(5)


def test_read_length_terms():
    c150 = verify.read_length_constants(150)
    assert c150 == {"L": 150, "eff": 0.65, "show": 98.0, "tail": 138.0, "stretch_max": 125, "far": 70}        # the constants the rules were shaped with
    assert verify.read_length_constants(151) == c150 and verify.read_length_constants(251) == c150   # longer reads: those of 150 bp
    c100 = verify.read_length_constants(100)
    assert c100["L"] == 100 and c100["eff"] == 0.48 and c100["show"] == 48.0 and c100["far"] == 47
    assert verify.read_length_constants(36)["L"] == 60


def _reference():
    """the rule file of the prototype as a module (made by benchmarks/dev-tools/rules-reference/port_vcall.py), if this is a source folder"""
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "benchmarks", "dev-tools", "rules-reference", "verify_ref.py")
    if not os.path.exists(p):
        return None
    spec = importlib.util.spec_from_file_location("verify_ref", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("case, paths", [("repeat", "16"), ("mixed", "16"), ("repeat", "64")])
def test_the_rules_give_the_lengths_of_the_prototype(data, cli, tmp_path, case, paths, monkeypatch):
    """argjacent/verify.py holds the parts of the prototype's rule file that its switches use, and since 0.4 one term more
    (the rise of the depth, in the flow rule) and the count from reads far beyond a place. Without these both must give the same length for every candidate on
    the records of the strict measure; and what argjacent reports is the longer of the lengths under the two measures."""
    ref = _reference()
    if ref is None:
        pytest.skip("the prototype's rule file is not in this folder")
    out = str(tmp_path / "k")
    import subprocess
    import sys
    r = subprocess.run([sys.executable, "-m", "argjacent", "-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, case + "_1.fq"),
                        "-2", os.path.join(data, case + "_2.fq"), "-o", out, "-t", "2", "--keep", "--max-paths", paths],
                       capture_output=True, text=True, env=dict(os.environ, AW_KEEP_RECORDS="1"))      # the records are kept for this test only
    assert r.returncode == 0, r.stdout + r.stderr
    opts = ref.frozen()
    col = ref.COLUMNS.index("V_r8")
    with monkeypatch.context() as mp:
        mp.setattr(verify, "ZR", 1e9)                      # no rise is ever beyond chance: the flow rule of 0.3.0
        mp.setattr(verify, "FARB", 10 ** 9)                # ... and no read counts as far beyond a place
        mp.setattr(verify, "BSIG", 0.0)                    # ... and the budget needs no margin over its charges
        mp.setattr(verify, "FCAP", False)                  # ... and the flow rule is not held to the budget
        mp.setattr(verify, "SOFTW", 0)                     # ... and reads that leave in one base are not charged
        mp.setattr(verify, "INSIDE", 100)                  # ... and joins are looked for in the gene's last 100 bases only
        mp.setattr(verify, "GINSIDE", False)               # ... and pairs "from the gene" are those that begin in the gene part
        mp.setattr(verify, "UNEV", 0.0)                    # ... and depth is taken to be even
        mp.setattr(verify, "AFTERG", False)                # ... and reads that leave in the gene's last 20 bases are among the charges
        mp.setattr(verify, "BRIDGE", 0)                    # ... and the certain part ends at the first position nothing vouches for
        mine = verify.verify_gene(out + ".tmp/v/0")
    theirs = ref.verify_gene(out + ".tmp/v/0", opts)
    assert len(mine) == len(theirs) > 1
    for m, t in zip(mine, theirs):
        assert m["name"].split("|")[:3] == list(t[:3])
        assert (m["flank_len"], m["verified"], m["joins"], m["open_joins"]) == (t[3], t[col], t[-3], t[-2]), m["name"]
        assert abs(m["depth"] - t[-1]) < 1e-9
    doc = json.load(open(out + ".flanks.json"))["genes"]["demo_gene"]
    certain = sorted(p["context_len"] for sd in ("left", "right") for p in doc[sd]["paths"])
    strict = verify.verify_gene(out + ".tmp/v/0")
    tolerant = verify.verify_gene(out + ".tmp/vt/0", patches=True)
    assert [m["name"] for m in strict] == [m["name"] for m in tolerant]
    assert certain == sorted(max(0, min(max(a["verified"], b["verified"]), a["flank_len"]) - verify.CUT) for a, b in zip(strict, tolerant))
    with monkeypatch.context() as mp:
        mp.setattr(verify, "BRIDGE", 0)
        unbridged = verify.verify_gene(out + ".tmp/v/0")
    assert all(a["verified"] <= m["verified"] for a, m in zip(unbridged, mine))     # the terms since 0.3.0 can only shorten ...
    assert all(a["verified"] <= b["verified"] for a, b in zip(unbridged, strict))   # ... and bridging only lengthen


@pytest.mark.parametrize("case", ["mixed", "shared"])
def test_the_terms_of_0_4_work_in_the_direction_they_should(data, cli, tmp_path, case, monkeypatch):
    """Under both measures and for every candidate: the allowance for uneven depth can only shorten what is certain, bridging
    can only lengthen it and less so when held to stretches where little was charged, and the end-anchored claim is never
    shorter than the claim."""
    out = str(tmp_path / "d")
    import subprocess
    import sys
    r = subprocess.run([sys.executable, "-m", "argjacent", "-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, case + "_1.fq"),
                        "-2", os.path.join(data, case + "_2.fq"), "-o", out, "-t", "2", "--keep"],
                       capture_output=True, text=True, env=dict(os.environ, AW_KEEP_RECORDS="1"))
    assert r.returncode == 0, r.stdout + r.stderr
    for sub, tol in (("v", False), ("vt", True)):
        def lengths(**kw):
            with monkeypatch.context() as mp:
                for k, v in kw.items():
                    mp.setattr(verify, k, v)
                return verify.verify_gene("%s.tmp/%s/0" % (out, sub), patches=tol)
        full = lengths()
        assert len(full) > 1
        for m, even, unbridged, anygap in zip(full, lengths(UNEV=0.0), lengths(BRIDGE=0), lengths(BRQ=1e9)):
            assert unbridged["verified"] <= m["verified"] <= even["verified"], (sub, m["name"])
            assert m["verified"] <= anygap["verified"], (sub, m["name"])
            assert m["verified_end"] >= m["verified"] and m["shared_end"] >= 0, (sub, m["name"])
            assert unbridged["bridged"] == 0
