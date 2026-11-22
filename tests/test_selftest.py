"""The built-in self-test: six simulated cases with a known answer."""
import json
import os


def test_the_cases_match_the_known_answer(data, cli, tmp_path):
    from argjacent import selftest
    for c in selftest.CASES:
        cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, c + "_1.fq"), "-2", os.path.join(data, c + "_2.fq"),
             "-o", os.path.join(data, c), "-t", "2"])
    lines = []
    assert selftest.check(data, echo=lines.append), "\n".join(lines)


def test_gene_in_two_genomes_gives_two_certain_contexts(data, cli):
    out = os.path.join(data, "mixed2")
    cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "mixed_1.fq"), "-2", os.path.join(data, "mixed_2.fq"),
         "-o", out, "-t", "1"])
    res = json.load(open(out + ".flanks.json"))["genes"]["demo_gene"]
    truth = json.load(open(os.path.join(data, "truth.json")))["mixed"]
    for side in ("left", "right"):
        paths = res[side]["paths"]
        assert len(paths) == 2
        assert all(p["context_len"] >= 1000 for p in paths)
        assert {p["flank"] for p in paths} == {truth["A"][side], truth["B"][side]}


def test_the_self_test_command(tmp_path):
    import subprocess
    import sys
    r = subprocess.run([sys.executable, "-m", "argjacent.selftest", str(tmp_path / "st")], capture_output=True, text=True)
    assert r.returncode == 0 and "All checks passed." in r.stdout, r.stdout[-2000:] + r.stderr[-2000:]
