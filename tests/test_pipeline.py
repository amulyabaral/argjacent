"""The steps between the kept reads and the rules: gene ends, the two runs of the proposer and how they are merged."""
import os
import random

from argjacent import pipeline

COMP = str.maketrans("ACGT", "TGCA")


def _seq(R, n):
    return "".join(R.choice("ACGT") for _ in range(n))


def _mutate(R, s, places):
    s = list(s)
    for i in places:
        s[i] = R.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def test_gene_ends_run_outward(tmp_path):
    R = random.Random(1)
    genes = [("long", _seq(R, 900)), ("absent", _seq(R, 500)), ("short", _seq(R, 250))]
    p = str(tmp_path / "seeds.fa")
    pipeline.write_seeds(genes, [0, 2], p)
    recs = open(p).read().split("\n")
    assert recs[0] == ">0|right|400" and recs[1] == genes[0][1][-400:]
    assert recs[2] == ">0|left|400" and recs[3] == genes[0][1][:400].translate(COMP)[::-1]
    assert recs[4] == ">2|right|250" and recs[5] == genes[2][1]              # a gene shorter than 400 bp is its own gene end
    assert not any(l.startswith(">1|") for l in recs)                        # a gene that was not found gets no gene end


def _write_run(prefix, redo, cands):
    with open("%s.r%d.fa" % (prefix, redo), "w") as fa, open("%s.r%d.tsv" % (prefix, redo), "w") as tsv:
        tsv.write("gene\tside\tcand\tflank\tstop\tborn_at\tpenalty\tvariants\tforks\tone_read\tjoin_first\tmass_end\n")
        for i, side, n, gl, seq, stop in cands:
            fa.write(">%d|%s|%d|%d|0|0\n%s\n" % (i, side, n, gl, seq))
            tsv.write("%d\t%s\t%d\t%d\t%s\t0\t0\t0\t0\t0\t-99999\t1.0000\n" % (i, side, n, len(seq) - gl, stop))


def test_candidates_from_the_samples_own_gene_end_are_taken_if_it_is_near_the_gene(tmp_path):
    R = random.Random(2)
    gene = _seq(R, 800)
    genes = [("g", gene)]
    end_r, end_l = gene[-400:], gene[:400].translate(COMP)[::-1]
    fl = [_seq(R, 500) for _ in range(6)]
    near = _mutate(R, end_r, range(110, 400, 30))            # 10 differences in the last 300 bases: the sample's version of the gene
    far = _mutate(R, end_r, range(105, 400, 7))              # 43 differences: a relative of the gene
    pre, cdir = str(tmp_path / "x"), str(tmp_path / "c")
    os.makedirs(cdir)
    _write_run(pre, 0, [(0, "right", 1, 400, end_r + fl[0], "reached target"), (0, "left", 1, 400, end_l + fl[1], "ran out of reads"),
                        (0, "left", 2, 400, end_l + fl[2], "loop")])
    _write_run(pre, 300, [(0, "right", 1, 400, far + fl[3], "reached target"), (0, "right", 2, 400, near + fl[4], "reached target"),
                          (0, "left", 1, 400, _mutate(R, end_l, range(101, 400, 6)) + fl[5], "reached target")])
    cands, st = pipeline.merge(genes, pre, cdir)
    right = [c for c in cands[0] if c["side"] == "right"]
    left = [c for c in cands[0] if c["side"] == "left"]
    # right: the run from the sample's gene end has one candidate near enough; the one that drifted to a relative is dropped
    assert [(c["cand"], c["gene_end"], c["gene_end_diffs"], c["run_cand"]) for c in right] == [(1, "sample", 10, 2)]
    assert right[0]["seq"] == near + fl[4]
    # left: every candidate of that run is too far from the gene, so the candidates that start from the gene as given are taken
    assert [(c["cand"], c["gene_end"], c["stop"]) for c in left] == [(1, "database", "ran out of reads"), (2, "database", "loop")]
    assert st["sides_sample_end"] == 1 and st["sides_database_end"] == 1 and st["sides_sample_end_too_far"] == 1
    assert st["dropped_too_far"] == 1 and st["candidates"] == 3
    names = [l[1:].strip() for l in open(os.path.join(cdir, "0.cands.fa")) if l.startswith(">")]
    assert names == ["0|right|1|400|0|0", "0|left|1|400|0|0", "0|left|2|400|0|0"]


def test_gene_names_never_reach_file_names_or_field_separators(data, cli, tmp_path):
    """a name with the separators of other tools in it (CARD's headers have '|'), a slash and a colon"""
    seq = "".join(l.strip() for l in open(os.path.join(data, "gene.fa")) if not l.startswith(">"))
    fa = str(tmp_path / "odd.fa")
    name = "gb|X1.1|+|0-1000|ARO:3000/x:y|demo"
    open(fa, "w").write(">%s\n%s\n" % (name, seq))
    out = str(tmp_path / "o")
    cli(["-g", fa, "-1", os.path.join(data, "clean_1.fq"), "-2", os.path.join(data, "clean_2.fq"), "-o", out, "-t", "2"])
    rows = [l.rstrip("\n").split("\t") for l in open(out + ".summary.tsv")][1:]
    assert [r[0] for r in rows] == [name, name] and [r[2] for r in rows] == ["complete", "complete"]
