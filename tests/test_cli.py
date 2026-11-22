"""The command: outputs, errors, repeatability, saved maps."""
import json
import os

import pytest


def args(data, case, out, *more):
    return ["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, case + "_1.fq"), "-2", os.path.join(data, case + "_2.fq"),
            "-o", out, "-t", "1"] + list(more)


def test_version(cli):
    rc, out, _ = cli(["--version"])
    assert out.startswith("argjacent ")


def test_outputs_and_columns(data, cli, tmp_path):
    out = str(tmp_path / "r")
    cli(args(data, "repeat", out))
    for ext in (".summary.tsv", ".contexts.tsv", ".flanks.fa", ".flanks.json", ".log"):
        assert os.path.getsize(out + ext) > 0
    assert not os.path.exists(out + ".tmp")
    head = open(out + ".summary.tsv").readline().rstrip("\n").split("\t")
    assert head == ["gene", "side", "status", "contexts", "certain_bp", "longest_bp", "contexts_certain_1kb", "read_depth", "gene_end",
                    "certain_end_bp"]
    rows = [l.rstrip("\n").split("\t") for l in open(out + ".summary.tsv")][1:]
    status = {r[1]: r[2] for r in rows}
    assert status == {"left": "complete", "right": "partial"}
    assert all(len(r) == len(head) for r in rows)
    ctx = [l.rstrip("\n").split("\t") for l in open(out + ".contexts.tsv")][1:]
    assert all(0 <= int(r[4]) <= int(r[3]) for r in ctx)                              # certain <= length
    by_side = {sd: [r for r in ctx if r[1] == sd] for sd in ("left", "right")}
    for sd, rs in by_side.items():                                                    # the best context first, named c1, c2, ...
        assert [r[2] for r in rs] == ["c%d" % (i + 1) for i in range(len(rs))]
        assert [int(r[4]) for r in rs] == sorted((int(r[4]) for r in rs), reverse=True)
    assert by_side["left"][0][6:8] == ["reached_target", "end_of_candidate"]
    assert by_side["right"][0][7] == "another_locus_on_the_path"                      # certainty ends at the far end of the repeat
    n_fa = sum(1 for l in open(out + ".flanks.fa") if l.startswith(">"))
    assert n_fa == len(ctx)
    head = open(out + ".contexts.tsv").readline().rstrip("\n").split("\t")
    assert head == ["gene", "side", "context", "length_bp", "certain_bp", "read_depth", "stop", "certain_stop", "other_genes", "gene_end", "main",
                    "certain_end_bp", "shared_end_bp"]
    assert all(int(r[4]) <= int(r[11]) <= int(r[3]) and int(r[12]) >= 0 for r in ctx)        # the end claim is never shorter than the claim
    for sd in ("left", "right"):                                                   # one main context per gene side
        assert [r[10] for r in ctx if r[1] == sd].count("yes") == 1
    assert sum(1 for l in open(out + ".flanks.fa") if l.startswith(">") and l.rstrip().endswith("main=yes")) == 2


def test_max_paths(data, cli, tmp_path):
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    cli(args(data, "repeat", a))
    cli(args(data, "repeat", b, "--max-paths", "2"))
    na = sum(1 for l in open(a + ".contexts.tsv")) - 1
    nb = sum(1 for l in open(b + ".contexts.tsv")) - 1
    assert nb <= na and nb <= 2 + 1                                                # 2 on the right side at most, 1 on the left
    rc, _, err = cli(args(data, "repeat", b, "--max-paths", "0"), expect=None)
    assert rc != 0 and "--max-paths" in err


def test_same_input_same_output(data, cli, tmp_path):
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    cli(args(data, "mixed", a))
    cli(args(data, "mixed", b, "-t", "2"))
    ja, jb = json.load(open(a + ".flanks.json")), json.load(open(b + ".flanks.json"))
    assert ja["genes"] == jb["genes"]


def test_saved_map_gives_the_same_result(data, cli, tmp_path):
    a, b, m = str(tmp_path / "a"), str(tmp_path / "b"), str(tmp_path / "sample.awmap")
    cli(args(data, "repeat", a, "--save-map", m))
    assert os.path.getsize(m) > 0
    cli(args(data, "repeat", b, "--map", m))
    assert json.load(open(a + ".flanks.json"))["genes"] == json.load(open(b + ".flanks.json"))["genes"]


def test_gene_not_in_the_reads(data, cli, tmp_path):
    import random
    R = random.Random(7)
    g = str(tmp_path / "absent.fa")
    open(g, "w").write(">absent\n" + "".join(R.choice("ACGT") for _ in range(900)) + "\n")
    out = str(tmp_path / "x")
    cli(["-g", g, "-1", os.path.join(data, "clean_1.fq"), "-2", os.path.join(data, "clean_2.fq"), "-o", out])
    rows = [l.rstrip("\n").split("\t") for l in open(out + ".summary.tsv")][1:]
    assert [r[2] for r in rows] == ["not_found", "not_found"]


def test_follow_everything_still_works(data, cli, tmp_path):
    out = str(tmp_path / "all")
    cli(args(data, "clean", out, "--min-share", "0"))
    res = json.load(open(out + ".flanks.json"))["genes"]["demo_gene"]
    assert all(res[s]["paths"][0]["context_len"] >= 9900 for s in ("left", "right"))


def test_reads_out(data, cli, tmp_path):
    out = str(tmp_path / "ro")
    cli(args(data, "clean", out, "--reads-out"))
    n1 = sum(1 for _ in open(out + ".near_1.fq"))
    assert n1 > 0 and n1 % 4 == 0 and n1 == sum(1 for _ in open(out + ".near_2.fq"))


@pytest.mark.parametrize("bad, text", [
    (["-g", "nope.fa", "-1", "a", "-2", "b", "-o", "x"], "file not found"),
])
def test_missing_file(cli, bad, text):
    rc, _, err = cli(bad, expect=1)
    assert text in err


def test_same_file_twice(data, cli, tmp_path):
    f = os.path.join(data, "clean_1.fq")
    rc, _, err = cli(["-g", os.path.join(data, "gene.fa"), "-1", f, "-2", f, "-o", str(tmp_path / "x")], expect=1)
    assert "same file" in err


def test_duplicate_gene_names(data, cli, tmp_path):
    g = str(tmp_path / "dup.fa")
    open(g, "w").write(">a\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n>a\nACGT\n")
    rc, _, err = cli(["-g", g, "-1", os.path.join(data, "clean_1.fq"), "-2", os.path.join(data, "clean_2.fq"),
                      "-o", str(tmp_path / "x")], expect=1)
    assert "occurs twice" in err


def test_gzip_input_gives_the_same_result(data, cli, tmp_path):
    import gzip
    import shutil
    for m in ("1", "2"):
        with open(os.path.join(data, "clean_%s.fq" % m), "rb") as fi, gzip.open(str(tmp_path / ("c_%s.fq.gz" % m)), "wb") as fo:
            shutil.copyfileobj(fi, fo)
    a, b = str(tmp_path / "plain"), str(tmp_path / "gz")
    cli(args(data, "clean", a, "-b", "24"))
    cli(["-g", os.path.join(data, "gene.fa"), "-1", str(tmp_path / "c_1.fq.gz"), "-2", str(tmp_path / "c_2.fq.gz"), "-o", b,
         "-t", "1", "-b", "24"])
    assert json.load(open(a + ".flanks.json"))["genes"] == json.load(open(b + ".flanks.json"))["genes"]


def _run_env(data, tmp_path, env, name="o"):
    import subprocess
    import sys
    e = dict(os.environ, **env)
    out = str(tmp_path / name / "x")
    r = subprocess.run([sys.executable, "-m", "argjacent", "-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "clean_1.fq"),
                        "-2", os.path.join(data, "clean_2.fq"), "-o", out, "-t", "2"], capture_output=True, text=True, env=e)
    return r, out


def test_a_failing_gene_is_reported_not_fatal(data, tmp_path):
    r, out = _run_env(data, tmp_path, {"AW_TEST_FAIL": "raise"})
    assert r.returncode == 1
    assert "failed" in r.stdout
    rows = [l.split("\t") for l in open(out + ".summary.tsv")][1:]
    assert rows and all(row[2] == "failed" for row in rows)
    assert not os.path.exists(out + ".tmp")


def test_a_dying_worker_gives_a_message_and_no_hang(data, tmp_path):
    r, out = _run_env(data, tmp_path, {"AW_TEST_FAIL": "die"})
    assert r.returncode == 1
    assert "worker process died" in r.stderr and "Traceback" not in r.stderr
    assert not os.path.exists(out + ".tmp")


def test_output_folder_that_cannot_be_written(data, tmp_path, cli):
    ro = tmp_path / "ro"
    ro.mkdir()
    os.chmod(ro, 0o555)
    try:
        if os.access(str(ro), os.W_OK):          # running as root: nothing to test
            return
        rc, out, err = cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "clean_1.fq"),
                            "-2", os.path.join(data, "clean_2.fq"), "-o", str(ro / "x")], expect=1)
        assert "cannot read or write" in err and "Traceback" not in err
    finally:
        os.chmod(ro, 0o755)


def test_same_output_whatever_the_hash_seed(data, tmp_path):
    """Python orders sets differently in every process (PYTHONHASHSEED); the result must not depend on it"""
    outs = []
    for seed in ("1", "2", "3"):
        r, out = _run_env(data, tmp_path, {"PYTHONHASHSEED": seed}, name="seed" + seed)
        assert r.returncode == 0, r.stderr
        outs.append(json.load(open(out + ".flanks.json"))["genes"])
    assert outs[0] == outs[1] == outs[2]


def test_many_genes_do_not_run_out_of_open_files(data, tmp_path):
    """300 genes with the limit of open files at 64 (macOS allows 256 by default)"""
    import resource
    import subprocess
    import sys
    gene = "".join(l.strip() for l in open(os.path.join(data, "gene.fa")) if not l.startswith(">"))
    fa = str(tmp_path / "many.fa")
    with open(fa, "w") as f:
        for i in range(300):
            f.write(">g%03d\n%s\n" % (i, gene[i % 50:500 + i % 50]))

    def low_limit():
        resource.setrlimit(resource.RLIMIT_NOFILE, (64, resource.getrlimit(resource.RLIMIT_NOFILE)[1]))
    out = str(tmp_path / "m" / "x")
    r = subprocess.run([sys.executable, "-m", "argjacent", "-g", fa, "-1", os.path.join(data, "clean_1.fq"),
                        "-2", os.path.join(data, "clean_2.fq"), "-o", out, "-t", "2", "-D", "300"],
                       capture_output=True, text=True, preexec_fn=low_limit)
    assert r.returncode == 0, r.stderr[-500:]
    rows = [l.split("\t") for l in open(out + ".summary.tsv")][1:]
    assert len(rows) == 600 and all(row[2] != "failed" for row in rows)


def test_blanks_in_the_gene_file_change_nothing(data, cli, tmp_path):
    """blanks and tabs at line ends, Windows line ends: aw and the Python part must see the same gene"""
    lines = open(os.path.join(data, "gene.fa")).read().split("\n")
    seq = "".join(l for l in lines if l and not l.startswith(">"))
    fa = str(tmp_path / "blank.fa")
    with open(fa, "w", newline="") as f:
        f.write(">demo_gene some words\r\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i + 60] + " \t\r\n")
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    cli(args(data, "clean", a))
    cli(["-g", fa] + args(data, "clean", b)[2:])
    assert json.load(open(a + ".flanks.json"))["genes"] == json.load(open(b + ".flanks.json"))["genes"]


def test_gene_names_that_are_not_plain_ascii_are_refused(data, cli, tmp_path):
    fa = str(tmp_path / "odd.fa")
    seq = "".join(l.strip() for l in open(os.path.join(data, "gene.fa")) if not l.startswith(">"))
    with open(fa, "wb") as f:
        f.write(b">demo_gene\xc2\xa0variant\n" + seq.encode() + b"\n")
    rc, out, err = cli(["-g", fa] + args(data, "clean", str(tmp_path / "o"))[2:], expect=1)
    assert "plain ASCII" in err and "Traceback" not in err


def test_a_folder_named_like_the_temporary_one_is_not_deleted(data, cli, tmp_path):
    out = str(tmp_path / "s")
    os.makedirs(out + ".tmp/my_notes")
    rc, o, err = cli(args(data, "clean", out), expect=1)
    assert "was not made by argjacent" in err
    assert os.path.isdir(out + ".tmp/my_notes")


def test_old_results_do_not_survive_a_failed_rerun(data, cli, tmp_path):
    out = str(tmp_path / "r")
    cli(args(data, "clean", out))
    assert os.path.exists(out + ".summary.tsv")
    rc, o, err = cli(args(data, "clean", out, "--map", os.path.join(data, "gene.fa")), expect=1)     # not a map: aw stops
    assert not os.path.exists(out + ".summary.tsv") and not os.path.exists(out + ".tmp")


def test_junk_and_adapter_reads_do_not_spoil_the_flanks(data, cli, tmp_path):
    """Real reads carry runs of one base and adapter read-through. With both added to a clean sample the run must end,
    and what it calls certain must still be the true flank."""
    import random
    R = random.Random(5)
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"
    out1, out2 = str(tmp_path / "j_1.fq"), str(tmp_path / "j_2.fq")
    for src, dst in ((os.path.join(data, "clean_1.fq"), out1), (os.path.join(data, "clean_2.fq"), out2)):
        lines = open(src).read().split("\n")
        with open(dst, "w") as f:
            for i in range(0, len(lines) - 3, 4):
                seq = lines[i + 1]
                if R.random() < 0.15:                       # the fragment was shorter than the read: adapter, then poly-A
                    cut = R.randrange(40, 110)
                    seq = (seq[:cut] + adapter + "A" * 150)[:len(seq)]
                f.write("%s\n%s\n+\n%s\n" % (lines[i], seq, lines[i + 3]))
            for j in range(400):                            # junk pairs: runs of one base, and a 2-base repeat
                s = R.choice(["G" * 150, "A" * 150, "AC" * 75])
                f.write("@junk%d\n%s\n+\n%s\n" % (j, s, "I" * 150))
    out = str(tmp_path / "junk")
    cli(["-g", os.path.join(data, "gene.fa"), "-1", out1, "-2", out2, "-o", out, "-t", "2"])
    truth = json.load(open(os.path.join(data, "truth.json")))["clean"]["A"]
    res = json.load(open(out + ".flanks.json"))["genes"]["demo_gene"]
    for side in ("left", "right"):
        for p in res[side]["paths"]:
            e = max(p["context_len"], p["end_len"])             # the end-anchored claim must be the true flank as well
            if side == "right":
                assert p["flank"][:e] == truth[side][:e]
            else:
                assert p["flank"][len(p["flank"]) - e:] == truth[side][len(truth[side]) - e:]
        best = res[side]["paths"][0]
        # (0.3.0: 1,000 or more on both sides. Since 0.4.0 the budget allows for uneven depth: adapter reads that come on two
        # at a place count as other loci with an upper limit of about ten reads each, and at this depth, about 30x, the
        # tolerant budget cannot stand clear of them for long: 797 on the left. The end-anchored claim still reaches 2.7 kb.)
        assert best["context_len"] >= 700, (side, best["context_len"], best["certain_stop"])
        assert best["end_len"] >= 1000, (side, best["end_len"])
        assert adapter[:25] not in best["flank"]


def test_the_json_says_what_the_tables_say(data, cli, tmp_path):
    out = str(tmp_path / "j")
    cli(args(data, "mixed", out))
    doc = json.load(open(out + ".flanks.json"))
    assert doc["version"] and doc["flank_bp"] == 10000 and doc["parameters"]["read_len"] == 150
    assert doc["parameters"]["certain_cut_bp"] == 30
    ctx = [l.rstrip("\n").split("\t") for l in open(out + ".contexts.tsv")][1:]
    fa = {}
    for l in open(out + ".flanks.fa"):
        if l.startswith(">"):
            name = l[1:].split()[0]
            fa[name] = ""
        else:
            fa[name] += l.strip()
    n = 0
    for sd in ("left", "right"):
        side = doc["genes"]["demo_gene"][sd]
        assert side["status"] in ("complete", "partial", "ended")
        for p in side["paths"]:
            row = [r for r in ctx if r[1] == sd and r[2] == p["id"]][0]
            assert int(row[3]) == p["flank_len"] == len(p["flank"]) and int(row[4]) == p["context_len"]
            assert fa["demo_gene|%s|%s" % (sd, p["id"])] == p["flank"]
            assert p["gene_end"] in ("sample", "database") and isinstance(p["main"], bool)
            n += 1
    assert n == len(ctx) == len(fa)


def test_other_genes_on_a_context_are_listed(data, cli, tmp_path):
    """a second gene 2 kb downstream of the first: each lists the other on its context, with direction and place"""
    truth = json.load(open(os.path.join(data, "truth.json")))["clean"]["A"]
    gene = "".join(l.strip() for l in open(os.path.join(data, "gene.fa")) if not l.startswith(">"))
    second = truth["right"][2000:2900]
    comp = str.maketrans("ACGT", "TGCA")
    fa = str(tmp_path / "two.fa")
    open(fa, "w").write(">first\n%s\n>second_rc\n%s\n" % (gene, second.translate(comp)[::-1]))
    out = str(tmp_path / "two")
    cli(["-g", fa, "-1", os.path.join(data, "clean_1.fq"), "-2", os.path.join(data, "clean_2.fq"), "-o", out, "-t", "2"])
    ctx = {(r[0], r[1]): r for r in (l.rstrip("\n").split("\t") for l in open(out + ".contexts.tsv")) if r[2] == "c1"}
    og = ctx[("first", "right")][8]
    assert og.startswith("second_rc[-]@") and "uncertain" not in og
    lo, hi = (int(x) for x in og.split("@")[1].split("(")[0].split("-"))
    assert 1950 <= lo <= 2100 and 2750 <= hi <= 2900
    assert ctx[("second_rc", "right")][8].startswith("first[-]@")          # seen from the second gene, the first lies downstream, reversed
    assert ctx[("first", "left")][8] == "-"


def test_options_of_the_old_caller_are_accepted_and_do_nothing(data, cli, tmp_path):
    a, b = str(tmp_path / "a"), str(tmp_path / "b")
    cli(args(data, "clean", a))
    rc, out, err = cli(args(data, "clean", b, "--strict", "--min-link", "2"))
    assert "no effect" in err
    assert json.load(open(a + ".flanks.json"))["genes"] == json.load(open(b + ".flanks.json"))["genes"]


def test_short_reads_are_told_apart(data, cli, tmp_path):
    """reads cut to 60 bp: too short for the rules; candidates are still given, nothing is called certain, and the run says so"""
    r1, r2 = str(tmp_path / "s_1.fq"), str(tmp_path / "s_2.fq")
    for src, dst in ((os.path.join(data, "clean_1.fq"), r1), (os.path.join(data, "clean_2.fq"), r2)):
        lines = open(src).read().split("\n")
        with open(dst, "w") as f:
            for i in range(0, len(lines) - 3, 4):
                f.write("%s\n%s\n+\n%s\n" % (lines[i], lines[i + 1][:60], lines[i + 3][:60]))
    out = str(tmp_path / "short")
    rc, o, err = cli(["-g", os.path.join(data, "gene.fa"), "-1", r1, "-2", r2, "-o", out, "-t", "2"])
    assert "too short" in o
    rows = [l.rstrip("\n").split("\t") for l in open(out + ".contexts.tsv")][1:]
    assert rows and all(r[4] == "0" and r[7] == "reads_too_short" for r in rows)


def test_keep_keeps_the_reads_and_candidates_not_the_records(data, cli, tmp_path):
    out = str(tmp_path / "k")
    cli(args(data, "repeat", out, "--keep"))
    kept = set(os.listdir(out + ".tmp"))
    assert {"aw.reads.tsv", "aw.genes.tsv", "aw.stats.txt", "seeds.fa", "c", "v"} <= kept
    assert os.path.getsize(os.path.join(out + ".tmp", "c", "0.cands.fa")) > 0
    assert [f for f in os.listdir(os.path.join(out + ".tmp", "v")) if f != "jobs.tsv"] == []


def test_hub_anchors_and_duplicate_gene_ends_give_identical_results(data, cli, tmp_path):
    """High-degree hub anchors (>16 links from chimeric poly-C + adapter tails) and multiple alleles sharing 400-bp
    gene ends give byte-identical maps and flanks across 1 and 4 threads."""
    import random
    R = random.Random(17)
    gene = "".join(l.strip() for l in open(os.path.join(data, "gene.fa")) if not l.startswith(">"))
    allele = gene[:400] + "".join(R.choice("ACGT") for _ in range(200)) + gene[-400:]
    gfa = str(tmp_path / "alleles.fa")
    open(gfa, "w").write(">g1\n%s\n>g2\n%s\n" % (gene, allele))
    hub = "CCCCCCCCCCCCCCCCCCCCCCTTTACAACCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    r1, r2 = str(tmp_path / "h_1.fq"), str(tmp_path / "h_2.fq")
    for src, dst in ((os.path.join(data, "strains_1.fq"), r1), (os.path.join(data, "strains_2.fq"), r2)):
        with open(src) as fi, open(dst, "w") as fo:
            fo.write(fi.read())
            for j in range(500):
                s = "".join(R.choice("ACGT") for _ in range(85)) + hub
                fo.write("@hub%d\n%s\n+\n%s\n" % (j, s[:150], "I" * 150))
    a, b = str(tmp_path / "t1"), str(tmp_path / "t4")
    m1, m4 = str(tmp_path / "m1.map"), str(tmp_path / "m4.map")
    cli(["-g", gfa, "-1", r1, "-2", r2, "-o", a, "-t", "1", "--save-map", m1])
    cli(["-g", gfa, "-1", r1, "-2", r2, "-o", b, "-t", "4", "--save-map", m4])
    assert open(m1, "rb").read() == open(m4, "rb").read()
    for ext in (".summary.tsv", ".contexts.tsv", ".flanks.fa"):
        assert open(a + ext, "rb").read() == open(b + ext, "rb").read()
    ja = json.load(open(a + ".flanks.json"))["genes"]
    assert ja == json.load(open(b + ".flanks.json"))["genes"]
    assert max(p["context_len"] for p in ja["g1"]["right"]["paths"]) >= 9900
    assert max(p["context_len"] for p in ja["g2"]["right"]["paths"]) >= 9900
