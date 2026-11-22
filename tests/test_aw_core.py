"""The C program aw on its own: the two read passes, the map of the sample in memory, in parts and in a file."""
import os
import struct
import subprocess

from argjacent import cli


def test_found_and_reports_version():
    for name in ("aw", "awx", "awv"):
        p = cli.find_bin(name)
        assert p, "the C program %s is not built (source folder: sh install.sh)" % name
        assert subprocess.run([p, "-v"], capture_output=True, text=True).stdout.startswith(name + " ")


def test_anchors_are_the_same_on_both_strands(tmp_path):
    """a sequence and its reverse complement give the same anchors (canonical k-mers)"""
    import random
    R = random.Random(3)
    s = "".join(R.choice("ACGT") for _ in range(600))
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    fa = str(tmp_path / "s.fa")
    open(fa, "w").write(">fwd\n%s\n>rev\n%s\n" % (s, rc))
    out = subprocess.run([cli.find_aw(), "-A", fa], capture_output=True, text=True, check=True).stdout.split("\n")
    fwd = {l.split("\t")[2] for l in out if l.startswith("fwd")}
    rev = {l.split("\t")[2] for l in out if l.startswith("rev")}
    assert fwd and fwd == rev
    assert 600 / 25 < len(fwd) < 600 / 5                    # about one anchor per 10 bases at w=20


def test_usage_error():
    assert subprocess.run([cli.find_aw()], capture_output=True).returncode == 1
    assert subprocess.run([cli.find_aw(), "-Z"], capture_output=True).returncode == 1


def _aw(args, cwd):
    return subprocess.run([cli.find_aw()] + args, capture_output=True, text=True, cwd=cwd)


OUT = ["o.genes.tsv", "o.reads.tsv", "o_1.fq", "o_2.fq", "o.map"]
RESULT = ["o.genes.tsv", "o.reads.tsv", "o_1.fq", "o_2.fq"]             # what does not depend on how the map was built


def _run(data, tmp_path, name, threads, r1=None, r2=None, expect=0, more=(), save=True):
    d = tmp_path / name
    d.mkdir()
    r = _aw(["-t", str(threads), "-f"] + (["-S", "o.map"] if save else []) + list(more) +
            ["-g", os.path.join(data, "gene.fa"), "-o", "o", r1 or os.path.join(data, "mixed_1.fq"), r2 or os.path.join(data, "mixed_2.fq")], str(d))
    assert r.returncode == expect, r.stderr
    return d, r.stderr


def test_any_number_of_threads_gives_the_same_files(data, tmp_path):
    """the map, the walk and the kept reads are the same, byte for byte, with 1, 2, 3 and 9 threads"""
    ref, _ = _run(data, tmp_path, "t1", 1)
    for f in OUT:
        assert os.path.getsize(ref / f) > 0, f
    for t in (2, 3, 9):
        d, _ = _run(data, tmp_path, "t%d" % t, t)
        for f in OUT:
            assert (d / f).read_bytes() == (ref / f).read_bytes(), "%s differs with %d threads" % (f, t)


def test_the_kept_reads_are_pairs_with_their_qualities(data, tmp_path):
    d, _ = _run(data, tmp_path, "k", 2)
    rows = [l.rstrip("\n").split("\t") for l in open(d / "o.reads.tsv")]
    assert rows and len(rows) % 2 == 0
    for a, b in zip(rows[0::2], rows[1::2]):
        assert a[0] == b[0] and (a[1], b[1]) == ("1", "2")
        assert len(a) == 5 and a[3] == "" and len(a[2]) == len(a[4]) == 150
    genes = [l.rstrip("\n").split("\t") for l in open(d / "o.genes.tsv")]
    assert genes[0][:4] == ["gene", "length", "anchors", "found"] and genes[1][0] == "demo_gene" and int(genes[1][3]) > 50


def test_a_map_built_in_parts_gives_the_same_reads(data, tmp_path):
    """the map built in one piece in memory, saved to a file, and built in 2, 3 and 7 parts: the walk reaches the same anchors"""
    mem, _ = _run(data, tmp_path, "mem", 2, save=False)
    assert not (mem / "o.map").exists()
    for parts in (1, 2, 3, 7):
        d, err = _run(data, tmp_path, "p%d" % parts, 3, more=["-P", str(parts)], save=parts == 1)
        for f in RESULT:
            assert (d / f).read_bytes() == (mem / f).read_bytes(), "%s differs with the map in %d parts" % (f, parts)
        assert ("map_parts=%d" % parts) in err
        assert parts == 1 or not (d / "o.map").exists()        # the temporary map file is gone
        walk = [l for l in err.split("\n") if l.startswith("walk ")]
        assert walk == [l for l in _run(data, tmp_path, "w%d" % parts, 1, save=False)[1].split("\n") if l.startswith("walk ")]


def test_a_memory_limit_makes_parts_and_changes_nothing(data, tmp_path):
    """-m: when the map does not fit, it is built in parts by itself; too little for even the Bloom filter is refused"""
    ref, _ = _run(data, tmp_path, "ref", 2, more=["-b", "24"], save=False)
    d, err = _run(data, tmp_path, "tight", 2, more=["-b", "24", "-m", "3.2"], save=False)
    assert "building it in" in err and "map_parts=1 " not in err
    for f in RESULT:
        assert (d / f).read_bytes() == (ref / f).read_bytes(), f
    roomy, err = _run(data, tmp_path, "roomy", 2, more=["-b", "24", "-m", "500"], save=False)
    assert "map_parts=1 " in err and (roomy / "o.reads.tsv").read_bytes() == (ref / "o.reads.tsv").read_bytes()
    _, err = _run(data, tmp_path, "none", 2, more=["-b", "24", "-m", "2"], save=False, expect=2)
    assert "too little" in err
    _, err = _run(data, tmp_path, "forced", 2, more=["-b", "24", "-m", "3.2", "-P", "1"], save=False, expect=2)
    assert "does not fit" in err


def test_gzip_and_plain_input_give_the_same_files(data, tmp_path):
    import gzip
    import shutil
    gz = []
    for m in "12":
        gz.append(str(tmp_path / ("r%s.fq.gz" % m)))
        with open(os.path.join(data, "mixed_%s.fq" % m), "rb") as i, gzip.open(gz[-1], "wb", compresslevel=1) as o:
            shutil.copyfileobj(i, o)
    ref, _ = _run(data, tmp_path, "plain", 1)
    for t in (1, 4):
        d, _ = _run(data, tmp_path, "gz%d" % t, t, gz[0], gz[1])
        for f in OUT:
            assert (d / f).read_bytes() == (ref / f).read_bytes(), f


def test_cut_off_gzip_file_is_reported_and_survived(data, tmp_path):
    """a download that stopped half way: warn, use the pairs before that point, the same with and without threads"""
    import gzip
    import shutil
    gz = []
    for m in "12":
        gz.append(str(tmp_path / ("r%s.fq.gz" % m)))
        with open(os.path.join(data, "mixed_%s.fq" % m), "rb") as i, gzip.open(gz[-1], "wb", compresslevel=1) as o:
            shutil.copyfileobj(i, o)
    size = os.path.getsize(gz[1])
    with open(gz[1], "r+b") as f:
        f.truncate(size * 2 // 3)
    a, err1 = _run(data, tmp_path, "cut1", 1, gz[0], gz[1])
    b, err4 = _run(data, tmp_path, "cut4", 4, gz[0], gz[1])
    c, err3 = _run(data, tmp_path, "cut3", 3, gz[0], gz[1], more=["-P", "3"], save=False)
    assert "truncated" in err1 and "truncated" in err4 and err3.count("truncated") == 1     # said once, not once per part
    for f in OUT:
        assert (a / f).read_bytes() == (b / f).read_bytes(), f
    assert (a / "o.reads.tsv").read_bytes() == (c / "o.reads.tsv").read_bytes()
    assert 0 < os.path.getsize(a / "o.reads.tsv")


def test_files_of_different_length(data, tmp_path):
    """mate files with different numbers of reads: the pairs present in both are used, with a warning"""
    lines = open(os.path.join(data, "mixed_2.fq")).readlines()
    short = str(tmp_path / "short_2.fq")
    open(short, "w").writelines(lines[:4 * 3000])
    a, err1 = _run(data, tmp_path, "s1", 1, None, short)
    b, err4 = _run(data, tmp_path, "s4", 4, None, short)
    assert "aw:" in err1 and "aw:" in err4
    for f in OUT:
        assert (a / f).read_bytes() == (b / f).read_bytes(), f


def test_empty_reads(data, tmp_path):
    e1, e2 = str(tmp_path / "e_1.fq"), str(tmp_path / "e_2.fq")
    open(e1, "w").close()
    open(e2, "w").close()
    for t in (1, 4):
        d, _ = _run(data, tmp_path, "e%d" % t, t, e1, e2)
        assert os.path.getsize(d / "o.reads.tsv") == 0
    d, _ = _run(data, tmp_path, "ep", 2, e1, e2, more=["-P", "3"], save=False)
    assert os.path.getsize(d / "o.reads.tsv") == 0


def test_a_saved_map_is_checked_before_use(data, tmp_path):
    """a map of other reads, a damaged map, a map of another version and a cut-off map are refused, not used"""
    d, _ = _run(data, tmp_path, "m", 1)
    gene = os.path.join(data, "gene.fa")
    reads = [os.path.join(data, "mixed_1.fq"), os.path.join(data, "mixed_2.fq")]
    ok = _aw(["-G", "o.map", "-g", gene, "-o", "again"] + reads, str(d))
    assert ok.returncode == 0, ok.stderr
    assert (d / "again.reads.tsv").read_bytes() == (d / "o.reads.tsv").read_bytes()
    other = _aw(["-G", "o.map", "-g", gene, "-o", "other", os.path.join(data, "clean_1.fq"), os.path.join(data, "clean_2.fq")], str(d))
    assert other.returncode == 1 and "not the same reads" in other.stderr
    wrong_k = _aw(["-k", "25", "-G", "o.map", "-g", gene, "-o", "k"] + reads, str(d))
    assert wrong_k.returncode == 1 and "was made with -k 31" in wrong_k.stderr
    raw = bytearray((d / "o.map").read_bytes())
    off_dir = struct.unpack_from("<Q", raw, 48)[0]
    bad = bytearray(raw)
    struct.pack_into("<Q", bad, off_dir, 1 << 40)                                # the first part claims far more anchors than the file holds
    (d / "bad.map").write_bytes(bytes(bad))
    r = _aw(["-G", "bad.map", "-g", gene, "-o", "bad"] + reads, str(d))
    assert r.returncode == 1 and "damaged" in r.stderr
    (d / "cut.map").write_bytes(bytes(raw[:len(raw) // 2]))
    r = _aw(["-G", "cut.map", "-g", gene, "-o", "cut"] + reads, str(d))
    assert r.returncode == 1 and "shorter" in r.stderr
    (d / "long.map").write_bytes(bytes(raw) + b"\0" * 64)
    r = _aw(["-G", "long.map", "-g", gene, "-o", "long"] + reads, str(d))
    assert r.returncode == 1 and "longer" in r.stderr
    old = bytearray(raw)
    old[:8] = b"AWMAP01\0"
    (d / "old.map").write_bytes(bytes(old))
    r = _aw(["-G", "old.map", "-g", gene, "-o", "old"] + reads, str(d))
    assert r.returncode == 1 and "another version" in r.stderr
    r = _aw(["-G", gene, "-g", gene, "-o", "nomap"] + reads, str(d))
    assert r.returncode == 1 and "not a map" in r.stderr


def test_windows_line_ends_in_the_reads(data, tmp_path):
    """CRLF reads give the same anchors and the same kept reads (without the CR) as LF reads"""
    crlf = []
    for m in "12":
        p = str(tmp_path / ("crlf_%s.fq" % m))
        with open(os.path.join(data, "mixed_%s.fq" % m), "rb") as i, open(p, "wb") as o:
            o.write(i.read().replace(b"\n", b"\r\n"))
        crlf.append(p)
    ref, _ = _run(data, tmp_path, "lf", 1)
    d, _ = _run(data, tmp_path, "crlf", 2, crlf[0], crlf[1])
    for f in ["o.genes.tsv", "o.reads.tsv", "o.map"]:
        assert (d / f).read_bytes() == (ref / f).read_bytes(), f
