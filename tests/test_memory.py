"""Memory. The map of the sample is built in parts when it does not fit in what the run may use (--max-mem); that may
not change a result. What the machine, a container or a batch job allows is found out."""
import json
import os
import subprocess
import sys

from argjacent import cli as c


def test_what_the_machine_has():
    mb = c.machine_mb()
    assert mb is None or mb > 100


def test_a_container_or_batch_job_that_allows_less_than_the_machine_has(tmp_path):
    def put(path, text):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        open(path, "w").write(text)
    v2 = str(tmp_path / "v2")
    put(v2 + "/proc/self/cgroup", "0::/system.slice/job_17/step_0\n")
    put(v2 + "/sys/fs/cgroup/system.slice/job_17/step_0/memory.max", "max\n")
    put(v2 + "/sys/fs/cgroup/system.slice/job_17/memory.max", "%d\n" % (512 << 20))
    put(v2 + "/sys/fs/cgroup/memory.max", "max\n")
    assert c.cgroup_limits(v2) == [512 << 20]
    assert c.machine_mb(v2) in (512.0, c.machine_mb(str(tmp_path / "nothing")))      # the smaller of the two
    assert c.machine_mb(v2) <= 512.0
    v1 = str(tmp_path / "v1")
    put(v1 + "/proc/self/cgroup", "5:cpu,cpuacct:/docker/abc\n4:memory:/docker/abc\n")
    put(v1 + "/sys/fs/cgroup/memory/docker/abc/memory.limit_in_bytes", "%d\n" % (256 << 20))
    put(v1 + "/sys/fs/cgroup/memory/memory.limit_in_bytes", "9223372036854771712\n")        # v1: no limit
    assert c.cgroup_limits(v1) == [256 << 20] and c.machine_mb(v1) <= 256.0
    assert c.cgroup_limits(str(tmp_path / "nothing")) == []


def _args(data, case, out, *more):
    return ["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, case + "_1.fq"), "-2", os.path.join(data, case + "_2.fq"),
            "-o", out, "-t", "2"] + list(more)


def test_the_memory_limit_changes_no_result(data, cli, tmp_path):
    """with room for everything, and with so little that the map of the sample is built in parts"""
    a, b, d = str(tmp_path / "a"), str(tmp_path / "b"), str(tmp_path / "d")
    cli(_args(data, "mixed", a))
    cli(_args(data, "mixed", b, "--max-mem", "1"))
    rc, out, err = cli(_args(data, "mixed", d, "--max-mem", "0.0034"))             # 3.5 MB: 2 for the Bloom filter, too little for the map in one piece
    assert "building it in" in out
    for ext in (".summary.tsv", ".contexts.tsv", ".flanks.fa"):
        ref = open(a + ext).read()
        assert open(b + ext).read() == ref and open(d + ext).read() == ref, ext
    ja = json.load(open(a + ".flanks.json"))
    assert json.load(open(d + ".flanks.json"))["genes"] == ja["genes"]
    assert int(json.load(open(d + ".flanks.json"))["aw_stats"]["map_parts"]) > 1 and ja["aw_stats"]["map_parts"] == "1"


def test_a_map_in_parts_forced_through_the_command(data, tmp_path):
    outs = []
    for parts in ("", "3"):
        out = str(tmp_path / ("p" + parts) / "x")
        env = dict(os.environ)
        if parts:
            env["AW_TEST_PARTS"] = parts
        r = subprocess.run([sys.executable, "-m", "argjacent"] + _args(data, "repeat", out), capture_output=True, text=True, env=env)
        assert r.returncode == 0, r.stderr
        outs.append(json.load(open(out + ".flanks.json")))
    assert outs[0]["genes"] == outs[1]["genes"] and outs[1]["aw_stats"]["map_parts"] == "3"


def test_the_memory_limit_must_be_above_zero(data, cli, tmp_path):
    for bad in ("0", "-1"):
        rc, out, err = cli(_args(data, "clean", str(tmp_path / "x"), "--max-mem=" + bad), expect=1)
        assert "--max-mem" in err


def test_too_little_memory_is_said_plainly(data, cli, tmp_path):
    rc, out, err = cli(_args(data, "clean", str(tmp_path / "x"), "--max-mem", "0.001"), expect=1)
    assert "too little" in out + err and "Traceback" not in err
    assert not os.path.exists(str(tmp_path / "x.tmp"))
