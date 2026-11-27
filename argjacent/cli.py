#!/usr/bin/env python3
"""ARGjacent: certified genomic contexts of genes from paired short reads, without assembly.

Three C programs and a verification module do the work:
  aw    reads the reads twice: a thin map of the whole sample, a walk outward from every gene on that map, and the read
        pairs that lie near the genes
  awx   proposes contexts: from each gene end it extends base by base with the kept reads, following up to N arms
        where they disagree
  awv   lays the kept reads on every proposed context and records what they say at each position
  argjacent/verify.py   says from those records how far each context is certainly a context of the gene

and this file writes:
  PREFIX.summary.tsv    one line per gene side: status, number of contexts, how far the best one is certain
  PREFIX.contexts.tsv   one line per context (candidate flank): length, certain part, read depth, why it stops
  PREFIX.flanks.fa      every context (left = upstream of the gene as given, right = downstream)
  PREFIX.flanks.json    everything (including 100-bp depth profiles and certified two-sided pairs)
  PREFIX.pairs.tsv      certified two-sided [left + gene + right] flank pairs (unless --no-pair)
  PREFIX.gfa            lossless compacted gene-centered GFA v1.0 graph (unless --no-gfa)
  PREFIX.log            details of the run

Use `argjacent-cohort` (or `argjacent cohort`) across S >= 2 samples to resolve tied depths and post-repeat exits.
"""
import argparse
import concurrent.futures
import errno
import json
import math
import multiprocessing
import os
import re
import resource
import shutil
import signal
import subprocess
import sys
import time

from . import __version__ as VERSION

HERE = os.path.dirname(os.path.abspath(__file__))
PATHS_DEFAULT = 16       # candidate contexts followed per gene side and gene-end setting


def find_bin(name):
    """one of the C programs aw, awx, awv: $AW_BIN (for aw), $AWX_BIN, $AWV_BIN if set, else the one built into this
    package (argjacent/bin/)"""
    for p in (os.environ.get(name.upper() + "_BIN"), os.path.join(HERE, "bin", name)):
        if p and os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return None


def find_aw():
    return find_bin("aw")


def fail(msg):
    sys.stderr.write("argjacent: " + msg + "\n")
    sys.exit(1)


ASCII_SPACE = re.compile(r"[ \t\v\f\r\n]+")              # what C's sscanf("%s") stops at; Python's split() also stops at no-break spaces


def read_genes(path):
    """[(name, seq)] using the first word of each header, as aw does; checks names."""
    recs, name, seq = [], None, []
    with open(path, encoding="latin-1") as f:              # any byte is a character: no decoding error on odd headers
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(seq)))
                words = ASCII_SPACE.split(line[1:].strip(" \t\v\f\r\n"))
                if not words or not words[0]:
                    fail("%s: a FASTA header has no name" % path)
                name, seq = words[0], []
                if len(name) > 1000 or not all(33 <= ord(ch) <= 126 for ch in name):
                    fail("%s: the gene name '%s' has a character that is not plain ASCII (accents, no-break spaces and the "
                         "like). Names end up in tables and file names: please rename the gene" % (path, name[:60]))
            elif line.strip():
                if name is None:
                    fail("%s does not look like FASTA (no '>' header before the first sequence)" % path)
                seq.append("".join(line.split()))
    if name is not None:
        recs.append((name, "".join(seq)))
    if not recs:
        fail("%s has no sequences" % path)
    seen = set()
    for n, _ in recs:
        if n in seen:
            fail("gene name '%s' occurs twice in %s (names = first word of the header; they must be unique)" % (n, path))
        seen.add(n)
    return recs


def auto_bloom_bits(files, w):
    """Size the Bloom filter (which drops anchors seen only once) from the input size.
    Rough: gzip FASTQ ~4x smaller than plain; ~45% of FASTQ bytes are bases; ~2/(w+1) anchors per base;
    assume up to half of all anchor sightings are distinct; ~9.8 anchors per 64-bit word keeps errors ~5%."""
    raw = sum(os.path.getsize(f) * (4.0 if f.endswith(".gz") else 1.0) for f in files)
    anchors = raw * 0.45 * 2.0 / (w + 1)
    words = max(anchors * 0.5 / 9.8, 1.0)
    return min(max(6 + int(math.ceil(math.log(words, 2))), 24), 33)


def peak_mb(who=resource.RUSAGE_SELF):
    r = resource.getrusage(who).ru_maxrss
    return r / 1048576.0 if sys.platform == "darwin" else r / 1024.0


def parse_stats(path):
    st = {}
    with open(path) as f:
        for line in f:
            for tok in line.split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    st[k] = v
    return st


def cgroup_limits(root="/"):
    """memory limits, in bytes, of the control groups this process is in (Linux: a container, or a job of a batch
    system such as Slurm; they allow less than the machine has). Empty where there are none."""
    out = []
    try:
        with open(os.path.join(root, "proc/self/cgroup")) as f:
            lines = f.read().splitlines()
    except OSError:
        return out
    for line in lines:
        parts = line.split(":", 2)
        if len(parts) != 3:
            continue
        if parts[0] == "0" and parts[1] == "":                          # cgroup v2
            base, name = os.path.join(root, "sys/fs/cgroup"), "memory.max"
        elif "memory" in parts[1].split(","):                           # cgroup v1
            base, name = os.path.join(root, "sys/fs/cgroup/memory"), "memory.limit_in_bytes"
        else:
            continue
        path = parts[2].strip("/")
        while True:                                                     # the group and every group above it
            try:
                with open(os.path.join(base, path, name)) as f:
                    v = f.read().strip()
                if v.isdigit() and 0 < int(v) < (1 << 60):              # "max" and v1's huge number mean no limit
                    out.append(int(v))
            except OSError:
                pass
            if not path:
                break
            path = os.path.dirname(path)
    return out


def machine_mb(root="/"):
    """the memory this run can have, in MB: the machine's, or what a container or batch job is allowed if that is
    less. None if it cannot be told."""
    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        total = None
    for limit in cgroup_limits(root):
        if total is None or limit < total:
            total = limit
    return total / 1048576.0 if total else None


COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"reached target": "reached_target", "ran out of reads": "ran_out_of_reads", "work cap": "work_limit", "loop": "loop"}


def _merge_measures(strict_res, tol_res):
    by_name = {}
    for r in strict_res:
        r["measure"] = "strict"
        by_name[r["name"]] = r
    for r in tol_res:
        r["measure"] = "tolerant"
        old = by_name.get(r["name"])
        if old is None:
            by_name[r["name"]] = r
        else:
            be = max((old, r), key=lambda z: z["verified_end"])
            end = (be["verified_end"], be["stop_end"], max(old["shared_end"], r["shared_end"]))
            if r["verified"] > old["verified"]:
                r["depth"] = old["depth"]
                r["gene_depth"] = old["gene_depth"]
                r["depth_bins"] = old["depth_bins"]
                by_name[r["name"]] = r
            kept = by_name[r["name"]]
            kept["verified_end"], kept["stop_end"], kept["shared_end"] = end
    return by_name


def rules_both(job):
    """worker: strict and tolerant rules on one gene's records -> (gene number, {cand_name: result})"""
    number, v_pre, vt_pre, read_len, keep = job
    hook = os.environ.get("AW_TEST_FAIL")                  # tests only: "raise" or "die" in the worker
    if hook == "die":
        os._exit(9)
    if hook == "raise":
        raise RuntimeError("test failure")
    from . import verify
    s_res = verify.verify_gene(v_pre, read_len, patches=False)
    t_res = verify.verify_gene(vt_pre, read_len, patches=True)
    if not keep:
        for pre in (v_pre, vt_pre):
            for ext in (".idx", ".i32", ".ent", ".lev", ".tpl"):
                try:
                    os.remove(pre + ext)
                except OSError:
                    pass
    return number, _merge_measures(s_res, t_res)


def rules_one(job):
    """worker: the rules on one gene's records -> (gene number, [result per candidate])"""
    number, prefix, read_len, keep, patches = job
    hook = os.environ.get("AW_TEST_FAIL")                  # tests only: "raise" or "die" in the worker
    if hook == "die":
        os._exit(9)
    if hook == "raise":
        raise RuntimeError("test failure")
    from . import verify
    res = verify.verify_gene(prefix, read_len, patches=patches)
    if not keep:
        for ext in (".idx", ".i32", ".ent", ".lev", ".tpl"):
            try:
                os.remove(prefix + ext)
            except OSError:
                pass
    return number, res


def other_genes(aw, k, w, genes_fa, names, cands, tmp):
    """Other genes of the input that lie on a candidate: for each, where (flank bp), in which direction relative to the
    gene the candidate belongs to, and how many of its anchors are there. Needs 3 anchors and a fifth of the gene's
    anchors. -> {(gene number, side, cand): [..]}"""
    where, total = {}, {}
    out = subprocess.run([aw, "-k", str(k), "-w", str(w), "-A", genes_fa], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    if out.returncode != 0:
        return {}
    idx = {n: i for i, n in enumerate(names)}
    for line in out.stdout.splitlines():
        g, pos, h, st = line.split("\t")
        where.setdefault(h, []).append((idx[g], int(st)))
        total[idx[g]] = total.get(idx[g], 0) + 1
    allfa = os.path.join(tmp, "candidates.fa")
    with open(allfa, "w") as f:
        for i, cs in cands.items():
            for c in cs:
                f.write(">%d|%s|%d\n%s\n" % (i, c["side"], c["cand"], c["seq"]))
    out = subprocess.run([aw, "-k", str(k), "-w", str(w), "-A", allfa], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    os.remove(allfa)
    if out.returncode != 0:
        return {}
    part = {(i, c["side"], c["cand"]): c["gene_part"] for i, cs in cands.items() for c in cs}
    hits = {}
    for line in out.stdout.splitlines():
        name, pos, h, st = line.split("\t")
        if h not in where:
            continue
        i, side, n = name.split("|")
        key = (int(i), side, int(n))
        for g, gst in where[h]:
            if g != key[0]:
                same = (int(st) == gst) == (side == "right")       # a left candidate runs against the gene's direction
                hits.setdefault(key, {}).setdefault(g, []).append((max(0, int(pos) - part[key]), same))
    res = {}
    for key, per in hits.items():
        found = []
        for g, v in per.items():
            if len(v) >= 3 and len(v) >= 0.2 * total[g]:
                found.append({"gene": names[g], "from_bp": min(b for b, _ in v), "to_bp": max(b for b, _ in v), "anchors": len(v),
                              "of": total[g], "strand": "+" if 2 * sum(1 for _, s in v if s) >= len(v) else "-"})
        if found:
            res[key] = sorted(found, key=lambda x: x["from_bp"])
    return res


def status_of(paths, flank, cut):
    """complete: a context is certain over the whole flank length; partial: certain up to a point, candidates beyond it;
    ended: the best context ends before the flank length because the reads ran out, and is certain up to its end (or
    there is no context at all on this side)"""
    if not paths:
        return "ended", None
    best = max(paths, key=lambda p: (p["context_len"], p["flank_len"]))
    if best["context_len"] >= flank - cut:
        return "complete", best
    if best["flank_len"] < flank and best["context_len"] >= best["flank_len"] - cut:
        return "ended", best
    return "partial", best


def main(argv=None):
    raw_argv = sys.argv[1:] if argv is None else list(argv)
    if raw_argv and raw_argv[0] == "cohort":
        from . import cohort
        return cohort.main(raw_argv[1:])

    ap = argparse.ArgumentParser(
        prog="argjacent",
        description="Certified genomic contexts of genes from paired short reads of a metagenome or isolate "
                    "(anchor walk: two passes over the reads, no assembly; every context says how far it is certain).",
        epilog="Example: argjacent -g genes.fa -1 sample_R1.fq.gz -2 sample_R2.fq.gz -o results/sample "
               "(multi-sample cohort deconvolution: argjacent-cohort / argjacent cohort -h)")
    ap.add_argument("-g", "--genes", required=True,
                    help="FASTA of genes to extend, e.g. the genes of a resistance gene database that a read mapper found in the sample")
    ap.add_argument("-1", dest="r1", required=True, help="reads, mate 1 (FASTQ, plain or .gz)")
    ap.add_argument("-2", dest="r2", required=True, help="reads, mate 2 (FASTQ, plain or .gz)")
    ap.add_argument("-o", "--out", required=True, help="output prefix, e.g. results/sample1")
    ap.add_argument("-D", "--flank", type=int, default=10000, help="flank length per side in bp (default 10000)")
    ap.add_argument("-t", "--threads", type=int, default=min(8, os.cpu_count() or 1),
                    help="threads (default: up to 8). The results do not depend on it")
    ap.add_argument("--max-mem", type=float, default=None, metavar="GB",
                    help="memory the run may use (default: 60 %% of this machine's memory, or of what a container or batch "
                         "job is allowed). If the map of the sample does not fit, it is built in parts, with one more "
                         "pass over the reads per part and a temporary file. The results do not depend on it")
    ap.add_argument("--max-paths", type=int, default=PATHS_DEFAULT, metavar="N",
                    help="candidate contexts followed per gene side (default %d). More: rarer contexts are found, "
                         "the run takes longer" % PATHS_DEFAULT)
    ap.add_argument("--min-share", type=float, default=0.002,
                    help="the walk that chooses the reads follows a branch only if it carries at least this share of the "
                         "gene's reads still on their way at that distance (default 0.002; 0 = follow everything)")
    ap.add_argument("--save-map", metavar="FILE", help="save the map of the sample made in the first read pass")
    ap.add_argument("--map", metavar="FILE", help="use a saved map instead of the first read pass (same -k and -w)")
    ap.add_argument("-w", type=int, default=11,
                    help="anchor window (default 11). Larger windows (e.g. 20 or 40) use less memory in the read passes")
    ap.add_argument("-k", type=int, default=21, help="k-mer size of the anchors, 11-32 (default 21)")
    ap.add_argument("-b", type=int, default=None,
                    help="log2 of Bloom filter bits (default: chosen from input size; aw warns if too small)")
    ap.add_argument("--max-anchors", type=int, default=2000000, help="hard limit of anchors reached per gene side (default 2,000,000)")
    ap.add_argument("--reads-out", action="store_true",
                    help="also write the read pairs found near the genes (PREFIX.near_1.fq, PREFIX.near_2.fq), e.g. to assemble them")
    ap.add_argument("--no-gfa", action="store_true",
                    help="do not write PREFIX.gfa (compacted gene-centered GFA v1.0 graph)")
    ap.add_argument("--no-pair", action="store_true",
                    help="do not compute certified two-sided [left + gene + right] flank pairs (PREFIX.pairs.tsv)")
    ap.add_argument("--pair-bp", type=int, default=500, metavar="BP",
                    help="certified flank radius in bp used for two-sided flank pairing (default 500)")
    ap.add_argument("--keep", action="store_true",
                    help="keep intermediate files in PREFIX.tmp/ (the reads kept near the genes, the candidates, the statistics of the read pass)")
    ap.add_argument("--strict", action="store_true", help=argparse.SUPPRESS)                     # options of 0.2 that the new
    ap.add_argument("--min-link", default=None, choices=["auto", "1", "2"], help=argparse.SUPPRESS)    # method has no use for
    ap.add_argument("-v", "--version", action="version", version="argjacent " + VERSION)
    a = ap.parse_args(raw_argv)

    bins = {n: find_bin(n) for n in ("aw", "awx", "awv")}
    missing = [n for n, p in bins.items() if not p]
    if missing:
        fail("the program%s %s not found next to this package. In a source folder run install.sh first; "
             "or set AW_BIN, AWX_BIN and AWV_BIN to their paths." % ("s" if len(missing) > 1 else "", " and ".join("'%s'" % n for n in missing) + (" were" if len(missing) > 1 else " was")))
    for f in (a.genes, a.r1, a.r2) + ((a.map,) if a.map else ()):
        if not os.path.isfile(f):
            fail("file not found: " + f)
    if os.path.abspath(a.r1) == os.path.abspath(a.r2):
        fail("-1 and -2 are the same file; argjacent needs the two files of a read pair")
    if not 100 <= a.flank <= 1000000 or not 11 <= a.k <= 32 or a.w < 1 or not 0 <= a.min_share < 1 or a.threads < 1:
        fail("check -D (100 to 1,000,000), -k (11-32), -w (>= 1), --min-share (0 to below 1) and -t (>= 1)")
    if a.max_mem is not None and not a.max_mem > 0:
        fail("--max-mem must be a number of GB above 0")
    if not 1 <= a.max_paths <= 256:
        fail("--max-paths must be between 1 and 256")
    if a.pair_bp < 100:
        fail("--pair-bp must be at least 100")
    try:
        import numpy  # noqa: F401
    except ImportError:
        fail("numpy is missing. Install it (pip install numpy), or in a source folder run install.sh.")
    if a.strict or a.min_link is not None:
        sys.stderr.write("argjacent: note: --strict and --min-link were options of the caller of version 0.2; they have no effect now\n")

    genes = read_genes(a.genes)
    tmp = a.out + ".tmp"
    mark = os.path.join(tmp, ".made-by-argjacent")

    def tidy():
        if not a.keep and os.path.exists(mark):
            shutil.rmtree(tmp, ignore_errors=True)

    def on_term(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, on_term)
    try:
        out_dir = os.path.dirname(os.path.abspath(a.out))
        os.makedirs(out_dir, exist_ok=True)
        if os.path.lexists(tmp):
            if not os.path.exists(mark):                   # not a left-over of ours: never delete what we did not make
                fail("%s exists and was not made by argjacent. Move it away or choose another output prefix." % tmp)
            shutil.rmtree(tmp)
        os.makedirs(os.path.join(tmp, "c"))
        os.makedirs(os.path.join(tmp, "v"))
        open(mark, "w").close()
        for ext in (".summary.tsv", ".contexts.tsv", ".flanks.fa", ".flanks.json", ".gfa", ".pairs.tsv", ".near_1.fq", ".near_2.fq"):
            if os.path.exists(a.out + ext):                # results of an earlier run with this prefix: not to be taken for new ones
                os.remove(a.out + ext)
        return run(a, bins, genes, tmp)
    except SystemExit:
        tidy()
        raise
    except KeyboardInterrupt:
        tidy()
        sys.stderr.write("argjacent: stopped by the user\n")
        return 130
    except OSError as e:
        tidy()
        if e.errno == errno.ENOSPC:
            fail("the disk is full (while writing %s). The temporary files in %s need about as much room as the reads "
                 "near the genes take uncompressed; free some space or use an output prefix on another disk."
                 % (e.filename or "a file", tmp))
        fail("cannot read or write %s: %s" % (e.filename or "a file", e.strerror or e))
    except MemoryError:
        tidy()
        fail("out of memory. Try fewer threads (-t); for the first pass, -w 40 makes a smaller map of the sample.")


def run(a, bins, genes, tmp):
    """the work of one argjacent call; `tmp` exists"""
    from . import gfa, pair, pipeline, verify
    names = [n for n, _ in genes]
    short = [n for n, s in genes if len(s) < a.k]
    pre = os.path.join(tmp, "aw")
    log_path = a.out + ".log"
    bits = a.b if a.b else auto_bloom_bits([a.r1, a.r2], a.w)
    budget = a.max_mem * 1024.0 if a.max_mem else ((machine_mb() or 0.0) * 0.6 or None)        # MB
    t_start = time.time()

    print("argjacent %s: %d genes, %d bp flanks" % (VERSION, len(genes), a.flank))
    if short:
        print("  note: %d genes are shorter than k (%d bp) and cannot be found: %s" % (len(short), a.k, ", ".join(short[:5])))
    print("[1/4] Reading the reads %s (map, walk, fetch) ..." % ("once, map from file" if a.map else "twice"), flush=True)
    no_bloom = os.environ.get("ARGJ_NO_BLOOM") == "1"
    cmd = [bins["aw"], "-k", str(a.k), "-w", str(a.w), "-D", str(a.flank)]
    if no_bloom:
        cmd += ["-B"]
    else:
        cmd += ["-b", str(bits)]
    cmd += ["-p", str(a.min_share), "-M", str(a.max_anchors), "-t", str(a.threads)]
    if budget and not a.map:
        # the map and its Bloom filter may take the memory allowed but for a tenth of it (256 MB at most), which is left
        # for the reads on their way through; a map that does not fit is built in parts
        if not no_bloom:
            while bits > 24 and (1 << (bits - 3)) / 1048576.0 > 0.35 * budget and not a.b:
                bits -= 1                                       # a filter that would take over a third of the memory: a smaller one
            cmd[cmd.index("-b") + 1] = str(bits)
        cmd += ["-m", "%.1f" % (budget - min(0.1 * budget, 256.0))]
    if os.environ.get("ARGJ_AW_EXTRA"):
        import shlex
        cmd += shlex.split(os.environ["ARGJ_AW_EXTRA"])
    if os.environ.get("AW_TEST_PARTS"):                    # tests only: build the map in this many parts
        cmd += ["-P", os.environ["AW_TEST_PARTS"]]
    if a.map:
        cmd += ["-G", a.map]
    if a.save_map:
        cmd += ["-S", a.save_map]
    if a.reads_out:
        cmd += ["-f"]
    cmd += ["-g", a.genes, "-o", pre, a.r1, a.r2]
    with open(log_path, "w") as log:
        log.write("command: %s\n\n" % " ".join(cmd))
        log.flush()
        rc = subprocess.call(cmd, stderr=log)
    with open(log_path, errors="replace") as log:
        lines = log.readlines()
    for w in [l.rstrip() for l in lines if l.startswith("aw: ") or l.startswith("    ")]:
        print("  " + w)
    if rc != 0:
        fail("aw stopped with an error (exit %d). See %s" % (rc, log_path))
    capped = sum(1 for l in lines if l.startswith("walk ") and "(CAPPED)" in l)
    if capped:
        print("  note: %d gene sides hit the limit of %d anchors (--max-anchors); their walks are cut short" % (capped, a.max_anchors))
    st = parse_stats(pre + ".stats.txt")
    t_aw = time.time() - t_start
    pairs = int(st.get("reads", "0")) // 2
    # a map built in parts is read from its file by the walk; those pages are the system's to take back, and are not counted
    print("  %d read pairs; kept %s near the genes. %.1f s, peak memory %s MB" % (
        pairs, st.get("recruited_pairs", "?"), t_aw, st.get("peak_rss_mb_after_pass1" if st.get("map_from_file") == "1" else "peak_rss_mb", "?")))
    if a.reads_out:
        for m in ("1", "2"):
            shutil.move("%s_%s.fq" % (pre, m), "%s.near_%s.fq" % (a.out, m))
    present = []
    with open(pre + ".genes.tsv") as f:
        f.readline()
        for i, line in enumerate(f):
            if int(line.rstrip("\n").split("\t")[3]) > 0:
                present.append(i)
    read_len = int(st.get("read_len_median", "0") or 0)
    reads = pre + ".reads.tsv"

    print("[2/4] Proposing contexts (up to %d per gene side) ..." % a.max_paths, flush=True)
    t1 = time.time()
    cands, mst = {}, {}
    kept = int(st.get("kept_reads", "0") or 0)
    # awx holds the kept reads (a byte per base), 16 bytes per read and an index of 14 bytes per 8 bases of read
    need = kept * (float(st.get("read_len_mean", "150") or 150) * (1 + 14.0 / 8) + 16 + 28) / 1048576.0 + 100
    if budget and need > budget:
        print("  note: this step holds all reads kept near the genes and an index of them: about %.1f GB here, more than the "
              "%.1f GB the run was given (--max-mem)" % (need / 1024.0, budget / 1024.0))
    if present:
        seeds = os.path.join(tmp, "seeds.fa")
        pipeline.write_seeds(genes, present, seeds)
        with open(log_path, "a") as log:
            log.write("\n")
            rc = pipeline.propose(bins["awx"], seeds, os.path.join(tmp, "x"), reads, a.max_paths, a.flank, a.threads, log)
        if rc != 0:
            fail("awx stopped with an error (exit %d%s). See %s" % (
                rc, ": out of memory? It holds all reads kept near the genes" if rc < 0 or rc == 137 else "", log_path))
        cands, mst = pipeline.merge(genes, os.path.join(tmp, "x"), os.path.join(tmp, "c"))
    t_prop = time.time() - t1
    print("  %d candidates on %d gene sides (%d sides start from the sample's own version of the gene end). %.1f s" % (
        mst.get("candidates", 0), mst.get("sides_sample_end", 0) + mst.get("sides_database_end", 0), mst.get("sides_sample_end", 0), t_prop))

    print("[3/4] Laying the reads on the candidates and judging how far each is certain ...", flush=True)
    t2 = time.time()
    rc_ = verify.read_length_constants(read_len if read_len else 150)
    if read_len and read_len < 70:
        print("  note: the reads are %d bp long: too short to show where other DNA joins a context; nothing is called certain" % read_len)
    elif read_len and read_len < 140:
        print("  note: the reads are %d bp long. The rules were held out with reads of 150 bp and more; for shorter reads "
              "their read-length terms are scaled, and that has not been held out" % read_len)
    ver, failed = {}, {}
    if cands:
        keep_rec = bool(os.environ.get("AW_KEEP_RECORDS"))
        no_dedup = keep_rec or bool(os.environ.get("AW_TEST_FAIL"))
        canon_of, uniq_genes, sig_to_gene = {}, [], {}
        for i in sorted(cands):
            if no_dedup:
                canon_of[i] = i
                uniq_genes.append(i)
            else:
                sig = tuple((c["side"], c["cand"], c["gene_part"], c["seq"]) for c in cands[i])
                c0 = sig_to_gene.get(sig)
                if c0 is None:
                    sig_to_gene[sig] = i
                    canon_of[i] = i
                    uniq_genes.append(i)
                else:
                    canon_of[i] = c0
        ctx = multiprocessing.get_context("spawn")
        workers = max(1, min(a.threads, len(uniq_genes)))
        pool = concurrent.futures.ProcessPoolExecutor(max_workers=workers, mp_context=ctx)
        died = "a worker process died while judging the candidates of %s (out of memory?). Try fewer threads (-t)."
        vdir = os.path.join(tmp, "v")
        vtdir = os.path.join(tmp, "vt")
        os.makedirs(vdir, exist_ok=True)
        os.makedirs(vtdir, exist_ok=True)
        try:
            running = {}
            with open(log_path, "a") as log:
                log.write("\n")
                for i in pipeline.records(bins["awv"], os.path.join(tmp, "c"), vdir, uniq_genes, reads, a.threads, log, vtdir=vtdir):
                    try:
                        fut = pool.submit(rules_both, (i, os.path.join(vdir, str(i)), os.path.join(vtdir, str(i)),
                                                       read_len if read_len else 150, keep_rec))
                    except concurrent.futures.process.BrokenProcessPool:
                        fail(died % ", ".join(sorted(names[x] for x in running.values()) or [names[i]]))
                    running[fut] = i
            for fut in concurrent.futures.as_completed(list(running)):
                i = running[fut]
                try:
                    _, merged = fut.result()
                    ver[i] = merged
                except concurrent.futures.process.BrokenProcessPool:
                    fail(died % names[i])
                except (OSError, MemoryError):
                    raise
                except Exception as e:                     # one gene must not take the others down with it
                    failed[names[i]] = "%s: %s" % (type(e).__name__, e)
            pool.shutdown()
            for i, c0 in canon_of.items():
                if i == c0:
                    continue
                if names[c0] in failed:
                    failed[names[i]] = failed[names[c0]]
                elif c0 in ver:
                    p_old, p_new = "%d|" % c0, "%d|" % i
                    ver[i] = {p_new + k[len(p_old):]: dict(v, name=p_new + v["name"][len(p_old):])
                              for k, v in ver[c0].items()}
        except RuntimeError as e:                          # awv failed
            pool.shutdown(wait=False)
            fail("%s. See %s" % (e, log_path))
        except BaseException:                              # Ctrl-C, a dead worker, a full disk: do not wait for the other genes
            procs = list(getattr(pool, "_processes", {}).values())
            try:
                pool.shutdown(wait=False, cancel_futures=True)
            except TypeError:                              # Python 3.8
                pool.shutdown(wait=False)
            for pr in procs:
                try:
                    pr.terminate()
                except Exception:
                    pass
            raise
    for g, why in sorted(failed.items()):
        print("  warning: judging the candidates failed for %s (%s); it is marked 'failed' in the summary" % (g, why))
    og = other_genes(bins["aw"], a.k, a.w, a.genes, names, cands, tmp) if len(genes) > 1 and cands else {}
    t_ver = time.time() - t2
    print("  %.1f s" % t_ver)

    print("[4/4] Writing results ...", flush=True)
    cut = verify.CUT
    no_claims = bool(read_len) and read_len < 70
    res, counts = {}, {}
    n_ctx = n_sure = 0
    kb = min(1000, a.flank - cut)           # a whole flank that is certain is given as its length less the cut
    for i in present:
        sides = {"left": [], "right": []}
        for c in cands.get(i, []):
            v = ver.get(i, {}).get("%d|%s|%d|%d|0|0" % (i, c["side"], c["cand"], c["gene_part"]))
            fl = c["seq"][c["gene_part"]:]
            term_ok = (v is not None and not no_claims
                       and c.get("gene_end_diffs_near", 0) <= pipeline.MAXNEAR
                       and c.get("gene_end_diffs_tip", 0) <= pipeline.MAXTIP)
            certain = max(0, min(v["verified"], len(fl)) - cut) if term_ok else 0
            # The end claim says "a context of a gene that ENDS like this one": only for a candidate whose own last 100 gene bases
            # are like the gene's. A sample's gene end may differ from the gene in up to 30 of its last 300 bases; where those
            # crowd into the last 100 (or >= 10 into the terminal 30-bp boundary, e.g. cfxA6 vs cfxA5 C-terminal divergence),
            # the candidate carries a relative's end (two-strain E. coli case, 20 Nov 2026: a relative 83 % like the gene over
            # its last 100 bases had its flank called certain for 1,821 bases under the end claim).
            end_ok = term_ok and c.get("gene_end_diffs_near", 0) <= pipeline.NEARDIFF
            end = max(certain, max(0, min(v["verified_end"], len(fl)) - cut)) if end_ok else certain
            p = {"flank": fl if c["side"] == "right" else fl.translate(COMP)[::-1], "flank_len": len(fl), "context_len": certain,
                 "end_len": end, "shared_end": v["shared_end"] if v else 0,
                 "end_stop": (v["stop_end"] if not no_claims else "reads_too_short") if v else "failed",
                 "bridged": v["bridged"] if v else 0,
                 "depth": round(v["depth"], 1) if v else 0.0,
                 "gene_depth": round(v["gene_depth"], 1) if v else 0.0,
                 "depth_bins": v["depth_bins"] if v else [],
                 "stop": STOPS.get(c["stop"], c["stop"].replace(" ", "_")),
                 "certain_stop": (v["stop"] if not no_claims else "reads_too_short") if v else "failed",
                 "main": c["cand"] == 1, "gene_end": c["gene_end"], "gene_end_diffs": c["gene_end_diffs"],
                 "gene_end_diffs_near": c.get("gene_end_diffs_near", 0),
                 "joins": v["joins"] if v else 0, "open_joins": v["open_joins"] if v else 0,
                 "genes": [dict(x, certain=x["to_bp"] <= certain) for x in og.get((i, c["side"], c["cand"]), [])]}
            sides[c["side"]].append(p)
        res[names[i]] = {}
        for sd in ("left", "right"):
            paths = sorted(sides[sd], key=lambda p: (-p["context_len"], -p["end_len"], -p["flank_len"]))
            for rank, p in enumerate(paths, 1):
                p["id"] = "c%d" % rank
            res[names[i]][sd] = {"paths": paths, "read_depth": max((p["depth"] for p in paths), default=0.0),
                                 "gene_end": paths[0]["gene_end"] if paths else "-"}
    with open(a.out + ".summary.tsv", "w") as tsv, open(a.out + ".contexts.tsv", "w") as ctsv, open(a.out + ".flanks.fa", "w") as fa:
        tsv.write("gene\tside\tstatus\tcontexts\tcertain_bp\tlongest_bp\tcontexts_certain_1kb\tread_depth\tgene_end\tcertain_end_bp\n")
        ctsv.write("gene\tside\tcontext\tlength_bp\tcertain_bp\tread_depth\tstop\tcertain_stop\tother_genes\tgene_end\tmain\tcertain_end_bp\tshared_end_bp\n")
        for name, _ in genes:
            for sd in ("left", "right"):
                if name in failed:
                    stt = "failed"
                elif name not in res:
                    stt = "not_found"
                else:
                    stt, best = status_of(res[name][sd]["paths"], a.flank, cut)
                    res[name][sd]["status"] = stt
                counts[stt] = counts.get(stt, 0) + 1
                if stt in ("failed", "not_found") or not res[name][sd]["paths"]:
                    tsv.write("%s\t%s\t%s\t0\t0\t0\t0\t0\t-\t0\n" % (name, sd, stt))
                    continue
                side = res[name][sd]
                paths = side["paths"]
                tsv.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.0f\t%s\t%d\n" % (
                    name, sd, stt, len(paths), best["context_len"], max(p["flank_len"] for p in paths),
                    sum(1 for p in paths if p["context_len"] >= kb), side["read_depth"], side["gene_end"],
                    max(p["end_len"] for p in paths)))
                for p in paths:
                    n_ctx += 1
                    n_sure += p["context_len"] >= kb
                    ogs = ",".join("%s[%s]@%d-%d(%d%%%s)" % (x["gene"], x["strand"], x["from_bp"], x["to_bp"],
                                                             round(100.0 * x["anchors"] / x["of"]), "" if x["certain"] else ",uncertain")
                                   for x in p["genes"]) or "-"
                    ctsv.write("%s\t%s\t%s\t%d\t%d\t%.0f\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n" % (
                        name, sd, p["id"], p["flank_len"], p["context_len"], p["depth"], p["stop"], p["certain_stop"], ogs, p["gene_end"],
                        "yes" if p["main"] else "-", p["end_len"], p["shared_end"]))
                    fa.write(">%s|%s|%s len=%d certain=%d depth=%.0f stop=%s%s%s\n" % (
                        name, sd, p["id"], p["flank_len"], p["context_len"], p["depth"], p["stop"],
                        " certain_end=%d" % p["end_len"] if p["end_len"] > p["context_len"] else "", " main=yes" if p["main"] else ""))
                    for j in range(0, len(p["flank"]), 80):
                        fa.write(p["flank"][j:j + 80] + "\n")
    doc_genes = {name: res[name] for name, _ in genes if name in res and name not in failed}
    genes_dict = dict(genes)
    pair_w = min(a.pair_bp, max(100, a.flank - cut))
    pairs_by_gene = {}
    if not a.no_pair:
        if doc_genes and not no_claims:
            reads_src = reads if os.path.isfile(reads) else None
            pairs_by_gene = pair.pair_sample(genes_dict, doc_genes, W=pair_w, reads_source=reads_src)
            for gname, prs in pairs_by_gene.items():
                doc_genes[gname]["pairs"] = prs
        pair.write_pairs_tsv(pairs_by_gene, a.out + ".pairs.tsv")

    with open(a.out + ".flanks.json", "w") as f:
        json.dump({"version": VERSION, "flank_bp": a.flank, "aw_stats": st,
                   "parameters": {"k": a.k, "w": a.w, "min_share": a.min_share, "max_paths": a.max_paths, "read_len": read_len,
                                  "read_len_used_by_rules": rc_["L"], "certain_cut_bp": cut, "pair_bp": pair_w},
                   # genes in the order of the input file
                   "genes": doc_genes}, f, indent=1)

    gfa_st = None
    if not a.no_gfa:
        gfa_st = gfa.write_gfa(
            genes_dict, doc_genes, a.out + ".gfa",
            pairs_by_gene=pairs_by_gene, verify=True,
            read_len=(read_len if read_len else 150),
        )

    n_pairs = sum(len(v) for v in pairs_by_gene.values())
    with open(log_path, "a") as log:
        log.write("\nproposing: %.1f s (%s)\njudging: %.1f s\nstatus counts: %s\n" % (t_prop, json.dumps(mst), t_ver, json.dumps(counts)))
        if not a.no_pair:
            log.write("pairing: %d certified two-sided pairs across %d genes (W=%d bp)\n" % (n_pairs, len(pairs_by_gene), pair_w))
        if gfa_st is not None:
            log.write("gfa: %s\n" % json.dumps(gfa_st))
    if not a.keep:
        shutil.rmtree(tmp)

    order = ["complete", "partial", "ended", "not_found", "failed"]
    print("Done in %.1f s. Gene sides: %s. Contexts: %d, of which certain for %s: %d%s" % (
        time.time() - t_start, ", ".join("%s %d" % (k, counts[k]) for k in order if k in counts), n_ctx,
        "1 kb or more" if kb == 1000 else "the whole %d bp" % a.flank, n_sure,
        "" if a.no_pair else "; certified pairs: %d" % n_pairs))
    outs = [a.out + ext for ext in (".summary.tsv", ".contexts.tsv", ".flanks.fa", ".flanks.json")]
    if not a.no_pair:
        outs.append(a.out + ".pairs.tsv")
    if not a.no_gfa:
        outs.append(a.out + ".gfa")
    print("Results: " + "  ".join(outs))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
