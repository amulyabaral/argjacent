"""From the reads kept near the genes to candidate contexts and the records the rules read.

  propose   the C program awx extends every gene end, base by base, with the reads that match the end of the path, and
            follows up to N arms where they disagree. A read that carries another version of a few bases (one base, a small
            insertion or deletion, a patch the arms come back from) stays on the path and is no arm of its own. It runs
            twice in one go: starting from the gene end as given, and with the gene's last 300 bases taken from the
            sample's own reads.
  merge     per gene side, the candidates that start from the sample's own gene end are taken if that end is at least
            90 % identical to the gene as given; else those that start from the gene as given.
  records   the C program awv lays all kept reads on a gene's candidates and writes, per position, what the reads say
            (depth of matching reads, reads that come on or leave, with what their mates do). It is run twice: with the
            strict measure (a read matches without any difference) and with the tolerant one (TOLERANT below).

Inside the temporary folder a gene is known by its number in the gene file, so that no gene name can break a file
name or a field separator.
"""
import os
import subprocess

COMP = str.maketrans("ACGTN", "TGCAN")
SEED = 400        # a candidate begins with the gene's last 400 bases (or the whole gene, if shorter)
REDO = 300        # of which the last 300 may be the sample's own
MAXDIFF = 30      # ... if at most 30 of them differ from the gene as given (90 % identity)
NEAR = 100        # the end claim ("a context of a gene that ENDS like this one") is made only for a candidate whose own last
NEARDIFF = 10     # NEAR bases of the gene part differ from the gene as given in at most NEARDIFF places
TOLERANT = ["-T", "50", "-B", "1", "-I", "12", "-g", "30", "-w", "9"]   # the tolerant measure of awv: single differences, small patches, one small insertion or deletion


def write_seeds(genes, present, path):
    """genes: [(name, sequence)]; present: numbers of the genes found in the reads -> the gene ends, running outward"""
    with open(path, "w") as f:
        for i in present:
            s = genes[i][1].upper()
            g = min(SEED, len(s))
            f.write(">%d|right|%d\n%s\n>%d|left|%d\n%s\n" % (i, g, s[-g:], i, g, s[:g].translate(COMP)[::-1]))


def propose(awx, seeds, prefix, reads, max_cands, flank, threads, log, stdin=None):
    """run awx; writes PREFIX.r0.fa, PREFIX.r0.tsv, PREFIX.r300.fa, PREFIX.r300.tsv -> exit status"""
    cmd = [awx, "-c", str(max_cands), "-D", str(flank), "-r", "0,%d" % REDO, "-t", str(threads), seeds, prefix, reads]
    log.write("command: %s\n" % " ".join(cmd))
    log.flush()
    return subprocess.call(cmd, stderr=log, stdin=stdin)


def _read_run(prefix, redo):
    """the candidates of one run -> {gene number: [(side, number in the run, gene part, sequence, row of the table)]}"""
    info = {}
    with open("%s.r%d.tsv" % (prefix, redo)) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            t = dict(zip(head, line.rstrip("\n").split("\t")))
            info[(t["gene"], t["side"], t["cand"])] = t
    out, name = {}, None
    with open("%s.r%d.fa" % (prefix, redo)) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip().split("|")
            else:
                i, side, n, gl = name[:4]
                out.setdefault(int(i), []).append((side, n, int(gl), line.strip(), info[(i, side, n)]))
    return out


def _end_diffs(gene_seq, side, gl, seq, last=None):
    """in how many places the candidate's gene part differs from the gene as given, over its last REDO bases (or `last`)"""
    last = REDO if last is None else last
    s = gene_seq.upper()
    g = min(gl, len(s))
    seed = s[-g:] if side == "right" else s[:g].translate(COMP)[::-1]
    part = seq[:gl]
    if len(part) != len(seed):
        return 10 ** 6
    return sum(1 for x, y in zip(part[-last:], seed[-last:]) if x != y)


def merge(genes, prefix, cdir):
    """-> ({gene number: [candidate]}, counts). Writes cdir/<gene number>.cands.fa for the record program.
    A candidate is a dict: side, cand (its number on that side, from 1), gene_part, seq (gene part + flank, running
    outward), stop, gene_end ('sample' or 'database'), gene_end_diffs."""
    plain, redo = _read_run(prefix, 0), _read_run(prefix, REDO)
    out, st = {}, {"sides_sample_end": 0, "sides_database_end": 0, "sides_sample_end_too_far": 0, "dropped_too_far": 0, "candidates": 0}
    for i in sorted(set(plain) | set(redo)):
        p, r = plain.get(i, []), redo.get(i, [])
        sides = []
        for c in p + r:
            if c[0] not in sides:
                sides.append(c[0])
        cands = []
        for side in sides:
            rs = [(c, _end_diffs(genes[i][1], side, c[2], c[3])) for c in r if c[0] == side]
            keep = [(c, d, "sample") for c, d in rs if d <= MAXDIFF]
            if keep:
                st["sides_sample_end"] += 1
                st["dropped_too_far"] += len(rs) - len(keep)
            else:
                keep = [(c, _end_diffs(genes[i][1], side, c[2], c[3]), "database") for c in p if c[0] == side]
                st["sides_database_end"] += 1
                if rs:
                    st["sides_sample_end_too_far"] += 1
            for n, (c, d, src) in enumerate(keep, 1):
                cands.append({"side": side, "cand": n, "gene_part": c[2], "seq": c[3], "stop": c[4]["stop"],
                              "gene_end": src, "gene_end_diffs": d, "run_cand": int(c[1]),
                              "gene_end_diffs_near": _end_diffs(genes[i][1], side, c[2], c[3], NEAR)})
            st["candidates"] += len(keep)
        if cands:
            out[i] = cands
            with open(os.path.join(cdir, "%d.cands.fa" % i), "w") as f:
                for c in cands:
                    f.write(">%d|%s|%d|%d|0|0\n%s\n" % (i, c["side"], c["cand"], c["gene_part"], c["seq"]))
    return out, st


def records(awv, cdir, vdir, gene_numbers, reads, threads, log, stdin=None, tolerant=False, vtdir=None):
    """run awv on the candidates of every gene; yields each gene number when its records (vdir/<n>.idx, .i32, .ent,
    .lev, .tpl) are written. When vtdir is given, runs strict (vdir) and tolerant (vtdir) together in a single pass.
    Raises RuntimeError if awv fails."""
    jobs = os.path.join(vdir, "jobs.tsv")
    order = list(gene_numbers)
    with open(jobs, "w") as f:
        for i in order:
            f.write("%s\t%s\n" % (os.path.join(cdir, "%d.cands.fa" % i), os.path.join(vdir, str(i))))
    if vtdir is not None:
        vtjobs = os.path.join(vtdir, "jobs.tsv")
        with open(vtjobs, "w") as f:
            for i in order:
                f.write("%s\t%s\n" % (os.path.join(cdir, "%d.cands.fa" % i), os.path.join(vtdir, str(i))))
        cmd = [awv, "-t", str(threads), "-q", "-Z", vtjobs] + TOLERANT + [jobs, reads]
    else:
        cmd = [awv, "-t", str(threads)] + (TOLERANT if tolerant else []) + [jobs, reads]
    log.write("command: %s\n" % " ".join(cmd))
    log.flush()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log, stdin=stdin, text=True)
    try:
        for line in proc.stdout:
            line = line.strip()
            if line.isdigit():
                yield order[int(line)]
    finally:
        proc.stdout.close()
        rc = proc.wait()
    if rc != 0:
        raise RuntimeError("awv stopped with an error (exit %d)" % rc)
