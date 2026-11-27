#!/usr/bin/env python3
"""Extract +/-D bp flanking sequences around target genes from a pre-made metagenome assembly
(MGnify/ENA ERZ .fasta.gz, Logan .contigs.fa.zst, or uncompressed .fa/.fna).

Outputs (all oriented 5' -> 3' in the target gene's orientation, matching ARGjacent's flanks.json):
  - <out-prefix>.asm_flanks.json   : machine-readable flanks (compatible with compare_argjacent_vs_assembly.py)
  - <out-prefix>.asm_contexts.tsv  : per-gene-side hit table (contig, contig_len, pident, qcov, flank_bp, stop)
  - <out-prefix>.asm_flanks.fa     : extracted flanking DNA FASTA ready for downstream annotation

Supports three alignment backends via --aligner auto|minimap2|blastn|kmer:
  - blastn   : makeblastdb + blastn -task megablast (highest nucleotide sensitivity, matches score_cami.py)
  - minimap2 : minimap2 -c -x asm20 --secondary=yes -N 50 (very fast on large metagenome assemblies)
  - kmer     : zero-dependency pure-Python 21-mer seeded matcher (fallback when neither tool is on PATH)
"""
import argparse
import collections
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
MERGE_BP = 300


def rc(seq):
    return seq.translate(COMP)[::-1]


def open_text_stream(path):
    """Yield lines from plain, .gz, or .zst FASTA without extracting huge assemblies to disk."""
    if path.endswith(".gz"):
        fh = gzip.open(path, "rt", errors="replace")
        return fh, None
    if path.endswith(".zst"):
        zstd = shutil.which("zstd")
        if not zstd:
            sys.exit(f"extract_assembly_flanks: 'zstd' required to read {path}")
        proc = subprocess.Popen([zstd, "-dcf", path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return proc.stdout, proc
    return open(path, "r", errors="replace"), None


def iter_fasta(path):
    fh, proc = open_text_stream(path)
    name, buf = None, []
    try:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper()
                name = line[1:].strip().split()[0]
                buf = []
            else:
                s = line.strip()
                if s:
                    buf.append(s)
        if name is not None:
            yield name, "".join(buf).upper()
    finally:
        fh.close()
        if proc is not None:
            err = proc.stderr.read() if proc.stderr else ""
            rc_code = proc.wait()
            if rc_code != 0:
                sys.exit(f"extract_assembly_flanks: zstd decompression failed for {path} (rc={rc_code}): {err.strip()}")


def read_fasta(path, want=None):
    out = {}
    for name, seq in iter_fasta(path):
        if want is None or name in want:
            out[name] = seq
    return out


def materialize_plain_fasta(path, tmp_dir):
    """Decompress .gz or .zst into tmp_dir when needed by makeblastdb."""
    if not (path.endswith(".gz") or path.endswith(".zst")):
        return path
    dst = os.path.join(tmp_dir, "contigs.fa")
    if path.endswith(".gz"):
        with gzip.open(path, "rb") as src, open(dst, "wb") as out:
            shutil.copyfileobj(src, out, length=4 * 1024 * 1024)
    else:
        zstd = shutil.which("zstd")
        if not zstd:
            sys.exit(f"extract_assembly_flanks: 'zstd' required to decompress {path}")
        with open(dst, "wb") as out:
            subprocess.run([zstd, "-dcf", path], stdout=out, check=True)
    return dst


def loci_from_blast(genes_fa, contigs_path, glen, min_id, min_cov, threads, tmp_dir):
    plain = materialize_plain_fasta(contigs_path, tmp_dir)
    db = os.path.join(tmp_dir, "ctg_db")
    subprocess.run(
        ["makeblastdb", "-in", plain, "-dbtype", "nucl", "-out", db],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    out_tsv = os.path.join(tmp_dir, "blast.tsv")
    subprocess.run(
        [
            "blastn", "-task", "megablast", "-query", genes_fa, "-db", db,
            "-perc_identity", str(min_id), "-max_target_seqs", "200",
            "-num_threads", str(threads),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
            "-out", out_tsv,
        ],
        stderr=subprocess.DEVNULL, check=True,
    )

    by = collections.defaultdict(list)
    with open(out_tsv) as fh:
        for line in fh:
            if not line.strip():
                continue
            q, s, pid, alen, mm, gap, qs, qe, ss, se, ev, bit, ql, sl = line.rstrip("\n").split("\t")
            qs, qe, ss, se = int(qs), int(qe), int(ss), int(se)
            if q not in glen or float(pid) < min_id or (qe - qs + 1) < 0.10 * glen[q]:
                continue
            by[(q, s, ss <= se)].append((min(ss, se), max(ss, se), float(bit), float(pid), qs, qe, ss, se, int(alen)))

    out = collections.defaultdict(list)
    for (q, s, plus), hs in by.items():
        hs.sort(key=lambda h: -h[2])
        kept = []
        for h in hs:
            if all(min(h[1], k[1]) - max(h[0], k[0]) < 0.5 * max(1, h[1] - h[0]) for k in kept):
                kept.append(h)
        kept.sort()
        groups = []
        for h in kept:
            g = groups[-1] if groups else None
            if g and h[0] - max(x[1] for x in g) <= MERGE_BP and (
                h[4] > max(x[5] for x in g) - 50 if plus else h[5] < min(x[4] for x in g) + 50
            ):
                g.append(h)
            else:
                groups.append([h])
        for g in groups:
            qcov = sum(x[5] - x[4] + 1 for x in g) / float(glen[q])
            if qcov < min_cov / 100.0:
                continue
            first = min(g, key=lambda x: x[4])
            last = max(g, key=lambda x: x[5])
            pid = sum(x[3] * x[8] for x in g) / max(1, sum(x[8] for x in g))
            # Convert 1-based BLAST coordinates to (contig, pid, qcov*100, qs1, qe1, ss1, se1)
            out[q].append((s, pid, min(100.0, qcov * 100.0), first[4], last[5], first[6], last[7]))
    return out


def loci_from_minimap2(genes_fa, contigs_path, glen, min_id, min_cov, threads, tmp_dir):
    # minimap2 natively reads .fa and .fa.gz; decompress .zst only if needed
    target = materialize_plain_fasta(contigs_path, tmp_dir) if contigs_path.endswith(".zst") else contigs_path
    paf = os.path.join(tmp_dir, "mm2.paf")
    with open(paf, "w") as out_fh:
        subprocess.run(
            ["minimap2", "-c", "-x", "asm20", "--secondary=yes", "-N", "100", "-p", "0.5", "-t", str(threads), target, genes_fa],
            stdout=out_fh, stderr=subprocess.DEVNULL, check=True,
        )
    by = collections.defaultdict(list)
    with open(paf) as fh:
        for line in fh:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            q, ql, qs0, qe0, strand, s, sl, ss0, se0, nmatch, blen = (
                p[0], int(p[1]), int(p[2]), int(p[3]), p[4], p[5], int(p[6]), int(p[7]), int(p[8]), int(p[9]), int(p[10])
            )
            if q not in glen or blen <= 0:
                continue
            pid = 100.0 * nmatch / blen
            qcov = 100.0 * (qe0 - qs0) / glen[q]
            if pid < min_id or qcov < min_cov:
                continue
            # Convert PAF 0-based [qs0, qe0), [ss0, se0) into 1-based BLAST-like endpoints
            qs1, qe1 = qs0 + 1, qe0
            if strand == "+":
                ss1, se1 = ss0 + 1, se0
            else:
                ss1, se1 = se0, ss0 + 1
            by[(q, s)].append((pid * qcov, s, pid, min(100.0, qcov), qs1, qe1, ss1, se1))

    out = collections.defaultdict(list)
    for (q, s), hits in by.items():
        hits.sort(reverse=True)
        kept = []
        for h in hits:
            a0, a1 = min(h[6], h[7]), max(h[6], h[7])
            if all(min(a1, max(k[6], k[7])) - max(a0, min(k[6], k[7])) < 0.5 * max(1, a1 - a0) for k in kept):
                kept.append(h)
        for h in kept:
            out[q].append(h[1:])
    return out


def loci_from_kmer(genes, contigs_path, min_id, min_cov):
    """Pure-Python 21-mer seeded matcher (fallback when blastn/minimap2 are unavailable)."""
    K = 21
    STEP = 7
    index = collections.defaultdict(list)
    for g, seq in genes.items():
        L = len(seq)
        if L < K:
            continue
        for pos in range(0, L - K + 1, STEP):
            kmer = seq[pos:pos + K]
            if "N" not in kmer:
                index[kmer].append((g, pos, True))
                index[rc(kmer)].append((g, pos, False))

    out = collections.defaultdict(list)
    for cname, cseq in iter_fasta(contigs_path):
        clen = len(cseq)
        if clen < K:
            continue
        # Vote on (gene, strand, projected_gene_start_on_contig)
        votes = collections.defaultdict(int)
        for pos in range(0, clen - K + 1, 5):
            hits = index.get(cseq[pos:pos + K])
            if not hits:
                continue
            for g, gpos, plus in hits:
                # Where gene index 0 falls on contig (+ strand: pos - gpos; - strand: pos + gpos + K)
                diag = (pos - gpos) if plus else (pos + gpos + K)
                votes[(g, plus, diag // 16)] += 1

        if not votes:
            continue
        # Group nearby diagonals per (gene, strand)
        cand_by_gs = collections.defaultdict(list)
        for (g, plus, dbin), cnt in votes.items():
            if cnt >= 2:
                cand_by_gs[(g, plus)].append((cnt, dbin * 16))

        for (g, plus), cands in cand_by_gs.items():
            gseq = genes[g]
            L = len(gseq)
            cands.sort(reverse=True)
            seen_pos = []
            for _, est in cands:
                if any(abs(est - sp) < L // 2 for sp in seen_pos):
                    continue
                # Refine exact offset within +/-24 bp
                best = None
                for delta in range(-24, 25):
                    if plus:
                        g0 = est + delta
                        g1 = g0 + L
                        o0, o1 = max(0, g0), min(clen, g1)
                        if o1 - o0 < (min_cov / 100.0) * L:
                            continue
                        sub = cseq[o0:o1]
                        qsub = gseq[o0 - g0:o1 - g0]
                    else:
                        g1 = est + delta
                        g0 = g1 - L
                        o0, o1 = max(0, g0), min(clen, g1)
                        if o1 - o0 < (min_cov / 100.0) * L:
                            continue
                        sub = rc(cseq[o0:o1])
                        qsub = gseq[g1 - o1:g1 - o0]
                    matches = sum(1 for a, b in zip(sub, qsub) if a == b)
                    pid = 100.0 * matches / len(sub)
                    if best is None or matches > best[0]:
                        best = (matches, pid, 100.0 * len(sub) / L, o0, o1, g0, g1)
                if best and best[1] >= min_id and best[2] >= min_cov:
                    _, pid, qcov, o0, o1, g0, g1 = best
                    seen_pos.append(est)
                    if plus:
                        qs1 = (o0 - g0) + 1
                        qe1 = (o1 - g0)
                        ss1, se1 = o0 + 1, o1
                    else:
                        qs1 = (g1 - o1) + 1
                        qe1 = (g1 - o0)
                        ss1, se1 = o1, o0 + 1
                    out[g].append((cname, pid, min(100.0, qcov), qs1, qe1, ss1, se1))
    return out


def project_5to3(contig_seq, gene_len, qs1, qe1, ss1, se1, flank_bp):
    """Return (left_5to3, right_5to3, g0, g1) in the target gene's 5'->3' orientation.

    Matches ARGjacent's flanks.json orientation exactly:
      [left_5to3][GENE (5'->3')][right_5to3]
    """
    n = len(contig_seq)
    if ss1 <= se1:
        # Gene lies on '+' strand of contig
        g0 = ss1 - qs1                      # 0-based start of full gene
        g1 = se1 + (gene_len - qe1)         # 0-based exclusive end of full gene
        left = contig_seq[max(0, g0 - flank_bp):max(0, g0)]
        right = contig_seq[min(n, g1):min(n, g1 + flank_bp)]
    else:
        # Gene lies on '-' strand of contig
        g0 = se1 - 1 - (gene_len - qe1)
        g1 = ss1 + (qs1 - 1)
        left = rc(contig_seq[min(n, g1):min(n, g1 + flank_bp)])
        right = rc(contig_seq[max(0, g0 - flank_bp):max(0, g0)])
    return left, right, g0, g1


def choose_aligner(requested):
    if requested != "auto":
        if requested in ("blastn", "minimap2") and not shutil.which(requested):
            sys.exit(f"extract_assembly_flanks: requested aligner '{requested}' not found on PATH")
        return requested
    if shutil.which("minimap2"):
        return "minimap2"
    if shutil.which("blastn") and shutil.which("makeblastdb"):
        return "blastn"
    return "kmer"


def main():
    ap = argparse.ArgumentParser(description="Extract +/-D bp gene flanks from a pre-made metagenome assembly.")
    ap.add_argument("-g", "--genes", required=True, help="Target genes FASTA (same genes.fa given to argjacent)")
    ap.add_argument("-c", "--contigs", required=True, help="Pre-made assembly FASTA (.fa, .fa.gz, or .fa.zst)")
    ap.add_argument("-o", "--out", required=True, help="Output prefix (writes .asm_flanks.json, .asm_contexts.tsv, .asm_flanks.fa)")
    ap.add_argument("-D", "--flank-bp", type=int, default=10000, help="Flank length to extract per side (default: 10000)")
    ap.add_argument("--min-id", type=float, default=90.0, help="Minimum gene alignment identity %% (default: 90.0)")
    ap.add_argument("--min-cov", type=float, default=60.0, help="Minimum gene alignment coverage %% (default: 60.0)")
    ap.add_argument("--aligner", choices=("auto", "blastn", "minimap2", "kmer"), default="auto",
                    help="Alignment backend (default: auto -> minimap2 > blastn > kmer)")
    ap.add_argument("-t", "--threads", type=int, default=4, help="Threads for minimap2/blastn (default: 4)")
    ap.add_argument("--source-label", default="assembly", help="Metadata tag stored in JSON (e.g. 'ERZ12345' or 'logan')")
    args = ap.parse_args()

    genes = read_fasta(args.genes)
    if not genes:
        sys.exit(f"extract_assembly_flanks: no sequences found in {args.genes}")
    glen = {g: len(s) for g, s in genes.items()}
    aligner = choose_aligner(args.aligner)

    with tempfile.TemporaryDirectory(prefix="argj_asm_") as tmp_dir:
        if aligner == "blastn":
            loci = loci_from_blast(args.genes, args.contigs, glen, args.min_id, args.min_cov, args.threads, tmp_dir)
        elif aligner == "minimap2":
            loci = loci_from_minimap2(args.genes, args.contigs, glen, args.min_id, args.min_cov, args.threads, tmp_dir)
        else:
            loci = loci_from_kmer(genes, args.contigs, args.min_id, args.min_cov)

        want_contigs = {h[0] for hs in loci.values() for h in hs}
        contig_seqs = read_fasta(args.contigs, want=want_contigs) if want_contigs else {}

    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    json_genes = {}
    tsv_rows = []
    fa_records = []
    n_hits = 0

    for g in sorted(genes):
        hits = loci.get(g, [])
        if not hits:
            continue
        sides = {"left": [], "right": []}
        # Sort deterministically by (max flank, pident, contig)
        enriched = []
        for cname, pid, qcov, qs1, qe1, ss1, se1 in hits:
            cseq = contig_seqs.get(cname, "")
            if not cseq:
                continue
            left, right, g0, g1 = project_5to3(cseq, glen[g], qs1, qe1, ss1, se1, args.flank_bp)
            enriched.append((max(len(left), len(right)), pid, qcov, cname, len(cseq), left, right))
        enriched.sort(key=lambda x: (-x[0], -x[1], x[3]))

        for idx, (_, pid, qcov, cname, clen, left, right) in enumerate(enriched, 1):
            n_hits += 1
            both_1kb = (len(left) >= 1000 and len(right) >= 1000)
            for side, fseq in (("left", left), ("right", right)):
                stop = "reached_target" if len(fseq) >= args.flank_bp else "contig_end"
                cid = f"a{idx}"
                sides[side].append({
                    "id": cid,
                    "flank": fseq,
                    "flank_len": len(fseq),
                    "context_len": len(fseq),  # assemblers emit no certainty boundary; full contig flank is claimed
                    "stop": stop,
                    "contig": cname,
                    "contig_len": clen,
                    "pident": round(pid, 2),
                    "qcov": round(qcov, 2),
                    "has_both_1kb": both_1kb,
                })
                tsv_rows.append((
                    g, side, cid, cname, clen, f"{pid:.2f}", f"{qcov:.2f}",
                    len(fseq), stop, "1" if both_1kb else "0",
                ))
                if fseq:
                    fa_records.append((f"{g}|{side}|{cid}|len={len(fseq)}|ctg={cname}|stop={stop}", fseq))

        json_genes[g] = {
            "left": {"paths": sides["left"]},
            "right": {"paths": sides["right"]},
        }

    with open(f"{args.out}.asm_flanks.json", "w") as fh:
        json.dump({
            "version": "asm_flanks_v1",
            "source": args.source_label,
            "aligner": aligner,
            "flank_bp": args.flank_bp,
            "min_id": args.min_id,
            "min_cov": args.min_cov,
            "genes": json_genes,
        }, fh)

    with open(f"{args.out}.asm_contexts.tsv", "w") as fh:
        fh.write("gene\tside\tcontext\tcontig\tcontig_len\tpident\tqcov\tflank_bp\tstop\thas_both_1kb\n")
        for r in tsv_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    with open(f"{args.out}.asm_flanks.fa", "w") as fh:
        for hdr, seq in fa_records:
            fh.write(f">{hdr}\n{seq}\n")

    sys.stderr.write(
        f"extract_assembly_flanks ({aligner}): {len(json_genes)}/{len(genes)} genes matched in "
        f"{len(want_contigs)} contigs ({n_hits} total gene copies).\n"
    )


if __name__ == "__main__":
    main()
