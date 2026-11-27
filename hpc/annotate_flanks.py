#!/usr/bin/env python3
"""Unified downstream annotation and head-to-head biological comparison for ARGjacent certified flanks
vs. pre-made metagenome assemblies (ENA/MGnify ERZ or Logan contigs).

Implements the exact 5-layer annotation stack from paper/global_study_protocol.md on BOTH flank sets
identically so every biological comparison (reagent vectors, IS/MGE proximity, co-carried ARGs,
plasmid vs chromosomal calls, MGE-masked host taxonomy, and cross-biome HGT clusters) is apples-to-apples:

  1. Reagent / cloning vector screen FIRST (--univec-db): flags pMB1/ColE1/pUC/pET/lacI/rop vectors
  2. Co-localized ARGs (--arg-db) + IS elements (--isfinder-db) + ORF/mobileOG-db (--mobileog-db)
     with exact bp distance from the focal ARG boundary (respecting left vs right 5'->3' orientation!)
  3. Plasmid / viral classification (--genomad-db on >= 2 kb flanks + --plasmidfinder-db replicons)
  4. Host taxonomy on MGE/IS-masked sequences (--kraken2-db or --mmseqs-tax-db)
  5. Cross-biome HGT clustering (--cluster-hgt, >= 99% id, >= 90% cov, excluding reagent vectors)

Examples:
    # Run paired comparison on outputs of compare_argjacent_vs_assembly.py --export-fasta
    python3 hpc/annotate_flanks.py \\
        --argj-fa results/comparison/global_cmp_argjacent_certified.fa \\
        --asm-fa  results/comparison/global_cmp_assembly_flanks.fa \\
        --manifest global_2023_2024.tsv \\
        --univec-db /db/univec/UniVec.fa \\
        --isfinder-db /db/isfinder/IS.fna \\
        --arg-db /db/resfinder/resfinder.fa \\
        --mobileog-db /db/mobileOG/mobileOG.dmnd \\
        --plasmidfinder-db /db/plasmidfinder/plasmidfinder.fa \\
        --genomad-db /db/genomad_db \\
        --kraken2-db /db/kraken2_gtdb \\
        --cluster-hgt -t 16 -o results/annotation

    # Run on a single certified flank FASTA only
    python3 hpc/annotate_flanks.py --fasta certified.fa --arm-name argjacent \\
        --univec-db /db/univec/UniVec.fa --isfinder-db /db/isfinder/IS.fna -o results/annotation
"""
import argparse
import collections
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}

PLASMID_MOG_HINTS = re.compile(r"\b(transfer|conjug|mob[a-z0-9]|tra[a-z]|trb[a-z]|rep[a-z0-9]|par[ab]|plasmid|oriT)\b", re.I)
IS_MOG_HINTS = re.compile(r"\b(integration/excision|transpos|integras|recombinas|resolv|ins[a-z]|tnp[a-z]|intI)\b", re.I)


def rc(seq):
    return seq.translate(COMP)[::-1]


def read_fasta(path):
    records = []
    name, desc, buf = None, "", []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, desc, "".join(buf).upper()))
                raw = line[1:].strip()
                parts = raw.split(None, 1)
                name = parts[0]
                desc = parts[1] if len(parts) > 1 else ""
                buf = []
            else:
                s = line.strip()
                if s:
                    buf.append(s)
        if name is not None:
            records.append((name, desc, "".join(buf).upper()))
    return records


def write_fasta(records, path):
    with open(path, "w") as out:
        for rec_id, seq in records:
            out.write(f">{rec_id}\n{seq}\n")


def parse_flank_header(rec_id, seq_len, manifest_meta=None):
    """Parse headers from compare_argjacent_vs_assembly.py, argjacent, or extract_assembly_flanks.py.

    Expected primary format:
        sample|gene|side|ctx_id|certain_bp=2470|biome=human_gut
    Robust even if gene names contain '|' characters (parses side from the right).
    """
    manifest_meta = manifest_meta or {}
    parts = rec_id.split("|")
    sample = parts[0] if parts else rec_id
    gene = ""
    side = "right"
    ctx_id = "c1"
    biome = ""

    # Locate 'left' or 'right' token from the right so pipes inside gene names never break parsing
    side_idx = -1
    for i in range(len(parts) - 1, 0, -1):
        if parts[i] in ("left", "right"):
            side_idx = i
            break
    if side_idx >= 2:
        side = parts[side_idx]
        gene = "|".join(parts[1:side_idx])
        if side_idx + 1 < len(parts):
            ctx_id = parts[side_idx + 1]
    elif side_idx == 1:
        # Standalone format without sample prefix: gene|side|ctx_id|...
        sample = "sample"
        gene = parts[0]
        side = parts[1]
        if len(parts) > 2:
            ctx_id = parts[2]
    else:
        gene = parts[1] if len(parts) > 1 else rec_id

    for tok in parts[side_idx + 1:] if side_idx >= 0 else parts[1:]:
        if tok.startswith("biome="):
            biome = tok.split("=", 1)[1]

    m = manifest_meta.get(sample, {})
    return {
        "flank_id": rec_id,
        "sample": sample,
        "study_accession": m.get("study_accession", ""),
        "biome_group": biome or m.get("biome_group", ""),
        "country": m.get("country", ""),
        "assembler": m.get("assembler", ""),
        "gene": gene,
        "side": side,
        "ctx_id": ctx_id,
        "flank_len": seq_len,
    }


def distance_to_arg(side, flank_len, start0, end0):
    """Exact bp distance from a feature interval [start0, end0) to the focal ARG boundary.

    Because flanks are oriented 5'->3' in the focal gene's orientation:
      - right flank: index 0 touches the 3' edge of the focal ARG -> dist = start0
      - left flank : index flank_len touches the 5' edge of the focal ARG -> dist = flank_len - end0
    """
    if side == "left":
        return max(0, flank_len - end0)
    return max(0, start0)


def parse_is_family(subject_id, subject_desc=""):
    """Extract IS element name and family from standard ISfinder headers (e.g. IS26_IS6_IS6 or IS26|IS6)."""
    s = subject_id.strip()
    if "|" in s:
        p = s.split("|")
        return p[0], (p[1] if len(p) > 1 else "unknown")
    if "_" in s:
        p = s.split("_")
        # Standard ISfinder format: ISname_ISfamily_ISgroup
        return p[0], (p[1] if len(p) > 1 else "unknown")
    m = re.search(r"(IS\w+)", subject_desc)
    return s, (m.group(1) if m else "unknown")


def parse_mobileog_target(subject_id, subject_desc=""):
    """Extract (gene_name, major_category) from mobileOG-db headers."""
    # Standard format: mobileOG_000123|tnpA|integration/excision|transposase|ISfinder|...
    raw = subject_id + (" " + subject_desc if subject_desc else "")
    parts = subject_id.split("|")
    if len(parts) >= 3:
        return parts[1] or parts[0], parts[2] or "MGE"
    # Infer category from keywords if custom protein FASTA was passed
    if PLASMID_MOG_HINTS.search(raw):
        return parts[0], "transfer_or_replication"
    if IS_MOG_HINTS.search(raw):
        return parts[0], "integration/excision"
    if re.search(r"\b(phage|capsid|tail|terminase|holin|lysin)\b", raw, re.I):
        return parts[0], "phage"
    return parts[0], "MGE"


def load_manifest_meta(path):
    if not path:
        return {}
    meta = {}
    with open(path) as fh:
        hdr = [c.strip() for c in fh.readline().rstrip("\r\n").split("\t")]
        col = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            row = {c: (parts[i].strip() if i < len(parts) else "") for c, i in col.items()}
            sid = row.get("sample")
            if sid:
                meta[sid] = row
    return meta


# ---------------------------------------------------------------------------
# Nucleotide search backend (blastn -> minimap2 -> optional python fallback)
# ---------------------------------------------------------------------------

def run_nuc_search(query_fa, db_path, min_id, min_alen, min_scov, threads, tmp_dir,
                   tag, allow_fallback=False, aligner="auto"):
    """Align query flank sequences against a nucleotide reference FASTA/BLAST DB.

    Returns dict: query_id -> list of dicts {subject, sdesc, pident, alen, scov, qstart0, qend0}
    """
    hits = collections.defaultdict(list)
    if not db_path:
        return hits

    have_blast = bool(shutil.which("blastn") and shutil.which("makeblastdb"))
    have_mm2 = bool(shutil.which("minimap2"))
    is_blast_prefix = os.path.exists(db_path + ".nin") or os.path.exists(db_path + ".nal") or os.path.exists(db_path + ".nsq")

    if aligner == "auto":
        if have_blast and (is_blast_prefix or os.path.isfile(db_path)):
            backend = "blastn"
        elif have_mm2 and os.path.isfile(db_path):
            backend = "minimap2"
        elif allow_fallback and os.path.isfile(db_path):
            backend = "python"
        else:
            sys.exit(
                f"annotate_flanks [{tag}]: requires 'blastn'+'makeblastdb' or 'minimap2' on PATH "
                f"(or pass --allow-python-fallback for small offline tests)."
            )
    else:
        backend = aligner

    if backend == "blastn":
        if is_blast_prefix:
            db_pref = db_path
        else:
            db_pref = os.path.join(tmp_dir, f"db_{tag}")
            subprocess.run(
                ["makeblastdb", "-in", db_path, "-dbtype", "nucl", "-out", db_pref],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
            )
        out_tsv = os.path.join(tmp_dir, f"blast_{tag}.tsv")
        subprocess.run(
            [
                "blastn", "-task", "blastn", "-query", query_fa, "-db", db_pref,
                "-perc_identity", str(min_id), "-evalue", "1e-5",
                "-max_target_seqs", "50", "-num_threads", str(threads),
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send slen stitle",
                "-out", out_tsv,
            ],
            stderr=subprocess.DEVNULL, check=True,
        )
        with open(out_tsv) as fh:
            for line in fh:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                q, s, pid, alen, qs, qe, ss, se, slen = (
                    p[0], p[1], float(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[6]), int(p[7]), int(p[8])
                )
                stitle = p[9] if len(p) > 9 else ""
                scov = 100.0 * (abs(se - ss) + 1) / max(1, slen)
                if pid < min_id or alen < min_alen or scov < min_scov:
                    continue
                hits[q].append({
                    "subject": s, "sdesc": stitle, "pident": pid, "alen": alen, "scov": min(100.0, scov),
                    "qstart0": min(qs, qe) - 1, "qend0": max(qs, qe),
                })

    elif backend == "minimap2":
        paf = os.path.join(tmp_dir, f"mm2_{tag}.paf")
        with open(paf, "w") as out_fh:
            subprocess.run(
                ["minimap2", "-c", "-x", "asm20", "--secondary=yes", "-N", "50", "-t", str(threads), db_path, query_fa],
                stdout=out_fh, stderr=subprocess.DEVNULL, check=True,
            )
        with open(paf) as fh:
            for line in fh:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                q, qs0, qe0, s, slen, ss0, se0, nmatch, blen = (
                    p[0], int(p[2]), int(p[3]), p[5], int(p[6]), int(p[7]), int(p[8]), int(p[9]), int(p[10])
                )
                if blen <= 0:
                    continue
                pid = 100.0 * nmatch / blen
                scov = 100.0 * (se0 - ss0) / max(1, slen)
                if pid < min_id or (qe0 - qs0) < min_alen or scov < min_scov:
                    continue
                hits[q].append({
                    "subject": s, "sdesc": "", "pident": pid, "alen": qe0 - qs0, "scov": min(100.0, scov),
                    "qstart0": qs0, "qend0": qe0,
                })

    elif backend == "python":
        # Lightweight 19-mer seeded matcher for offline smoke tests
        K = 19
        STEP = 5
        refs = {}
        rdesc = {}
        idx = collections.defaultdict(list)
        for sid, desc, sseq in read_fasta(db_path):
            refs[sid] = sseq
            rdesc[sid] = desc
            if len(sseq) < K:
                continue
            for p in range(0, len(sseq) - K + 1, STEP):
                km = sseq[p:p + K]
                if "N" not in km:
                    idx[km].append((sid, p, True))
                    idx[rc(km)].append((sid, p, False))
        for qid, _, qseq in read_fasta(query_fa):
            ql = len(qseq)
            votes = collections.defaultdict(int)
            for p in range(0, ql - K + 1, 4):
                hs = idx.get(qseq[p:p + K])
                if not hs:
                    continue
                for sid, sp, plus in hs:
                    diag = (p - sp) if plus else (p + sp + K)
                    votes[(sid, plus, diag // 12)] += 1
            by_s = collections.defaultdict(list)
            for (sid, plus, dbin), cnt in votes.items():
                if cnt >= 2:
                    by_s[(sid, plus)].append((cnt, dbin * 12))
            for (sid, plus), cands in by_s.items():
                sseq = refs[sid]
                sl = len(sseq)
                cands.sort(reverse=True)
                used = []
                for _, est in cands:
                    if any(abs(est - u) < max(50, sl // 2) for u in used):
                        continue
                    best = None
                    for delta in range(-18, 19):
                        if plus:
                            q0, q1 = est + delta, est + delta + sl
                            o0, o1 = max(0, q0), min(ql, q1)
                            if o1 - o0 < min_alen:
                                continue
                            a = qseq[o0:o1]
                            b = sseq[o0 - q0:o1 - q0]
                        else:
                            q1, q0 = est + delta, est + delta - sl
                            o0, o1 = max(0, q0), min(ql, q1)
                            if o1 - o0 < min_alen:
                                continue
                            a = rc(qseq[o0:o1])
                            b = sseq[q1 - o1:q1 - o0]
                        m = sum(1 for x, y in zip(a, b) if x == y)
                        pid = 100.0 * m / len(a)
                        scov = 100.0 * len(a) / sl
                        if best is None or m > best[0]:
                            best = (m, pid, scov, o0, o1)
                    if best and best[1] >= min_id and best[2] >= min_scov:
                        used.append(est)
                        hits[qid].append({
                            "subject": sid, "sdesc": rdesc[sid], "pident": best[1],
                            "alen": best[4] - best[3], "scov": min(100.0, best[2]),
                            "qstart0": best[3], "qend0": best[4],
                        })

    # Deduplicate overlapping hits per query (keep highest pident*alen)
    dedup = collections.defaultdict(list)
    for q, hs in hits.items():
        hs.sort(key=lambda h: -(h["pident"] * h["alen"]))
        kept = []
        for h in hs:
            ov = any(
                min(h["qend0"], k["qend0"]) - max(h["qstart0"], k["qstart0"])
                >= 0.5 * max(1, min(h["qend0"] - h["qstart0"], k["qend0"] - k["qstart0"]))
                for k in kept
            )
            if not ov:
                kept.append(h)
        dedup[q] = kept
    return dedup


# ---------------------------------------------------------------------------
# ORF calling + Protein MGE search (pyrodigal/prodigal + diamond/mmseqs/blastp)
# ---------------------------------------------------------------------------

def translate_dna(dna):
    aa = []
    for i in range(0, len(dna) - 2, 3):
        c = CODON_TABLE.get(dna[i:i + 3], "X")
        if c == "*":
            break
        aa.append(c)
    return "".join(aa)


def call_orfs(query_fa, out_faa, threads, allow_fallback=False):
    """Predict metagenomic ORFs and return dict: orf_id -> (flank_id, start0, end0, strand)."""
    orf_map = {}

    # 1. Prefer pyrodigal Python package (fast, native multiprocessing, metagenomic mode)
    try:
        import pyrodigal  # type: ignore
        finder = pyrodigal.GeneFinder(meta=True)
        with open(out_faa, "w") as out:
            for fid, _, seq in read_fasta(query_fa):
                for idx, gene in enumerate(finder.find_genes(seq.encode("ascii")), 1):
                    oid = f"{fid}__orf{idx}"
                    s0, e0 = gene.begin - 1, gene.end
                    st = "+" if gene.strand == 1 else "-"
                    orf_map[oid] = (fid, s0, e0, st)
                    out.write(f">{oid} # {gene.begin} # {gene.end} # {gene.strand}\n{gene.translate()}\n")
        return orf_map, "pyrodigal_api"
    except ImportError:
        pass

    # 2. Try pyrodigal or prodigal CLI
    cli = shutil.which("pyrodigal") or shutil.which("prodigal")
    if cli:
        cmd = [cli, "-i", query_fa, "-a", out_faa, "-p", "meta", "-q"]
        if os.path.basename(cli) == "pyrodigal":
            cmd.extend(["-j", str(threads)])
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        for oid, desc, _ in read_fasta(out_faa):
            # Header: >flank_id_1 # start # end # strand # ...
            fid = oid.rsplit("_", 1)[0]
            m = re.match(r"#\s*(\d+)\s*#\s*(\d+)\s*#\s*(-?1)", desc)
            if m:
                s0, e0, st = int(m.group(1)) - 1, int(m.group(2)), ("+" if m.group(3) == "1" else "-")
            else:
                s0, e0, st = 0, 0, "+"
            orf_map[oid] = (fid, s0, e0, st)
        return orf_map, os.path.basename(cli)

    if not allow_fallback:
        sys.exit(
            "annotate_flanks [ORFs]: requires 'pyrodigal' (pip install pyrodigal) or 'prodigal' on PATH "
            "(or pass --allow-python-fallback for small offline tests)."
        )

    # 3. Lightweight 6-frame bacterial ORF finder (>= 70 aa) for offline testing
    min_nt = 210
    with open(out_faa, "w") as out:
        for fid, _, seq in read_fasta(query_fa):
            L = len(seq)
            oid_idx = 0
            for strand, s in (("+", seq), ("-", rc(seq))):
                for frame in range(3):
                    i = frame
                    while i <= L - min_nt:
                        if s[i:i + 3] in START_CODONS:
                            j = i + 3
                            while j <= L - 3 and s[j:j + 3] not in STOP_CODONS:
                                j += 3
                            if j - i >= min_nt:
                                aa = translate_dna(s[i:j])
                                if len(aa) >= 70:
                                    oid_idx += 1
                                    oid = f"{fid}__orf{oid_idx}"
                                    s0, e0 = (i, j) if strand == "+" else (L - j, L - i)
                                    orf_map[oid] = (fid, s0, e0, strand)
                                    out.write(f">{oid}\n{aa}\n")
                                i = j + 3
                                continue
                        i += 3
    return orf_map, "python_orf"


def run_prot_search(query_faa, orf_map, db_path, min_id, min_qcov, threads, tmp_dir, allow_fallback=False):
    """Search predicted ORFs against a protein MGE DB (mobileOG-db) via diamond, mmseqs, or blastp."""
    hits = collections.defaultdict(list)
    if not db_path or not os.path.isfile(query_faa) or os.path.getsize(query_faa) == 0:
        return hits

    have_diamond = bool(shutil.which("diamond"))
    have_mmseqs = bool(shutil.which("mmseqs"))
    have_blastp = bool(shutil.which("blastp") and shutil.which("makeblastdb"))

    if have_diamond:
        dmnd = db_path
        if not db_path.endswith(".dmnd"):
            dmnd = os.path.join(tmp_dir, "mobileog.dmnd")
            subprocess.run(
                ["diamond", "makedb", "--in", db_path, "-d", dmnd, "--quiet"],
                check=True,
            )
        out_tsv = os.path.join(tmp_dir, "diamond_mge.tsv")
        subprocess.run(
            [
                "diamond", "blastp", "-q", query_faa, "-d", dmnd,
                "--id", str(min_id), "--query-cover", str(min_qcov),
                "--evalue", "1e-5", "--max-target-seqs", "1",
                "--threads", str(threads), "--quiet",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "stitle",
                "-o", out_tsv,
            ],
            check=True,
        )
        with open(out_tsv) as fh:
            for line in fh:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                oid, sid, pid = p[0], p[1], float(p[2])
                stitle = p[5] if len(p) > 5 else ""
                if oid in orf_map:
                    fid, s0, e0, _ = orf_map[oid]
                    hits[fid].append({"orf_id": oid, "subject": sid, "sdesc": stitle, "pident": pid, "qstart0": s0, "qend0": e0})
        return hits

    if have_mmseqs and not db_path.endswith(".dmnd"):
        out_tsv = os.path.join(tmp_dir, "mmseqs_mge.tsv")
        mtmp = os.path.join(tmp_dir, "mmseqs_tmp")
        subprocess.run(
            [
                "mmseqs", "easy-search", query_faa, db_path, out_tsv, mtmp,
                "--min-seq-id", str(min_id / 100.0), "-c", str(min_qcov / 100.0), "--cov-mode", "2",
                "-e", "1e-5", "--max-seqs", "1", "--threads", str(threads), "-v", "0",
            ],
            check=True,
        )
        with open(out_tsv) as fh:
            for line in fh:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                oid, sid, pid = p[0], p[1], float(p[2]) * 100.0 if float(p[2]) <= 1.0 else float(p[2])
                if oid in orf_map:
                    fid, s0, e0, _ = orf_map[oid]
                    hits[fid].append({"orf_id": oid, "subject": sid, "sdesc": "", "pident": pid, "qstart0": s0, "qend0": e0})
        return hits

    if have_blastp and not db_path.endswith(".dmnd"):
        db_pref = db_path if os.path.exists(db_path + ".pin") else os.path.join(tmp_dir, "prot_db")
        if db_pref != db_path:
            subprocess.run(["makeblastdb", "-in", db_path, "-dbtype", "prot", "-out", db_pref],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        out_tsv = os.path.join(tmp_dir, "blastp_mge.tsv")
        subprocess.run(
            [
                "blastp", "-query", query_faa, "-db", db_pref, "-evalue", "1e-5",
                "-max_target_seqs", "1", "-num_threads", str(threads),
                "-outfmt", "6 qseqid sseqid pident length qlen stitle", "-out", out_tsv,
            ],
            stderr=subprocess.DEVNULL, check=True,
        )
        with open(out_tsv) as fh:
            for line in fh:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                oid, sid, pid, alen, qlen = p[0], p[1], float(p[2]), int(p[3]), int(p[4])
                if pid >= min_id and 100.0 * alen / max(1, qlen) >= min_qcov and oid in orf_map:
                    fid, s0, e0, _ = orf_map[oid]
                    hits[fid].append({"orf_id": oid, "subject": sid, "sdesc": p[5] if len(p) > 5 else "", "pident": pid, "qstart0": s0, "qend0": e0})
        return hits

    if not allow_fallback or db_path.endswith(".dmnd"):
        sys.exit(
            "annotate_flanks [mobileOG-db]: requires 'diamond' (recommended), 'mmseqs', or 'blastp' on PATH."
        )

    # Lightweight 7-mer AA fallback for small offline tests
    K = 7
    idx = collections.defaultdict(list)
    refs = {}
    rdesc = {}
    for sid, desc, seq in read_fasta(db_path):
        refs[sid] = seq
        rdesc[sid] = desc
        for i in range(0, len(seq) - K + 1, 2):
            idx[seq[i:i + K]].append(sid)
    for oid, _, aaseq in read_fasta(query_faa):
        if len(aaseq) < K or oid not in orf_map:
            continue
        cnt = collections.Counter()
        for i in range(0, len(aaseq) - K + 1, 2):
            for sid in idx.get(aaseq[i:i + K], ()):
                cnt[sid] += 1
        if not cnt:
            continue
        best_s, best_c = cnt.most_common(1)[0]
        frac = best_c / max(1, (len(aaseq) - K) // 2 + 1)
        est_id = 100.0 * (min(1.0, frac) ** (1.0 / K))
        if est_id >= min_id:
            fid, s0, e0, _ = orf_map[oid]
            hits[fid].append({"orf_id": oid, "subject": best_s, "sdesc": rdesc[best_s], "pident": est_id, "qstart0": s0, "qend0": e0})
    return hits


# ---------------------------------------------------------------------------
# geNomad plasmid/virus classifier + MGE-masked host taxonomy + HGT clustering
# ---------------------------------------------------------------------------

def run_genomad(records, genomad_db, min_len, threads, tmp_dir):
    out = {}
    if not genomad_db:
        return out
    if not shutil.which("genomad"):
        sys.exit("annotate_flanks [geNomad]: --genomad-db was set but 'genomad' is not on PATH.")
    eligible = [(fid, seq) for fid, _, seq in records if len(seq) >= min_len]
    if not eligible:
        return out
    in_fa = os.path.join(tmp_dir, "genomad_in.fa")
    gout = os.path.join(tmp_dir, "genomad_out")
    write_fasta(eligible, in_fa)
    subprocess.run(
        ["genomad", "end-to-end", "--cleanup", "--threads", str(threads), in_fa, gout, genomad_db],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
    )
    ps = os.path.join(gout, "genomad_in_summary", "genomad_in_plasmid_summary.tsv")
    vs = os.path.join(gout, "genomad_in_summary", "genomad_in_virus_summary.tsv")
    if os.path.isfile(ps):
        with open(ps) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            col = {c: i for i, c in enumerate(hdr)}
            for line in fh:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                fid = p[col["seq_name"]]
                out[fid] = {
                    "class": "plasmid",
                    "plasmid_score": p[col["plasmid_score"]] if "plasmid_score" in col else "",
                    "virus_score": "",
                    "conj_genes": p[col["conjugation_genes"]] if "conjugation_genes" in col else "",
                }
    if os.path.isfile(vs):
        with open(vs) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            col = {c: i for i, c in enumerate(hdr)}
            for line in fh:
                if not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                fid = p[col["seq_name"]]
                vscore = p[col["virus_score"]] if "virus_score" in col else ""
                if fid not in out or (vscore and float(vscore) > float(out[fid]["plasmid_score"] or 0)):
                    out[fid] = {"class": "virus", "plasmid_score": out.get(fid, {}).get("plasmid_score", ""), "virus_score": vscore, "conj_genes": ""}
    return out


def mask_mge_intervals(records, mask_intervals_by_id, min_unmasked_bp=300):
    """N-mask IS/MGE intervals before host taxonomy so promiscuous transposases don't hijack host calls."""
    masked = []
    unmasked_bp = {}
    for fid, _, seq in records:
        ivs = mask_intervals_by_id.get(fid)
        if not ivs:
            unmasked_bp[fid] = len(seq)
            if len(seq) >= min_unmasked_bp:
                masked.append((fid, seq))
            continue
        arr = list(seq)
        for s0, e0 in ivs:
            a, b = max(0, s0), min(len(arr), e0)
            for i in range(a, b):
                arr[i] = "N"
        mseq = "".join(arr)
        rem = len(mseq) - mseq.count("N")
        unmasked_bp[fid] = rem
        if rem >= min_unmasked_bp:
            masked.append((fid, mseq))
    return masked, unmasked_bp


def run_host_taxonomy(masked_records, kraken2_db, mmseqs_tax_db, threads, tmp_dir):
    out = {}
    if not masked_records or not (kraken2_db or mmseqs_tax_db):
        return out
    in_fa = os.path.join(tmp_dir, "masked_host.fa")
    write_fasta(masked_records, in_fa)

    if kraken2_db:
        if not shutil.which("kraken2"):
            sys.exit("annotate_flanks [Kraken2]: --kraken2-db was set but 'kraken2' is not on PATH.")
        kout = os.path.join(tmp_dir, "kraken2.out")
        subprocess.run(
            ["kraken2", "--db", kraken2_db, "--threads", str(threads), "--use-names", "--output", kout, in_fa],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
        )
        with open(kout) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3 and p[0] == "C":
                    # With --use-names, p[2] looks like "Escherichia coli (taxid 562)"
                    m = re.match(r"^(.*)\s+\(taxid\s+(\d+)\)$", p[2].strip())
                    if m:
                        out[p[1]] = (m.group(2), m.group(1))
                    else:
                        out[p[1]] = ("", p[2].strip())
        return out

    if mmseqs_tax_db:
        if not shutil.which("mmseqs"):
            sys.exit("annotate_flanks [MMseqs2 tax]: --mmseqs-tax-db was set but 'mmseqs' is not on PATH.")
        pref = os.path.join(tmp_dir, "mmtax")
        mtmp = os.path.join(tmp_dir, "mmtax_tmp")
        subprocess.run(
            ["mmseqs", "easy-taxonomy", in_fa, mmseqs_tax_db, pref, mtmp, "--threads", str(threads), "-v", "0"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
        )
        lca = f"{pref}_lca.tsv"
        if os.path.isfile(lca):
            with open(lca) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 4 and p[1] != "0":
                        out[p[0]] = (p[1], p[3])
    return out


def run_hgt_clustering(records, exclude_ids, min_len, min_id, min_cov, threads, tmp_dir, allow_fallback=False):
    """Cluster non-vector flanks >= min_len at >= min_id identity and >= min_cov coverage."""
    eligible = [(fid, seq) for fid, _, seq in records if fid not in exclude_ids and len(seq) >= min_len]
    if not eligible:
        return {}
    in_fa = os.path.join(tmp_dir, "hgt_in.fa")
    write_fasta(eligible, in_fa)

    if shutil.which("mmseqs"):
        pref = os.path.join(tmp_dir, "hgt_cl")
        mtmp = os.path.join(tmp_dir, "hgt_tmp")
        subprocess.run(
            [
                "mmseqs", "easy-linclust", in_fa, pref, mtmp,
                "--min-seq-id", str(min_id), "-c", str(min_cov), "--cov-mode", "0",
                "--threads", str(threads), "-v", "0",
            ],
            check=True,
        )
        rep_to_cid = {}
        fid_to_cid = {}
        with open(f"{pref}_cluster.tsv") as fh:
            for line in fh:
                rep, mem = line.rstrip("\n").split("\t")[:2]
                if rep not in rep_to_cid:
                    rep_to_cid[rep] = f"HGT_{len(rep_to_cid) + 1:05d}"
                fid_to_cid[mem] = rep_to_cid[rep]
        return fid_to_cid

    if shutil.which("vsearch"):
        uc = os.path.join(tmp_dir, "hgt.uc")
        subprocess.run(
            [
                "vsearch", "--cluster_fast", in_fa, "--id", str(min_id),
                "--query_cov", str(min_cov), "--target_cov", str(min_cov),
                "--strand", "both", "--threads", str(threads), "--uc", uc, "--quiet",
            ],
            check=True,
        )
        fid_to_cid = {}
        with open(uc) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if p[0] in ("S", "H"):
                    fid_to_cid[p[8]] = f"HGT_{int(p[1]) + 1:05d}"
        return fid_to_cid

    if not allow_fallback:
        sys.exit(
            "annotate_flanks [--cluster-hgt]: requires 'mmseqs' or 'vsearch' on PATH "
            "(or pass --allow-python-fallback for small offline tests)."
        )

    # Strand-invariant 21-mer Jaccard greedy clustering (offline smoke test fallback)
    K = 21
    kmer_thr = min_id ** K
    eligible.sort(key=lambda x: -len(x[1]))
    reps = []  # (cid, length, canonical_kmer_set)
    fid_to_cid = {}
    for fid, seq in eligible:
        rc_s = rc(seq)
        kms = set()
        for i in range(0, len(seq) - K + 1, 5):
            a = seq[i:i + K]
            b = rc_s[len(seq) - i - K:len(seq) - i]
            kms.add(a if a < b else b)
        placed = False
        for cid, rlen, rkms in reps:
            if min(len(seq), rlen) / max(len(seq), rlen) < min_cov:
                continue
            ov = len(kms & rkms) / max(1, min(len(kms), len(rkms)))
            if ov >= kmer_thr * 0.85:
                fid_to_cid[fid] = cid
                placed = True
                break
        if not placed:
            cid = f"HGT_{len(reps) + 1:05d}"
            reps.append((cid, len(seq), kms))
            fid_to_cid[fid] = cid
    return fid_to_cid


# ---------------------------------------------------------------------------
# Per-arm orchestration + Head-to-head ARGjacent vs Assembly comparison
# ---------------------------------------------------------------------------

def annotate_arm(arm_name, fasta_path, manifest_meta, args, out_dir):
    sys.stderr.write(f"\n=== Annotating arm '{arm_name}' from {fasta_path} ===\n")
    raw_recs = read_fasta(fasta_path)
    records = [(fid, desc, seq) for fid, desc, seq in raw_recs if len(seq) >= args.min_len]
    sys.stderr.write(f"  Loaded {len(records):,} / {len(raw_recs):,} flanks >= {args.min_len} bp.\n")

    if not records:
        return [], []

    parsed = {fid: parse_flank_header(fid, len(seq), manifest_meta) for fid, _, seq in records}

    with tempfile.TemporaryDirectory(prefix=f"argj_ann_{arm_name}_") as tmp_dir:
        filt_fa = os.path.join(tmp_dir, "flanks.fa")
        write_fasta([(fid, seq) for fid, _, seq in records], filt_fa)

        # Layer 1: Lab reagent / cloning vector screen (FIRST STEP!)
        vec_hits = {}
        if args.univec_db:
            sys.stderr.write("  [1/5] Screening cloning/reagent vectors (--univec-db)...\n")
            vec_hits = run_nuc_search(
                filt_fa, args.univec_db, min_id=90.0, min_alen=150, min_scov=0.0,
                threads=args.threads, tmp_dir=tmp_dir, tag="univec",
                allow_fallback=args.allow_python_fallback, aligner=args.nuc_aligner,
            )

        # Layer 2a: Co-localized ARGs (nucleotide)
        coarg_hits = {}
        if args.arg_db:
            sys.stderr.write("  [2/5] Annotating co-localized ARGs (--arg-db)...\n")
            coarg_hits = run_nuc_search(
                filt_fa, args.arg_db, min_id=90.0, min_alen=200, min_scov=60.0,
                threads=args.threads, tmp_dir=tmp_dir, tag="coarg",
                allow_fallback=args.allow_python_fallback, aligner=args.nuc_aligner,
            )

        # Layer 2b: IS elements (ISfinder nucleotide)
        is_hits = {}
        if args.isfinder_db:
            sys.stderr.write("  [2/5] Annotating IS elements (--isfinder-db)...\n")
            is_hits = run_nuc_search(
                filt_fa, args.isfinder_db, min_id=80.0, min_alen=200, min_scov=0.0,
                threads=args.threads, tmp_dir=tmp_dir, tag="isfinder",
                allow_fallback=args.allow_python_fallback, aligner=args.nuc_aligner,
            )

        # Layer 2c: ORF calling + mobileOG-db protein MGE search
        orf_counts = collections.Counter()
        mog_hits = {}
        if args.mobileog_db or args.call_orfs:
            sys.stderr.write("  [2/5] Predicting ORFs & searching mobileOG-db...\n")
            faa = os.path.join(tmp_dir, "orfs.faa")
            orf_map, caller = call_orfs(filt_fa, faa, args.threads, allow_fallback=args.allow_python_fallback)
            for oid, (fid, _, _, _) in orf_map.items():
                orf_counts[fid] += 1
            sys.stderr.write(f"        Called {len(orf_map):,} ORFs via {caller}.\n")
            if args.mobileog_db:
                mog_hits = run_prot_search(
                    faa, orf_map, args.mobileog_db, min_id=60.0, min_qcov=60.0,
                    threads=args.threads, tmp_dir=tmp_dir, allow_fallback=args.allow_python_fallback,
                )

        # Layer 3: Plasmid / viral classification (geNomad + PlasmidFinder)
        gn_hits = {}
        if args.genomad_db:
            sys.stderr.write(f"  [3/5] Running geNomad plasmid/virus classification (>= {args.genomad_min_len} bp)...\n")
            gn_hits = run_genomad(records, args.genomad_db, args.genomad_min_len, args.threads, tmp_dir)

        pf_hits = {}
        if args.plasmidfinder_db:
            sys.stderr.write("  [3/5] Screening PlasmidFinder replicons (--plasmidfinder-db)...\n")
            pf_hits = run_nuc_search(
                filt_fa, args.plasmidfinder_db, min_id=80.0, min_alen=150, min_scov=60.0,
                threads=args.threads, tmp_dir=tmp_dir, tag="plasmidfinder",
                allow_fallback=args.allow_python_fallback, aligner=args.nuc_aligner,
            )

        # Layer 4: MGE-masked host taxonomy
        mask_ivs = collections.defaultdict(list)
        for fid, hs in is_hits.items():
            for h in hs:
                mask_ivs[fid].append((h["qstart0"], h["qend0"]))
        for fid, hs in mog_hits.items():
            for h in hs:
                mask_ivs[fid].append((h["qstart0"], h["qend0"]))
        masked_recs, unmasked_bp = mask_mge_intervals(records, mask_ivs, min_unmasked_bp=300)
        host_tax = {}
        if args.kraken2_db or args.mmseqs_tax_db:
            sys.stderr.write(f"  [4/5] Running MGE-masked host taxonomy on {len(masked_recs):,} flanks...\n")
            host_tax = run_host_taxonomy(masked_recs, args.kraken2_db, args.mmseqs_tax_db, args.threads, tmp_dir)

        # Layer 5: Cross-biome HGT clustering (strictly excluding reagent vectors!)
        vec_ids = set(vec_hits.keys())
        fid_to_hgt = {}
        if args.cluster_hgt:
            sys.stderr.write(f"  [5/5] Clustering non-vector flanks for cross-biome HGT (>={args.hgt_id*100:.0f}% id, >={args.hgt_cov*100:.0f}% cov)...\n")
            fid_to_hgt = run_hgt_clustering(
                records, vec_ids, min_len=args.hgt_min_len, min_id=args.hgt_id, min_cov=args.hgt_cov,
                threads=args.threads, tmp_dir=tmp_dir, allow_fallback=args.allow_python_fallback,
            )

    # Pre-aggregate HGT cluster metadata (biomes, countries, studies, focal ARGs)
    hgt_members = collections.defaultdict(list)
    for fid, cid in fid_to_hgt.items():
        hgt_members[cid].append(fid)

    hgt_summary_rows = []
    hgt_info = {}
    for cid, fids in sorted(hgt_members.items()):
        biomes = sorted({parsed[f]["biome_group"] for f in fids if parsed[f]["biome_group"]})
        countries = sorted({parsed[f]["country"] for f in fids if parsed[f]["country"]})
        studies = sorted({parsed[f]["study_accession"] for f in fids if parsed[f]["study_accession"]})
        genes = sorted({parsed[f]["gene"] for f in fids if parsed[f]["gene"]})
        info = {
            "cluster_id": cid,
            "n_flanks": len(fids),
            "n_samples": len({parsed[f]["sample"] for f in fids}),
            "n_studies": len(studies),
            "n_biomes": len(biomes),
            "n_countries": len(countries),
            "is_cross_biome": 1 if len(biomes) >= 2 else 0,
            "biomes": ",".join(biomes),
            "countries": ",".join(countries[:10]),
            "genes": ",".join(genes[:10]),
        }
        hgt_info[cid] = info
        hgt_summary_rows.append(info)

    # Assemble per-flank unified annotation table
    ann_rows = []
    for fid, _, _ in records:
        p = parsed[fid]
        side, flen = p["side"], p["flank_len"]

        vh = vec_hits.get(fid, [])
        is_vec = 1 if vh else 0
        vec_names = ",".join(sorted({h["subject"] for h in vh})[:5])
        vec_cov = sum(h["alen"] for h in vh)
        vec_dist = min((distance_to_arg(side, flen, h["qstart0"], h["qend0"]) for h in vh), default="")

        ch = coarg_hits.get(fid, [])
        coarg_names = ",".join(sorted({h["subject"] for h in ch})[:8])
        coarg_dist = min((distance_to_arg(side, flen, h["qstart0"], h["qend0"]) for h in ch), default="")

        ih = is_hits.get(fid, [])
        is_parsed = [parse_is_family(h["subject"], h["sdesc"]) for h in ih]
        is_names = ",".join(sorted({x[0] for x in is_parsed})[:6])
        is_fams = ",".join(sorted({x[1] for x in is_parsed})[:6])
        is_dist = min((distance_to_arg(side, flen, h["qstart0"], h["qend0"]) for h in ih), default="")

        mh = mog_hits.get(fid, [])
        m_parsed = [parse_mobileog_target(h["subject"], h["sdesc"]) for h in mh]
        mog_genes = ",".join(sorted({x[0] for x in m_parsed})[:6])
        mog_cats = ",".join(sorted({x[1] for x in m_parsed})[:6])
        mog_dist = min((distance_to_arg(side, flen, h["qstart0"], h["qend0"]) for h in mh), default="")

        gn = gn_hits.get(fid, {})
        gn_cls = gn.get("class", "not_run" if not args.genomad_db else ("too_short" if flen < args.genomad_min_len else "chromosome_or_unclassified"))
        pf = pf_hits.get(fid, [])
        pf_reps = ",".join(sorted({h["subject"] for h in pf})[:5])

        # Hierarchical element classification (reagent_vector takes strict precedence!)
        has_plasmid_mog = any(PLASMID_MOG_HINTS.search(x[0] + " " + x[1]) for x in m_parsed)
        if is_vec:
            elem_cls = "reagent_vector"
        elif gn_cls == "plasmid" or pf or has_plasmid_mog:
            elem_cls = "plasmid"
        elif gn_cls == "virus":
            elem_cls = "virus"
        elif ih or mh:
            elem_cls = "mobile_unclassified"
        else:
            elem_cls = "chromosomal_or_other"

        htax_id, htax_name = host_tax.get(fid, ("", ""))
        cid = fid_to_hgt.get(fid, "")
        hi = hgt_info.get(cid, {})

        ann_rows.append({
            "arm": arm_name,
            "flank_id": fid,
            "sample": p["sample"],
            "study_accession": p["study_accession"],
            "biome_group": p["biome_group"],
            "country": p["country"],
            "assembler": p["assembler"],
            "gene": p["gene"],
            "side": side,
            "ctx_id": p["ctx_id"],
            "flank_len": flen,
            "element_class": elem_cls,
            "is_reagent_vector": is_vec,
            "vector_hits": vec_names,
            "vector_cov_bp": vec_cov,
            "vector_nearest_bp": vec_dist,
            "n_coargs": len(ch),
            "coarg_genes": coarg_names,
            "nearest_coarg_bp": coarg_dist,
            "has_is": 1 if ih else 0,
            "n_is": len(ih),
            "is_elements": is_names,
            "is_families": is_fams,
            "nearest_is_bp": is_dist,
            "n_orfs": orf_counts.get(fid, 0),
            "n_mge_orfs": len(mh),
            "mge_categories": mog_cats,
            "mge_genes": mog_genes,
            "nearest_mge_bp": mog_dist,
            "genomad_class": gn_cls,
            "genomad_plasmid_score": gn.get("plasmid_score", ""),
            "plasmid_replicons": pf_reps,
            "unmasked_bp": unmasked_bp.get(fid, flen),
            "host_taxid": htax_id,
            "host_taxon": htax_name,
            "hgt_cluster_id": cid,
            "hgt_cluster_size": hi.get("n_flanks", ""),
            "hgt_cluster_n_biomes": hi.get("n_biomes", ""),
            "is_cross_biome_hgt": hi.get("is_cross_biome", 0) if cid else 0,
        })

    # Write per-arm TSVs
    ann_path = os.path.join(out_dir, f"{arm_name}_annotations.tsv")
    cols = list(ann_rows[0].keys())
    with open(ann_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in ann_rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    if hgt_summary_rows:
        hgt_path = os.path.join(out_dir, f"{arm_name}_hgt_clusters.tsv")
        hcols = list(hgt_summary_rows[0].keys())
        with open(hgt_path, "w") as out:
            out.write("\t".join(hcols) + "\n")
            for r in sorted(hgt_summary_rows, key=lambda x: (-x["is_cross_biome"], -x["n_samples"], x["cluster_id"])):
                out.write("\t".join(str(r[c]) for c in hcols) + "\n")

    return ann_rows, hgt_summary_rows


def summarize_arm(label, rows, hgt_rows):
    n = len(rows)
    if n == 0:
        return {
            "group": label, "flanks": 0, "samples": 0, "reagent_vectors": 0, "reagent_vector_pct": "0.00",
            "has_is": 0, "median_nearest_is_bp": "", "has_coarg": 0, "has_mge_orf": 0,
            "plasmid_calls": 0, "host_assigned": 0, "cross_biome_hgt_flanks": 0, "cross_biome_hgt_clusters": 0,
        }
    is_dists = [int(r["nearest_is_bp"]) for r in rows if r["nearest_is_bp"] != ""]
    n_vec = sum(r["is_reagent_vector"] for r in rows)
    return {
        "group": label,
        "flanks": n,
        "samples": len({r["sample"] for r in rows}),
        "reagent_vectors": n_vec,
        "reagent_vector_pct": f"{100.0 * n_vec / n:.2f}",
        "has_is": sum(r["has_is"] for r in rows),
        "median_nearest_is_bp": int(statistics.median(is_dists)) if is_dists else "",
        "has_coarg": sum(1 for r in rows if r["n_coargs"] > 0),
        "has_mge_orf": sum(1 for r in rows if r["n_mge_orfs"] > 0),
        "plasmid_calls": sum(1 for r in rows if r["element_class"] == "plasmid"),
        "host_assigned": sum(1 for r in rows if r["host_taxon"]),
        "cross_biome_hgt_flanks": sum(r["is_cross_biome_hgt"] for r in rows),
        "cross_biome_hgt_clusters": sum(1 for h in hgt_rows if h["is_cross_biome"] == 1),
    }


def write_vector_audit(arm_rows_map, out_path):
    """Per-BioProject audit of lab reagent/cloning vector contamination across arms."""
    rows = []
    for arm, arows in arm_rows_map.items():
        by_study = collections.defaultdict(list)
        for r in arows:
            by_study[r["study_accession"] or "UNKNOWN"].append(r)
        for st, items in sorted(by_study.items()):
            vecs = [x for x in items if x["is_reagent_vector"] == 1]
            top_g = collections.Counter(x["gene"] for x in vecs).most_common(5)
            top_v = collections.Counter(x["vector_hits"] for x in vecs if x["vector_hits"]).most_common(3)
            rows.append({
                "arm": arm,
                "study_accession": st,
                "biomes": ",".join(sorted({x["biome_group"] for x in items if x["biome_group"]})),
                "total_flanks": len(items),
                "vector_flanks": len(vecs),
                "vector_pct": f"{100.0 * len(vecs) / max(1, len(items)):.2f}",
                "affected_samples": len({x["sample"] for x in vecs}),
                "total_samples": len({x["sample"] for x in items}),
                "top_vector_args": ",".join(f"{k}({v})" for k, v in top_g),
                "top_vectors": ",".join(f"{k}({v})" for k, v in top_v),
            })
    cols = [
        "arm", "study_accession", "biomes", "total_flanks", "vector_flanks", "vector_pct",
        "affected_samples", "total_samples", "top_vector_args", "top_vectors",
    ]
    with open(out_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: (-x["vector_flanks"], x["arm"], x["study_accession"])):
            out.write("\t".join(str(r[c]) for c in cols) + "\n")


def compare_arms_per_side(argj_rows, asm_rows, out_path):
    """Side-by-side annotation comparison keyed by (sample, gene, side)."""
    by_a = collections.defaultdict(list)
    by_m = collections.defaultdict(list)
    for r in argj_rows:
        by_a[(r["sample"], r["gene"], r["side"])].append(r)
    for r in asm_rows:
        by_m[(r["sample"], r["gene"], r["side"])].append(r)

    keys = sorted(set(by_a) | set(by_m))
    out_rows = []
    for k in keys:
        sample, gene, side = k
        aa = by_a.get(k, [])
        mm = by_m.get(k, [])
        ref = aa[0] if aa else mm[0]

        def agg(items):
            if not items:
                return {"n": 0, "cls": "absent", "vec": 0, "is": 0, "is_d": "", "coarg": 0, "host": "", "hgt": 0}
            clss = sorted({x["element_class"] for x in items})
            is_ds = [int(x["nearest_is_bp"]) for x in items if x["nearest_is_bp"] != ""]
            hosts = sorted({x["host_taxon"] for x in items if x["host_taxon"]})
            return {
                "n": len(items),
                "cls": "+".join(clss),
                "vec": 1 if any(x["is_reagent_vector"] for x in items) else 0,
                "is": 1 if any(x["has_is"] for x in items) else 0,
                "is_d": min(is_ds) if is_ds else "",
                "coarg": max(x["n_coargs"] for x in items),
                "host": ";".join(hosts[:3]),
                "hgt": 1 if any(x["is_cross_biome_hgt"] for x in items) else 0,
            }

        sa, sm = agg(aa), agg(mm)
        flags = []
        if sa["vec"] and not sm["vec"]:
            flags.append("vector_argj_only")
        elif sm["vec"] and not sa["vec"]:
            flags.append("vector_asm_only")
        if sa["is"] and not sm["is"]:
            flags.append("is_argj_only")
        elif sm["is"] and not sa["is"]:
            flags.append("is_asm_only")
        if sa["n"] >= 2 and sm["n"] <= 1 and "+" in sa["cls"]:
            flags.append("multi_class_argj_collapsed_in_asm")
        if sa["cls"] != "absent" and sm["cls"] != "absent" and sa["cls"] != sm["cls"]:
            flags.append("class_diff")
        if sa["host"] and sm["host"] and sa["host"] != sm["host"]:
            flags.append("host_diff")

        out_rows.append({
            "sample": sample,
            "study_accession": ref["study_accession"],
            "biome_group": ref["biome_group"],
            "assembler": ref["assembler"],
            "gene": gene,
            "side": side,
            "argj_flanks": sa["n"],
            "argj_class": sa["cls"],
            "argj_vector": sa["vec"],
            "argj_has_is": sa["is"],
            "argj_nearest_is_bp": sa["is_d"],
            "argj_max_coargs": sa["coarg"],
            "argj_host": sa["host"],
            "argj_cross_biome_hgt": sa["hgt"],
            "asm_flanks": sm["n"],
            "asm_class": sm["cls"],
            "asm_vector": sm["vec"],
            "asm_has_is": sm["is"],
            "asm_nearest_is_bp": sm["is_d"],
            "asm_max_coargs": sm["coarg"],
            "asm_host": sm["host"],
            "asm_cross_biome_hgt": sm["hgt"],
            "delta_flags": ",".join(flags) if flags else "concordant",
        })

    cols = list(out_rows[0].keys()) if out_rows else ["sample", "gene", "side", "delta_flags"]
    with open(out_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in out_rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    return out_rows


def main():
    ap = argparse.ArgumentParser(
        description="Unified 5-layer annotation and biological comparison for ARGjacent certified flanks vs. pre-made assemblies."
    )
    inp = ap.add_argument_group("Inputs (provide --argj-fa/--asm-fa pair, OR standalone --fasta)")
    inp.add_argument("--argj-fa", help="ARGjacent certified flanks FASTA (from compare_argjacent_vs_assembly.py --export-fasta)")
    inp.add_argument("--asm-fa", help="Pre-made assembly flanks FASTA (from compare_argjacent_vs_assembly.py --export-fasta)")
    inp.add_argument("--fasta", help="Single FASTA to annotate (standalone mode)")
    inp.add_argument("--arm-name", default="argjacent", help="Arm label for --fasta standalone mode (default: argjacent)")
    inp.add_argument("--manifest", help="Optional global manifest TSV (attaches biome, country, BioProject, assembler)")
    inp.add_argument("-o", "--out-dir", required=True, help="Output directory for annotation and comparison tables")

    dbs = ap.add_argument_group("Annotation databases (enable any subset)")
    dbs.add_argument("--univec-db", help="Layer 1: UniVec / cloning-vector nucleotide FASTA or BLAST DB prefix")
    dbs.add_argument("--arg-db", help="Layer 2a: ResFinder / CARD nucleotide FASTA or BLAST DB prefix (co-localized ARGs)")
    dbs.add_argument("--isfinder-db", help="Layer 2b: ISfinder nucleotide FASTA or BLAST DB prefix")
    dbs.add_argument("--mobileog-db", help="Layer 2c: mobileOG-db protein FASTA or DIAMOND .dmnd (triggers ORF calling)")
    dbs.add_argument("--call-orfs", action="store_true", help="Call ORFs even if --mobileog-db is not provided")
    dbs.add_argument("--genomad-db", help="Layer 3a: geNomad database directory")
    dbs.add_argument("--genomad-min-len", type=int, default=2000, help="Min flank length for geNomad (default: 2000 bp)")
    dbs.add_argument("--plasmidfinder-db", help="Layer 3b: PlasmidFinder replicon nucleotide FASTA or BLAST DB prefix")
    dbs.add_argument("--kraken2-db", help="Layer 4a: Kraken2 database directory (run on IS/MGE-masked sequences)")
    dbs.add_argument("--mmseqs-tax-db", help="Layer 4b: MMseqs2 taxonomy database prefix (alternative to Kraken2)")
    dbs.add_argument("--cluster-hgt", action="store_true", help="Layer 5: cluster non-vector flanks for cross-biome HGT")
    dbs.add_argument("--hgt-min-len", type=int, default=1000, help="Min flank length for HGT clustering (default: 1000 bp)")
    dbs.add_argument("--hgt-id", type=float, default=0.99, help="Min sequence identity for HGT clustering (default: 0.99)")
    dbs.add_argument("--hgt-cov", type=float, default=0.90, help="Min aligned coverage for HGT clustering (default: 0.90)")

    opt = ap.add_argument_group("Runtime & backend options")
    opt.add_argument("--min-len", type=int, default=1000, help="Minimum flank length to annotate (default: 1000 bp)")
    opt.add_argument("--nuc-aligner", choices=("auto", "blastn", "minimap2", "python"), default="auto",
                     help="Nucleotide search backend (default: auto -> blastn > minimap2)")
    opt.add_argument("--allow-python-fallback", action="store_true",
                     help="Allow pure-Python fallbacks for offline smoke tests (disabled by default on HPC so missing tools fail loudly)")
    opt.add_argument("-t", "--threads", type=int, default=8, help="Worker threads (default: 8)")
    args = ap.parse_args()

    if not (args.argj_fa or args.asm_fa or args.fasta):
        ap.error("Provide --argj-fa and/or --asm-fa, or standalone --fasta.")

    os.makedirs(args.out_dir, exist_ok=True)
    manifest_meta = load_manifest_meta(args.manifest)

    arms = []
    if args.fasta:
        arms.append((args.arm_name, args.fasta))
    else:
        if args.argj_fa:
            arms.append(("argjacent", args.argj_fa))
        if args.asm_fa:
            arms.append(("assembly", args.asm_fa))

    arm_rows = {}
    arm_hgt = {}
    summary_rows = []
    for name, path in arms:
        rows, hgts = annotate_arm(name, path, manifest_meta, args, args.out_dir)
        arm_rows[name] = rows
        arm_hgt[name] = hgts
        summary_rows.append(summarize_arm(f"{name}:ALL", rows, hgts))
        by_b = collections.defaultdict(list)
        for r in rows:
            if r["biome_group"]:
                by_b[r["biome_group"]].append(r)
        for b in sorted(by_b):
            summary_rows.append(summarize_arm(f"{name}:biome:{b}", by_b[b], []))

    # Write reagent vector audit
    if args.univec_db:
        write_vector_audit(arm_rows, os.path.join(args.out_dir, "reagent_vector_audit.tsv"))

    # Write overall + biome-stratified summary
    sum_path = os.path.join(args.out_dir, "annotation_comparison_summary.tsv")
    scols = list(summary_rows[0].keys())
    with open(sum_path, "w") as out:
        out.write("\t".join(scols) + "\n")
        for r in summary_rows:
            out.write("\t".join(str(r[c]) for c in scols) + "\n")

    # Head-to-head per-side annotation comparison when both arms ran
    delta_counts = collections.Counter()
    if "argjacent" in arm_rows and "assembly" in arm_rows:
        per_side = compare_arms_per_side(
            arm_rows["argjacent"], arm_rows["assembly"],
            os.path.join(args.out_dir, "annotation_comparison_per_side.tsv"),
        )
        for r in per_side:
            for fl in r["delta_flags"].split(","):
                delta_counts[fl] += 1

    print(f"\n=== Unified Flank Annotation Summary ({args.out_dir}) ===")
    for r in summary_rows:
        if not r["group"].endswith(":ALL"):
            continue
        print(
            f"  [{r['group']:<14s}] flanks={r['flanks']:,} | vectors={r['reagent_vectors']:,} ({r['reagent_vector_pct']}%) | "
            f"IS={r['has_is']:,} (med dist={r['median_nearest_is_bp'] or 'n/a'} bp) | co-ARG={r['has_coarg']:,} | "
            f"MGE_ORF={r['has_mge_orf']:,} | plasmid={r['plasmid_calls']:,} | "
            f"cross-biome HGT={r['cross_biome_hgt_flanks']:,} ({r['cross_biome_hgt_clusters']:,} clusters)"
        )
    if delta_counts:
        print("  Per-side delta flags (ARGjacent vs Assembly): " +
              ", ".join(f"{k}={v}" for k, v in delta_counts.most_common()))


if __name__ == "__main__":
    main()
