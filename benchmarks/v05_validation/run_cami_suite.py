#!/usr/bin/env python3
"""Storage-safe streaming runner for frozen ARGjacent v0.5.0 on CAMI II and CAMI III suites.

Storage guarantees (designed for <= 17 GiB free disk):
  - Never saves multi-GB read tarballs to disk: streams HTTP -> tar -> gzip -dc -> deint -> R1.fq.gz + R2.fq.gz.
  - Checks free disk before every sample (aborts cleanly if below `--min-free-gb`, default 7.5 GiB).
  - Deletes R1.fq.gz and R2.fq.gz IMMEDIATELY after KMA + ARGjacent finish, before scoring starts.
  - Deletes `run.tmp/` IMMEDIATELY after `score_v05_cami.py` finishes scoring `W=1000`.
  - Runs `argjacent-cohort` in pure JSON-only mode (`depth_bins` + `gene_depth` in `flanks.json`),
    then deletes the suite's temporary genomes/BLAST DB before moving to the next suite.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESFINDER_FA = "/private/tmp/claude-501/-Users-amulya-Documents-bigargs-anchorwalk-0-1/04d0a03e-7cdb-49b9-80f1-a80edf158510/scratchpad/resfinder_db/all.fsa"
KMA = "/opt/miniconda3/bin/kma"
KMA_INDEX = "/opt/miniconda3/bin/kma_index"
MINIMAP2 = "/opt/miniconda3/pkgs/minimap2-2.31-h6bd33b9_0/bin/minimap2"
BLASTN = "/opt/miniconda3/pkgs/blast-2.17.0-h81fc81c_1/bin/blastn"
MAKEBLASTDB = "/opt/miniconda3/pkgs/blast-2.17.0-h81fc81c_1/bin/makeblastdb"
ARGJACENT = os.path.join(REPO, "bin", "argjacent")
KMA_GENES = os.path.join(REPO, "benchmarks", "cami", "kma_genes.py")
SCORE_SAMPLE = os.path.join(REPO, "benchmarks", "v05_validation", "score_v05_cami.py")
SCORE_COHORT = os.path.join(REPO, "benchmarks", "v05_validation", "score_v05_cohort.py")
DEINT_C = os.path.join(REPO, "benchmarks", "cami", "deint.c")
PREP_CAMI3 = os.path.join(REPO, "benchmarks", "cami", "prep_genomes.py")
PSTREAM = os.path.join(REPO, "benchmarks", "cami", "pstream.py")

WORK_ROOT = "/tmp/argj_val/work"
RESULTS_ROOT = os.path.join(REPO, "benchmarks", "v05_validation", "results")


def free_gb(path="/tmp"):
    return shutil.disk_usage(path).free / (1024 ** 3)


def check_disk(min_gb, stage=""):
    fg = free_gb("/tmp")
    print(f"[disk] {stage}: {fg:.2f} GiB free", flush=True)
    if fg < min_gb:
        raise RuntimeError(f"Low disk ({fg:.2f} GiB free < {min_gb:.1f} GiB required) at stage: {stage}")


def ensure_tools():
    os.makedirs(WORK_ROOT, exist_ok=True)
    deint = os.path.join(WORK_ROOT, "deint")
    if not os.path.isfile(deint):
        subprocess.run(["cc", "-O2", "-o", deint, DEINT_C], check=True)
    kmadb = os.path.join(WORK_ROOT, "kmadb")
    if not os.path.isfile(kmadb + ".name"):
        subprocess.run([KMA_INDEX, "-i", RESFINDER_FA, "-o", kmadb], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return deint, kmadb


def stream_deinterleave(url, r1_out, r2_out, deint_bin):
    """Stream remote .tar.gz containing anonymous_reads.fq.gz directly into R1.fq.gz and R2.fq.gz (0 tar on disk)."""
    t0 = time.time()
    tmp_dir = os.path.dirname(r1_out)
    f1 = os.path.join(tmp_dir, "fifo_r1")
    f2 = os.path.join(tmp_dir, "fifo_r2")
    for p in (f1, f2, r1_out, r2_out):
        if os.path.exists(p):
            os.remove(p)
    os.mkfifo(f1)
    os.mkfifo(f2)

    p_gz1 = subprocess.Popen(f"gzip -1 < '{f1}' > '{r1_out}'", shell=True)
    p_gz2 = subprocess.Popen(f"gzip -1 < '{f2}' > '{r2_out}'", shell=True)
    p_deint = subprocess.Popen([deint_bin, f1, f2], stdin=subprocess.PIPE)
    p_gunzip = subprocess.Popen(["gzip", "-dc"], stdin=subprocess.PIPE, stdout=p_deint.stdin)
    p_deint.stdin.close()

    # Stream tar.gz over HTTP via curl (with retries), extract only anonymous_reads.fq.gz into gunzip
    p_curl = subprocess.Popen(["curl", "-sSL", "--retry", "5", "--retry-delay", "5", url], stdout=subprocess.PIPE)
    found = False
    try:
        with tarfile.open(fileobj=p_curl.stdout, mode="r|gz") as tf:
            for m in tf:
                if m.isfile() and m.name.endswith("anonymous_reads.fq.gz"):
                    found = True
                    src = tf.extractfile(m)
                    while True:
                        buf = src.read(1 << 20)
                        if not buf:
                            break
                        p_gunzip.stdin.write(buf)
                    break
    finally:
        try:
            p_gunzip.stdin.close()
        except Exception:
            pass
        p_curl.stdout.close()
        p_curl.terminate()
        p_curl.wait()

    rc_g = p_gunzip.wait()
    rc_d = p_deint.wait()
    rc_1 = p_gz1.wait()
    rc_2 = p_gz2.wait()
    for p in (f1, f2):
        if os.path.exists(p):
            os.remove(p)
    if not found or rc_g != 0 or rc_d != 0 or rc_1 != 0 or rc_2 != 0:
        raise RuntimeError(f"stream_deinterleave failed: found={found} gunzip={rc_g} deint={rc_d} gz1={rc_1} gz2={rc_2}")
    sz_gb = (os.path.getsize(r1_out) + os.path.getsize(r2_out)) / 1e9
    print(f"  streamed {sz_gb:.2f} GB (R1+R2) in {time.time() - t0:.1f}s", flush=True)


def prep_cami2_genomes(suite, samples, genomes_url, meta_dir, out_dir, existing_combined=None):
    """Build compact BLAST-safe genomes.fa (>G0001|c1) + genomes.tsv with exact per-sample coverages."""
    os.makedirs(out_dir, exist_ok=True)
    gfa = os.path.join(out_dir, "genomes.fa")
    gtsv = os.path.join(out_dir, "genomes.tsv")
    bdb = os.path.join(out_dir, "blastdb")
    if os.path.isfile(gfa) and os.path.isfile(gtsv) and os.path.isfile(bdb + ".nhr"):
        print(f"  reusing prepared genomes + BLAST DB in {out_dir}", flush=True)
        return gfa, gtsv, bdb

    # Map filename basename -> genome_ID
    f2g = {}
    g2f = {}
    with open(os.path.join(meta_dir, "genome_to_id.tsv")) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                gid = p[0].strip()
                bn = os.path.basename(p[1].strip())
                f2g[bn] = gid
                stem = os.path.splitext(bn)[0]
                f2g[stem] = gid
                g2f[gid] = bn

    # Load coverages across selected samples
    cov_by_g = {}
    for s in samples:
        cp = os.path.join(meta_dir, f"coverage{s}.tsv")
        with open(cp) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 2:
                    try:
                        v = float(p[1])
                    except ValueError:
                        continue
                    cov_by_g.setdefault(p[0].strip(), {})[s] = v

    # Keep only genomes with cov > 0 in at least one selected sample
    active_gids = {g for g, cmap in cov_by_g.items() if any(v > 0 for v in cmap.values())}
    print(f"  {suite}: {len(active_gids)} active genomes across samples {samples}", flush=True)

    rows = []
    short_map = {}
    n_written = 0
    total_bp = 0

    def emit_genome(orig_gid, fname, seq_list, out_fh):
        nonlocal n_written, total_bp
        n_written += 1
        sgid = f"G{n_written:04d}"
        short_map[orig_gid] = sgid
        g_bp = 0
        for ci, seq in enumerate(seq_list, 1):
            g_bp += len(seq)
            out_fh.write(f">{sgid}|c{ci}\n")
            for k in range(0, len(seq), 100):
                out_fh.write(seq[k:k + 100] + "\n")
        total_bp += g_bp
        rows.append((sgid, orig_gid, fname, len(seq_list), g_bp))

    with open(gfa, "w") as out_fh:
        if existing_combined and os.path.isfile(existing_combined):
            print(f"  converting existing {existing_combined} -> BLAST-safe {gfa}", flush=True)
            cur_stem = None
            cur_seqs = []
            buf = []
            for line in open(existing_combined):
                if line.startswith(">"):
                    if buf:
                        cur_seqs.append("".join(buf))
                        buf = []
                    stem = line[1:].split("|", 1)[0].strip()
                    if stem != cur_stem:
                        if cur_stem is not None:
                            og = f2g.get(cur_stem, cur_stem)
                            if og in active_gids:
                                emit_genome(og, g2f.get(og, cur_stem + ".fasta"), cur_seqs, out_fh)
                        cur_stem = stem
                        cur_seqs = []
                else:
                    buf.append(line.strip())
            if buf:
                cur_seqs.append("".join(buf))
            if cur_stem is not None:
                og = f2g.get(cur_stem, cur_stem)
                if og in active_gids:
                    emit_genome(og, g2f.get(og, cur_stem + ".fasta"), cur_seqs, out_fh)
        else:
            print(f"  streaming genomes from {genomes_url} ...", flush=True)
            p_curl = subprocess.Popen(["curl", "-sSL", "--retry", "5", genomes_url], stdout=subprocess.PIPE)
            with tarfile.open(fileobj=p_curl.stdout, mode="r|gz") as tf:
                for m in tf:
                    if not m.isfile():
                        continue
                    bn = os.path.basename(m.name)
                    if not (bn.endswith(".fasta") or bn.endswith(".fna") or bn.endswith(".fa")):
                        continue
                    stem = os.path.splitext(bn)[0]
                    og = f2g.get(bn) or f2g.get(stem)
                    if og is None or og not in active_gids:
                        continue
                    raw = tf.extractfile(m).read().decode("latin-1")
                    seqs, buf = [], []
                    for ln in raw.splitlines():
                        if ln.startswith(">"):
                            if buf:
                                seqs.append("".join(buf))
                                buf = []
                        elif ln.strip():
                            buf.append(ln.strip())
                    if buf:
                        seqs.append("".join(buf))
                    if seqs:
                        emit_genome(og, bn, seqs, out_fh)
            p_curl.wait()

    with open(gtsv, "w") as f:
        f.write("genome\torig_id\tfile\tn_seqs\tlength\t" + "\t".join(f"cov_s{s}" for s in samples) + "\n")
        for sgid, og, fname, ns, gbp in rows:
            covs = "\t".join(f"{cov_by_g.get(og, {}).get(s, 0.0):.4f}" for s in samples)
            f.write(f"{sgid}\t{og}\t{fname}\t{ns}\t{gbp}\t{covs}\n")

    print(f"  wrote {n_written} genomes ({total_bp/1e6:.1f} Mbp); building BLAST DB ...", flush=True)
    subprocess.run([MAKEBLASTDB, "-in", gfa, "-dbtype", "nucl", "-out", bdb, "-parse_seqids"],
                   stdout=subprocess.DEVNULL, check=True)
    return gfa, gtsv, bdb


def prep_cami3_genomes(samples, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    gfa = os.path.join(out_dir, "genomes.fa")
    gtsv = os.path.join(out_dir, "genomes.tsv")
    bdb = os.path.join(out_dir, "blastdb")
    if os.path.isfile(gfa) and os.path.isfile(gtsv) and os.path.isfile(bdb + ".nhr"):
        print(f"  reusing prepared CAMI III genomes + BLAST DB in {out_dir}", flush=True)
        return gfa, gtsv, bdb

    base_u = "https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy"
    cov_dir = os.path.join(out_dir, "coverages")
    if not os.path.isdir(cov_dir):
        subprocess.run(f"curl -sSL '{base_u}/coverages.tar.gz' | tar -xzf - -C '{out_dir}'", shell=True, check=True)
    s_arg = ",".join(str(s) for s in samples)
    print(f"  streaming CAMI III source_genomes.tar.gz for samples {s_arg} ...", flush=True)
    subprocess.run(
        f"curl -sSL --retry 5 '{base_u}/source_genomes.tar.gz' | python3 '{PREP_CAMI3}' - '{cov_dir}' '{out_dir}' '{s_arg}'",
        shell=True, check=True,
    )
    print("  building CAMI III BLAST DB ...", flush=True)
    subprocess.run([MAKEBLASTDB, "-in", gfa, "-dbtype", "nucl", "-out", bdb, "-parse_seqids"],
                   stdout=subprocess.DEVNULL, check=True)
    return gfa, gtsv, bdb


def run_one_sample(suite, s_num, reads_url, gfa, gtsv, bdb, deint_bin, kmadb, min_free_gb, extra_genes_fa=None):
    sid = f"{suite}_s{s_num}"
    res_dir = os.path.join(RESULTS_ROOT, suite, f"s{s_num}")
    score_json = os.path.join(res_dir, "score.v05_score.json")
    flanks_json = os.path.join(res_dir, "run.flanks.json")
    blast_tsv = os.path.join(res_dir, "genes_vs_genomes.tsv")
    genes_fa = os.path.join(res_dir, "genes.fa")

    if os.path.isfile(score_json) and os.path.isfile(flanks_json) and os.path.isfile(blast_tsv):
        print(f"[{sid}] already completed, skipping ({score_json})", flush=True)
        return {"sid": sid, "genes_fa": genes_fa, "flanks_json": flanks_json, "genomes_fa": gfa,
                "genomes_tsv": gtsv, "cov_col": f"cov_s{s_num}", "blast_tsv": blast_tsv}

    os.makedirs(res_dir, exist_ok=True)
    tmp_s = os.path.join(WORK_ROOT, f"tmp_{sid}")
    shutil.rmtree(tmp_s, ignore_errors=True)
    os.makedirs(tmp_s, exist_ok=True)

    check_disk(min_free_gb, f"before {sid}")
    r1 = os.path.join(tmp_s, "R1.fq.gz")
    r2 = os.path.join(tmp_s, "R2.fq.gz")
    t_all = time.time()
    print(f"\n================ [{sid}] ================", flush=True)
    stream_deinterleave(reads_url, r1, r2, deint_bin)

    # 1. Seed gene identification (KMA vs ResFinder, OR provided extra_genes_fa)
    if extra_genes_fa and os.path.isfile(extra_genes_fa):
        shutil.copy2(extra_genes_fa, genes_fa)
    else:
        kma_pre = os.path.join(tmp_s, "kma")
        t_k = time.time()
        subprocess.run(
            [KMA, "-ipe", r1, r2, "-o", kma_pre, "-t_db", kmadb, "-1t1", "-cge", "-ef", "-nf", "-t", "8"],
            stdout=subprocess.DEVNULL, stderr=open(os.path.join(res_dir, "kma.log"), "w"), check=True,
        )
        subprocess.run(
            [sys.executable, KMA_GENES, kma_pre, genes_fa],
            stdout=open(os.path.join(res_dir, "kma_genes.txt"), "w"), check=True,
        )
        shutil.copy2(kma_pre + ".res", os.path.join(res_dir, "kma.res"))
        n_genes = sum(1 for ln in open(genes_fa) if ln.startswith(">"))
        print(f"  KMA found {n_genes} ResFinder seed genes in {time.time() - t_k:.1f}s", flush=True)
        if n_genes == 0:
            shutil.rmtree(tmp_s, ignore_errors=True)
            print(f"  [{sid}] 0 ResFinder genes detected; skipping argjacent", flush=True)
            return None

    # 2. Run frozen ARGjacent v0.5.0
    run_pre = os.path.join(res_dir, "run")
    t_a = time.time()
    rc = subprocess.call(
        [ARGJACENT, "-g", genes_fa, "-1", r1, "-2", r2, "-o", run_pre, "-D", "10000", "-t", "8", "--max-mem", "7", "--keep"]
    )
    print(f"  ARGjacent v0.5.0 finished (rc={rc}) in {time.time() - t_a:.1f}s", flush=True)
    if rc != 0:
        shutil.rmtree(tmp_s, ignore_errors=True)
        raise RuntimeError(f"argjacent failed on {sid} with exit code {rc}")

    # 3. IMMEDIATELY DELETE raw FASTQ files before BLAST & scoring!
    for fdel in (r1, r2):
        if os.path.exists(fdel):
            os.remove(fdel)
    shutil.rmtree(tmp_s, ignore_errors=True)
    check_disk(min_free_gb - 2.0, f"after deleting FASTQs for {sid}")

    # 4. BLAST seed genes against suite genomes & score all 4 v0.5.0 outputs!
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    subprocess.run(
        [BLASTN, "-query", genes_fa, "-db", bdb, "-perc_identity", "85", "-max_target_seqs", "100000",
         "-evalue", "1e-20", "-num_threads", "8", "-outfmt", fmt],
        stdout=open(blast_tsv, "w"), check=True,
    )
    subprocess.run(
        [sys.executable, SCORE_SAMPLE, genes_fa, run_pre, gfa, gtsv, f"cov_s{s_num}", blast_tsv,
         os.path.join(res_dir, "score"), "--min-cov", "5.0", "--blastdb", bdb,
         "--minimap2", MINIMAP2, "--blastn", BLASTN, "--label", sid],
        check=True,
    )

    # 5. IMMEDIATELY DELETE run.tmp/ (leaving only compact JSON/TSV/GFA/FA outputs!)
    shutil.rmtree(run_pre + ".tmp", ignore_errors=True)
    # Also remove bulky flanks.fa? Wait — keep it if <5 MB, or check size; flanks.json already holds everything.
    if os.path.exists(run_pre + ".flanks.fa") and os.path.getsize(run_pre + ".flanks.fa") > 5_000_000:
        os.remove(run_pre + ".flanks.fa")
    print(f"[{sid}] DONE in {time.time() - t_all:.1f}s ({free_gb('/tmp'):.2f} GiB free)", flush=True)
    return {"sid": sid, "genes_fa": genes_fa, "flanks_json": flanks_json, "genomes_fa": gfa,
            "genomes_tsv": gtsv, "cov_col": f"cov_s{s_num}", "blast_tsv": blast_tsv}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", required=True, choices=["strain", "cami3", "marine", "plant"])
    ap.add_argument("--samples", required=True, help="comma-separated sample numbers, e.g. 0,10,20,30,40,50,60,70")
    ap.add_argument("--min-free-gb", type=float, default=6.5)
    ap.add_argument("--keep-work", action="store_true", help="keep suite genomes/BLAST DB in /tmp after suite finishes")
    ap.add_argument("--existing-strain-fa", default="/tmp/argj_val/strain_check/all_genomes.fa")
    a = ap.parse_args()

    samples = [int(x) for x in a.samples.split(",") if x.strip()]
    deint_bin, kmadb = ensure_tools()
    suite_work = os.path.join(WORK_ROOT, a.suite)
    suite_res = os.path.join(RESULTS_ROOT, a.suite)
    os.makedirs(suite_res, exist_ok=True)

    check_disk(a.min_free_gb, f"start suite {a.suite}")

    if a.suite == "strain":
        gfa, gtsv, bdb = prep_cami2_genomes(
            "strain", samples,
            "https://frl.publisso.de/data/frl:6425521/strain/strmgCAMI2_genomes.tar.gz",
            "/tmp/argj_val/meta/strain", suite_work,
            existing_combined=a.existing_strain_fa if os.path.isfile(a.existing_strain_fa) else None,
        )
        # Immediately remove original uncompressed strain_check duplicate once converted!
        if a.existing_strain_fa and os.path.isfile(a.existing_strain_fa):
            os.remove(a.existing_strain_fa)
        url_fn = lambda s: f"https://frl.publisso.de/data/frl:6425521/strain/short_read/strmgCAMI2_sample_{s}_reads.tar.gz"
    elif a.suite == "cami3":
        gfa, gtsv, bdb = prep_cami3_genomes(samples, suite_work)
        url_fn = lambda s: f"https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/sample_{s}_reads.tar.gz"
    elif a.suite == "marine":
        gfa, gtsv, bdb = prep_cami2_genomes(
            "marine", samples,
            "https://frl.publisso.de/data/frl:6425521/marine/marmgCAMI2_genomes.tar.gz",
            "/tmp/argj_val/meta/marine", suite_work,
        )
        url_fn = lambda s: f"https://frl.publisso.de/data/frl:6425521/marine/short_read/marmgCAMI2_sample_{s}_reads.tar.gz"
    elif a.suite == "plant":
        gfa, gtsv, bdb = prep_cami2_genomes(
            "plant", samples,
            "https://frl.publisso.de/data/frl:6425521/plant_associated/rhimgCAMI2_genomes.tar.gz",
            "/tmp/argj_val/meta/plant", suite_work,
        )
        url_fn = lambda s: f"https://frl.publisso.de/data/frl:6425521/plant_associated/short_read/rhimgCAMI2_sample_{s}_reads.tar.gz"

    manifest = []
    for s in samples:
        m = run_one_sample(a.suite, s, url_fn(s), gfa, gtsv, bdb, deint_bin, kmadb, a.min_free_gb)
        if m:
            manifest.append(m)

    # Run multi-sample `argjacent-cohort` (JSON-only mode!) at W=500 and W=1000 if >= 2 samples succeeded
    if len(manifest) >= 2:
        man_path = os.path.join(suite_res, "cohort_manifest.json")
        json.dump(manifest, open(man_path, "w"), indent=2)
        for W in (500, 1000):
            print(f"\n=== Running argjacent-cohort on {a.suite} (S={len(manifest)}, W={W}) ===", flush=True)
            subprocess.run(
                [sys.executable, SCORE_COHORT, "--manifest", man_path,
                 "--out-prefix", os.path.join(suite_res, f"cohort_W{W}"),
                 "--pair-bp", str(W), "--minimap2", MINIMAP2, "--label", f"{a.suite}_W{W}"],
                check=True,
            )

    if not a.keep_work:
        print(f"Cleaning temporary genome/BLAST work directory {suite_work} ...", flush=True)
        shutil.rmtree(suite_work, ignore_errors=True)
    check_disk(0.0, f"end suite {a.suite}")


if __name__ == "__main__":
    main()
