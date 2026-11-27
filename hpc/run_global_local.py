#!/usr/bin/env python3
"""Pipelined local runner for the ARGjacent Global Metagenome Benchmark (tuned for 16 GB Apple Silicon / laptops).

Architecture for a 16 GB RAM / 10-core MacBook Air M4:
  - 1 compute worker (--threads 8, --max-mem 8): keeps peak RAM safely under ~9-10 GB with zero swap.
  - 1 background prefetcher (--prefetch 2): downloads the next sample(s) while the current sample computes,
    gating every download behind a strict --min-free-gb disk safety floor (default 25 GB free).
  - Ephemeral scratch cleanup: deletes raw FASTQs right after fastp, trimmed FASTQs right after argjacent,
    and full assembly FASTAs right after extract_assembly_flanks.py (~1-15 MB retained per sample).
  - Fully resumable via <out-dir>/samples/<sample>/<sample>.done markers; checkpoints global comparison
    tables every --checkpoint-every samples and runs full cohort + comparison + annotation at completion.

Example:
    python3 hpc/run_global_local.py \\
        --manifest global_run/global_2023_2024_100.tsv \\
        --db-env global_run/dbs/db_paths.env \\
        --out global_run/out_2023_2024 \\
        --threads 8 --max-mem 8 --prefetch 2 --min-free-gb 25
"""
import argparse
import collections
import concurrent.futures
import fcntl
import json
import os
import queue
import shutil
import subprocess
import sys
import threading
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HPC_DIR = os.path.dirname(os.path.abspath(__file__))
ENA_SEM = threading.Semaphore(1)  # Strictly at most 1 concurrent TCP connection to ftp.sra.ebi.ac.uk to prevent 403/52 throttling


def build_env(extra_bins=None):
    env = dict(os.environ)
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    env["VECLIB_MAXIMUM_THREADS"] = "1"

    search_dirs = [
        os.path.join(REPO_ROOT, "bin"),
        "/opt/homebrew/bin",
        "/opt/miniconda3/bin",
        "/opt/miniconda3/envs/binner/bin",
    ]
    if extra_bins:
        search_dirs = list(extra_bins) + search_dirs
    cur_path = env.get("PATH", "").split(":")
    merged = []
    for d in search_dirs + cur_path:
        if d and os.path.isdir(d) and d not in merged:
            merged.append(d)
    env["PATH"] = ":".join(merged)
    return env


def load_db_env(path):
    vals = {}
    if not path or not os.path.isfile(path):
        return vals
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("export "):
                s = s[7:].strip()
            if "=" not in s:
                continue
            k, v = s.split("=", 1)
            v = v.strip().strip('"').strip("'")
            vals[k.strip()] = v
    return vals


def load_manifest(path):
    rows = []
    with open(path) as fh:
        hdr = [c.strip() for c in fh.readline().rstrip("\r\n").split("\t")]
        col = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            r = {c: (parts[i].strip() if i < len(parts) else "") for c, i in col.items()}
            if r.get("sample"):
                rows.append(r)
    return rows


def free_gb(path):
    return shutil.disk_usage(path).free / (1024 ** 3)


def download_file(url, dst, env, max_attempts=6):
    if os.path.isfile(url):
        shutil.copy2(url, dst)
        return
    tmp = dst + ".part"
    is_ena = "ftp.sra.ebi.ac.uk/" in url
    # Prefer plain HTTP first for ENA FTP mirrors (avoids ENA port-443 TLS drops & connection refusals)
    if url.startswith("https://ftp.sra.ebi.ac.uk/"):
        urls = ["http://" + url[len("https://"):], url]
    else:
        urls = [url]

    last_err = None
    for attempt in range(max_attempts):
        for u in urls:
            is_https_ena = u.startswith("https://ftp.sra.ebi.ac.uk/")
            cmd = [
                "curl", "-fL", "-C", "-",
                "--retry", "2" if is_https_ena else "4",
                "--retry-delay", "5", "--retry-all-errors",
                "--connect-timeout", "10" if is_https_ena else "30",
                "--speed-limit", "10000", "--speed-time", "60",
                "--max-time", "3600", "-sS", "-o", tmp, u,
            ]
            if is_ena:
                with ENA_SEM:
                    res = subprocess.run(cmd, env=env)
                    time.sleep(1.0)
            else:
                res = subprocess.run(cmd, env=env)
            if res.returncode == 0 and os.path.isfile(tmp) and os.path.getsize(tmp) > 0:
                os.replace(tmp, dst)
                return
            # Fast-fail 404 on S3 mirror so caller can fall back immediately
            if not is_ena and res.returncode == 22:
                shutil.rmtree(tmp, ignore_errors=True)
                if os.path.exists(tmp):
                    os.remove(tmp)
                raise subprocess.CalledProcessError(res.returncode, cmd)
            last_err = subprocess.CalledProcessError(res.returncode, cmd)
        time.sleep(min(90, 15 * (attempt + 1)))
    raise last_err


def extract_kma_genes(res_p, fsa_p, out_p, min_id, min_cov, ref_fa):
    kept = []
    if os.path.isfile(res_p):
        with open(res_p) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                p = [x.strip() for x in line.split("\t")]
                if len(p) < 6:
                    continue
                name, tid, tcov = p[0].split()[0], float(p[4]), float(p[5])
                if tid >= min_id and tcov >= min_cov:
                    kept.append(name)
    kept_set = set(kept)
    seqs = {}
    src = ref_fa if (ref_fa and os.path.isfile(ref_fa)) else fsa_p
    if kept_set and os.path.isfile(src):
        cur, buf = None, []
        with open(src) as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur in kept_set:
                        seqs[cur] = "".join(buf).upper().replace("-", "")
                    cur, buf = line[1:].strip().split()[0], []
                else:
                    buf.append(line.strip())
            if cur in kept_set:
                seqs[cur] = "".join(buf).upper().replace("-", "")
    n = 0
    with open(out_p, "w") as out:
        for name in kept:
            s = seqs.get(name, "")
            if len(s) >= 31:
                out.write(f">{name}\n{s}\n")
                n += 1
    return n


def atomic_copy(src, dst):
    part = dst + ".part"
    shutil.copy2(src, part)
    os.replace(part, dst)


def stage_sample(row, scratch_root, out_dir, args, env):
    """Stage SRA (via AWS S3 Open Data, falling back to single-stream ENA FASTQs) + pre-made assembly into scratch_root/<sample>/."""
    sid = row["sample"]
    run_acc = (row.get("run_accession") or sid).strip()
    s_tmp = os.path.join(scratch_root, sid)
    shutil.rmtree(s_tmp, ignore_errors=True)
    os.makedirs(s_tmp, exist_ok=True)

    # Wait if disk free space is below safety floor
    est_gb = float(row.get("fastq_gb") or 5.0)
    while free_gb(out_dir) - est_gb < args.min_free_gb and free_gb(out_dir) < args.min_free_gb + 12.0:
        time.sleep(3.0)

    t0 = time.time()
    r1_url = row.get("r1_url") or row.get("r1")
    r2_url = row.get("r2_url") or row.get("r2")
    r1_dst = os.path.join(s_tmp, f"{sid}_1.fastq.gz")
    r2_dst = os.path.join(s_tmp, f"{sid}_2.fastq.gz")
    sra_dst = ""

    fqdump = shutil.which("fasterq-dump", path=env["PATH"])
    used_s3 = False
    if fqdump and not args.no_s3 and run_acc.startswith(("SRR", "ERR", "DRR")):
        s3_url = f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{run_acc}/{run_acc}"
        cand_sra = os.path.join(s_tmp, f"{sid}.sra")
        try:
            download_file(s3_url, cand_sra, env, max_attempts=3)
            sra_dst = cand_sra
            r1_dst = r2_dst = ""
            used_s3 = True
        except Exception:
            if os.path.exists(cand_sra):
                os.remove(cand_sra)

    if not used_s3:
        download_file(r1_url, r1_dst, env)
        download_file(r2_url, r2_dst, env)

    asm_url = row.get("logan_url") if (args.use_logan and row.get("logan_url")) else (row.get("assembly_url") or row.get("assembly"))
    asm_label = f"logan:{sid}" if (args.use_logan and row.get("logan_url")) else (row.get("erz_accession") or "assembly")
    asm_dst = ""
    cached_asm_dir = ""
    if args.asm_cache_dir and not args.use_logan:
        cdir = os.path.join(args.asm_cache_dir, sid)
        if (os.path.isfile(os.path.join(cdir, f"{sid}.asm_flanks.json"))
                and os.path.isfile(os.path.join(cdir, f"{sid}.asm_flanks.fa"))
                and os.path.isfile(os.path.join(cdir, f"{sid}.asm_contexts.tsv"))):
            cached_asm_dir = cdir
    if not args.skip_assembly and asm_url and not cached_asm_dir:
        ext = ".fa.zst" if asm_url.endswith(".zst") else (".fa" if asm_url.endswith((".fa", ".fna", ".fasta")) else ".fasta.gz")
        asm_dst = os.path.join(s_tmp, f"{sid}{ext}")
        download_file(asm_url, asm_dst, env)

    dl_s = time.time() - t0
    return {
        "row": row,
        "tmp": s_tmp,
        "sra": sra_dst,
        "r1": r1_dst,
        "r2": r2_dst,
        "r1_url": r1_url,
        "r2_url": r2_url,
        "asm": asm_dst,
        "cached_asm_dir": cached_asm_dir,
        "asm_label": asm_label,
        "dl_s": dl_s,
    }


def process_staged(staged, out_dir, kma_db, kma_ref_fa, shared_genes, args, env):
    """Run fastp -> KMA -> ARGjacent -> extract_assembly_flanks -> compare, purging scratch immediately."""
    row = staged["row"]
    sid = row["sample"]
    s_tmp = staged["tmp"]
    s_out = os.path.join(out_dir, "samples", sid)
    os.makedirs(s_out, exist_ok=True)
    local_pre = os.path.join(s_tmp, sid)

    try:
        t_comp0 = time.time()
        r1, r2 = staged["r1"], staged["r2"]
        if staged.get("sra") and os.path.isfile(staged["sra"]):
            fqdump = shutil.which("fasterq-dump", path=env["PATH"])
            unpacked_ok = False
            if fqdump:
                while free_gb(out_dir) < max(16.0, args.min_free_gb - 6.0):
                    time.sleep(3.0)
                rc_fq = subprocess.run(
                    [fqdump, staged["sra"], "-O", s_tmp, "-t", s_tmp, "-e", str(args.threads), "--split-3", "--skip-technical"],
                    env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                )
                if os.path.exists(staged["sra"]):
                    os.remove(staged["sra"])
                unmated = os.path.join(s_tmp, f"{sid}.fastq")
                if os.path.exists(unmated):
                    os.remove(unmated)
                cand_r1 = os.path.join(s_tmp, f"{sid}_1.fastq")
                cand_r2 = os.path.join(s_tmp, f"{sid}_2.fastq")
                if rc_fq.returncode == 0 and os.path.isfile(cand_r1) and os.path.isfile(cand_r2) and os.path.getsize(cand_r1) > 0:
                    r1, r2 = cand_r1, cand_r2
                    unpacked_ok = True
                else:
                    for p in (cand_r1, cand_r2):
                        if os.path.exists(p):
                            os.remove(p)
            if not unpacked_ok:
                r1 = os.path.join(s_tmp, f"{sid}_1.fastq.gz")
                r2 = os.path.join(s_tmp, f"{sid}_2.fastq.gz")
                download_file(staged["r1_url"], r1, env)
                download_file(staged["r2_url"], r2, env)

        # 1. Optional fastp trimming (and immediate deletion of raw FASTQs)
        if not args.skip_fastp and shutil.which("fastp", path=env["PATH"]):
            tr1 = os.path.join(s_tmp, f"{sid}_trim_1.fq.gz")
            tr2 = os.path.join(s_tmp, f"{sid}_trim_2.fq.gz")
            subprocess.run(
                [
                    "fastp", "-i", r1, "-I", r2, "-o", tr1, "-O", tr2,
                    "--detect_adapter_for_pe", "-q", "20", "-l", "80",
                    "-w", str(args.threads),
                    "-j", f"{local_pre}.fastp.json", "-h", "/dev/null",
                ],
                env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
            )
            os.remove(r1)
            os.remove(r2)
            r1, r2 = tr1, tr2

        # 2. Target genes via KMA pre-screen (or shared FASTA)
        n_genes = -1
        if kma_db:
            kma_pre = os.path.join(s_tmp, "kma_out")
            subprocess.run(
                ["kma", "-ipe", r1, r2, "-o", kma_pre, "-t_db", kma_db, "-t", str(args.threads), "-1t1", "-cge", "-ef", "-nf"],
                env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
            )
            genes_fa = os.path.join(s_tmp, f"{sid}.genes.fa")
            n_genes = extract_kma_genes(
                f"{kma_pre}.res", f"{kma_pre}.fsa", genes_fa,
                args.kma_min_id, args.kma_min_cov, kma_ref_fa,
            )
            if os.path.isfile(f"{kma_pre}.res"):
                atomic_copy(f"{kma_pre}.res", os.path.join(s_out, f"{sid}.kma.res"))
            if n_genes == 0:
                # Clean zero-gene completion
                open(os.path.join(s_out, f"{sid}.summary.tsv"), "w").write(
                    "gene\tside\tstatus\tcontexts\tcertain_bp\tlongest_bp\tcontexts_certain_1kb\tread_depth\tgene_end\tcertain_end_bp\n"
                )
                open(os.path.join(s_out, f"{sid}.contexts.tsv"), "w").write(
                    "gene\tside\tcontext\tlength_bp\tcertain_bp\tread_depth\tstop\tcertain_stop\tother_genes\tgene_end\tmain\tcertain_end_bp\tshared_end_bp\n"
                )
                open(os.path.join(s_out, f"{sid}.flanks.fa"), "w").write("")
                open(os.path.join(s_out, f"{sid}.pairs.tsv"), "w").write(
                    "gene\tpair_id\tleft_context\tright_context\tleft_certain_bp\tright_certain_bp\tleft_depth\tright_depth\ttier\tscore\tnote\n"
                )
                json.dump({"version": "0.6.0", "flank_bp": args.flank_bp, "aw_stats": {}, "parameters": {}, "genes": {}},
                          open(os.path.join(s_out, f"{sid}.flanks.json"), "w"))
                json.dump({"version": "asm_flanks_v1", "source": staged["asm_label"], "flank_bp": args.flank_bp, "genes": {}},
                          open(os.path.join(s_out, f"{sid}.asm_flanks.json"), "w"))
                open(os.path.join(s_out, f"{sid}.log"), "w").write("argjacent local: 0 genes detected by KMA\n")
                open(os.path.join(s_out, f"{sid}.done"), "w").write(time.strftime("%Y-%m-%dT%H:%M:%SZ\n", time.gmtime()))
                return {
                    "sample": sid, "status": "ZERO_GENES", "dl_s": staged["dl_s"],
                    "argj_s": 0.0, "asm_s": 0.0, "tot_s": time.time() - t_comp0,
                    "n_genes": 0, "argj_1kb": 0, "asm_1kb": 0, "multi_wins": 0, "past_cliff": 0,
                }
            atomic_copy(genes_fa, os.path.join(s_out, f"{sid}.genes.fa"))
            active_genes = os.path.join(s_out, f"{sid}.genes.fa")
        else:
            active_genes = shared_genes

        # 3. Run ARGjacent and immediately delete FASTQs
        argj_bin = shutil.which("argjacent", path=env["PATH"]) or os.path.join(REPO_ROOT, "bin", "argjacent")
        t_a0 = time.time()
        cmd = [
            argj_bin, "-g", active_genes, "-1", r1, "-2", r2, "-o", local_pre,
            "-D", str(args.flank_bp), "-t", str(args.threads), "--max-mem", str(args.max_mem),
        ]
        if args.extra_args:
            cmd.extend(args.extra_args.split())
        subprocess.run(cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        argj_s = time.time() - t_a0
        for p in (r1, r2):
            if os.path.exists(p):
                os.remove(p)

        # 4. Extract pre-made assembly flanks + per-sample comparison (then delete assembly FASTA!)
        asm_s = 0.0
        argj_1kb = asm_1kb = multi_wins = past_cliff = 0
        has_asm_flanks = False
        if staged.get("cached_asm_dir"):
            t_m0 = time.time()
            for ext in (".asm_flanks.json", ".asm_contexts.tsv", ".asm_flanks.fa"):
                atomic_copy(os.path.join(staged["cached_asm_dir"], sid + ext), local_pre + ext)
            has_asm_flanks = True
        elif staged["asm"] and os.path.isfile(staged["asm"]):
            t_m0 = time.time()
            subprocess.run(
                [
                    sys.executable, os.path.join(HPC_DIR, "extract_assembly_flanks.py"),
                    "-g", active_genes, "-c", staged["asm"], "-o", local_pre,
                    "-D", str(args.flank_bp), "--aligner", args.asm_aligner,
                    "-t", str(args.threads), "--source-label", staged["asm_label"],
                ],
                env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
            )
            os.remove(staged["asm"])
            has_asm_flanks = True
        if has_asm_flanks:
            subprocess.run(
                [
                    sys.executable, os.path.join(HPC_DIR, "compare_argjacent_vs_assembly.py"),
                    "--argj-json", f"{local_pre}.flanks.json",
                    "--asm-json", f"{local_pre}.asm_flanks.json",
                    "--sample", sid, "-o", f"{local_pre}.cmp",
                ],
                env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
            )
            asm_s = time.time() - t_m0
            if os.path.isfile(f"{local_pre}.cmp_summary.tsv"):
                with open(f"{local_pre}.cmp_summary.tsv") as fh:
                    h = fh.readline().rstrip("\n").split("\t")
                    vals = fh.readline().rstrip("\n").split("\t")
                    if len(vals) == len(h):
                        d = dict(zip(h, vals))
                        argj_1kb = int(d.get("argj_1kb") or 0)
                        asm_1kb = int(d.get("asm_1kb") or 0)
                        multi_wins = int(d.get("argj_multi_vs_asm_single") or 0)
                        past_cliff = int(d.get("asm_past_repeat_cliff") or 0)

        # 5. Stage compact outputs atomically
        for ext in (
            ".summary.tsv", ".contexts.tsv", ".flanks.fa", ".flanks.json", ".pairs.tsv", ".gfa", ".log",
            ".fastp.json", ".asm_flanks.json", ".asm_contexts.tsv", ".asm_flanks.fa",
            ".cmp_per_side.tsv", ".cmp_summary.tsv",
        ):
            src = local_pre + ext
            if os.path.isfile(src):
                atomic_copy(src, os.path.join(s_out, sid + ext))

        open(os.path.join(s_out, f"{sid}.genes_path.txt"), "w").write(active_genes + "\n")
        open(os.path.join(s_out, f"{sid}.done"), "w").write(time.strftime("%Y-%m-%dT%H:%M:%SZ\n", time.gmtime()))
        return {
            "sample": sid, "status": "OK", "dl_s": staged["dl_s"],
            "argj_s": argj_s, "asm_s": asm_s, "tot_s": time.time() - t_comp0,
            "n_genes": n_genes, "argj_1kb": argj_1kb, "asm_1kb": asm_1kb,
            "multi_wins": multi_wins, "past_cliff": past_cliff,
        }
    finally:
        shutil.rmtree(s_tmp, ignore_errors=True)


def run_post_analysis(manifest_path, out_dir, db_vals, args, env, run_cohort=True, run_annotate=True):
    """Run summarize_run, per-BioProject argjacent-cohort, global assembly comparison, and annotation."""
    os.makedirs(os.path.join(out_dir, "summary"), exist_ok=True)
    os.makedirs(os.path.join(out_dir, "comparison"), exist_ok=True)

    subprocess.run(
        [sys.executable, os.path.join(HPC_DIR, "summarize_run.py"),
         "--samples", manifest_path, "--out-dir", out_dir, "--summary-dir", os.path.join(out_dir, "summary")],
        env=env, check=False,
    )

    if not args.skip_assembly:
        subprocess.run(
            [sys.executable, os.path.join(HPC_DIR, "compare_argjacent_vs_assembly.py"),
             "--manifest", manifest_path, "--out-dir", out_dir,
             "-o", os.path.join(out_dir, "comparison", "global_cmp"), "--export-fasta"],
            env=env, check=False,
        )

    if run_cohort:
        cohort_bin = shutil.which("argjacent-cohort", path=env["PATH"]) or os.path.join(REPO_ROOT, "bin", "argjacent-cohort")
        by_study = collections.defaultdict(list)
        for r in load_manifest(manifest_path):
            sid = r["sample"]
            st = r.get("study_accession") or "ALL"
            sdir = os.path.join(out_dir, "samples", sid)
            jp = os.path.join(sdir, f"{sid}.flanks.json")
            if not (os.path.isfile(jp) and os.path.getsize(jp) > 0):
                continue
            try:
                if not json.load(open(jp)).get("genes"):
                    continue
            except Exception:
                continue
            gp = ""
            gf = os.path.join(sdir, f"{sid}.genes_path.txt")
            if os.path.isfile(gf):
                cand = open(gf).read().strip()
                if cand and os.path.isfile(cand):
                    gp = cand
            by_study[st].append((sid, os.path.abspath(jp), os.path.abspath(gp) if gp else "", ""))
        for st, rows in sorted(by_study.items()):
            if len(rows) < args.cohort_min_samples:
                continue
            cdir = os.path.join(out_dir, "cohorts", st)
            os.makedirs(cdir, exist_ok=True)
            sheet = os.path.join(cdir, "cohort_samplesheet.tsv")
            with open(sheet, "w") as out:
                out.write("sample\tflanks_json\tgenes_fa\treads_tsv\n")
                for rr in rows:
                    out.write("\t".join(rr) + "\n")
            subprocess.run(
                [cohort_bin, "--sample-sheet", sheet, "-o", os.path.join(cdir, st), "--pair-bp", "500", "-t", str(min(4, args.threads))],
                env=env, check=False,
            )

    if run_annotate and db_vals:
        argj_fa = os.path.join(out_dir, "comparison", "global_cmp_argjacent_certified.fa")
        asm_fa = os.path.join(out_dir, "comparison", "global_cmp_assembly_flanks.fa")
        if os.path.isfile(argj_fa) and os.path.getsize(argj_fa) > 0:
            ann_dir = os.path.join(out_dir, "annotation")
            cmd = [
                sys.executable, os.path.join(HPC_DIR, "annotate_flanks.py"),
                "--argj-fa", argj_fa, "--manifest", manifest_path,
                "-o", ann_dir, "-t", str(args.threads), "--min-len", "1000",
            ]
            if os.path.isfile(asm_fa) and os.path.getsize(asm_fa) > 0:
                cmd.extend(["--asm-fa", asm_fa])
            for env_k, flag in (
                ("UNIVEC_DB", "--univec-db"),
                ("ARG_DB", "--arg-db"),
                ("ISFINDER_DB", "--isfinder-db"),
                ("PLASMIDFINDER_DB", "--plasmidfinder-db"),
                ("MOBILEOG_DB", "--mobileog-db"),
                ("GENOMAD_DB", "--genomad-db"),
                ("KRAKEN2_DB", "--kraken2-db"),
                ("MMSEQS_TAX_DB", "--mmseqs-tax-db"),
            ):
                val = db_vals.get(env_k, "")
                if val and os.path.exists(val):
                    cmd.extend([flag, val])
            if shutil.which("mmseqs", path=env["PATH"]) or shutil.which("vsearch", path=env["PATH"]):
                cmd.append("--cluster-hgt")
            else:
                # Allow lightweight HGT clustering on local laptop if < 2,500 flanks
                cmd.extend(["--cluster-hgt", "--allow-python-fallback"])
            subprocess.run(cmd, env=env, check=False)


def main():
    ap = argparse.ArgumentParser(description="Pipelined MacBook / local runner for ARGjacent global metagenome benchmarking.")
    ap.add_argument("--manifest", required=True, help="Manifest TSV from hpc/fetch_global_manifest.py")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--db-env", help="db_paths.env from hpc/setup_annotation_dbs.sh (provides KMA_DB + annotation DBs)")
    ap.add_argument("--kma-db", help="Override KMA database prefix")
    ap.add_argument("--kma-ref-fa", help="Override KMA reference FASTA")
    ap.add_argument("--genes", help="Shared genes FASTA (used only if no KMA DB)")
    ap.add_argument("--limit", type=int, default=0, help="Process at most first N manifest rows (0 = all)")
    ap.add_argument("--threads", type=int, default=8, help="Threads for fastp/kma/argjacent/minimap2 (default: 8)")
    ap.add_argument("--max-mem", type=float, default=8.0, help="ARGjacent --max-mem GB cap (default: 8.0 for 16 GB Mac)")
    ap.add_argument("--prefetch", type=int, default=2, help="Max staged samples queued ahead of compute (default: 2)")
    ap.add_argument("--download-workers", type=int, default=2, help="Concurrent sample download threads (default: 2)")
    ap.add_argument("--min-free-gb", type=float, default=25.0, help="Pause prefetch if SSD free space < N GB (default: 25)")
    ap.add_argument("--flank-bp", type=int, default=10000, help="Target flank bp (default: 10000)")
    ap.add_argument("--kma-min-id", type=float, default=90.0)
    ap.add_argument("--kma-min-cov", type=float, default=60.0)
    ap.add_argument("--asm-aligner", default="auto", choices=("auto", "minimap2", "blastn", "kmer"))
    ap.add_argument("--use-logan", action="store_true", help="Compare against Logan S3 contigs instead of ENA ERZ")
    ap.add_argument("--no-s3", action="store_true", help="Disable AWS S3 SRA Open Data staging (use ENA FASTQs directly)")
    ap.add_argument("--asm-cache-dir", default="", help="Optional existing samples/ dir with pre-extracted .asm_flanks.* to reuse")
    ap.add_argument("--skip-fastp", action="store_true")
    ap.add_argument("--skip-assembly", action="store_true")
    ap.add_argument("--cohort-min-samples", type=int, default=4)
    ap.add_argument("--checkpoint-every", type=int, default=10, help="Refresh global comparison tables every N samples (default: 10)")
    ap.add_argument("--extra-args", default="")
    args = ap.parse_args()

    env = build_env()
    db_vals = load_db_env(args.db_env)
    kma_db = args.kma_db or db_vals.get("KMA_DB", "")
    kma_ref_fa = args.kma_ref_fa or db_vals.get("KMA_REF_FA", "")
    shared_genes = args.genes or ""
    if not (kma_db or shared_genes):
        sys.exit("ERROR: provide --db-env, --kma-db, or --genes.")

    os.makedirs(args.out, exist_ok=True)
    lock_fh = open(os.path.join(args.out, ".runner.lock"), "w")
    try:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit(f"ERROR: another run_global_local.py instance is already active on {args.out}")
    lock_fh.write(f"{os.getpid()}\n")
    lock_fh.flush()

    scratch_root = os.path.join(args.out, ".scratch")
    # Purge any stale scratch left from a previous interrupted run
    shutil.rmtree(scratch_root, ignore_errors=True)
    os.makedirs(scratch_root, exist_ok=True)

    rows = load_manifest(args.manifest)
    if args.limit > 0:
        rows = rows[:args.limit]

    todo = []
    already = 0
    for r in rows:
        sid = r["sample"]
        sdir = os.path.join(args.out, "samples", sid)
        if os.path.isfile(os.path.join(sdir, f"{sid}.done")) and os.path.isfile(os.path.join(sdir, f"{sid}.flanks.json")):
            already += 1
        else:
            todo.append(r)

    print(f"=== ARGjacent Local Pipelined Global Runner ===", flush=True)
    print(f"  Manifest      : {args.manifest} ({len(rows)} rows; {already} already done, {len(todo)} remaining)", flush=True)
    print(f"  Output dir    : {args.out} (free disk: {free_gb(args.out):.1f} GiB, safety floor: {args.min_free_gb:.0f} GiB)", flush=True)
    print(f"  Compute config: 1 worker x {args.threads} threads, --max-mem {args.max_mem} GB, prefetch={args.prefetch}, dl_workers={args.download_workers}", flush=True)
    print(f"  Tools found   : fastp={bool(shutil.which('fastp', path=env['PATH']))}, kma={bool(shutil.which('kma', path=env['PATH']))}, minimap2={bool(shutil.which('minimap2', path=env['PATH']))}", flush=True)

    if not todo:
        print("All samples already complete! Running final aggregation + annotation...", flush=True)
        run_post_analysis(args.manifest, args.out, db_vals, args, env, run_cohort=True, run_annotate=True)
        return

    q = queue.Queue(maxsize=max(1, args.prefetch))
    stop_ev = threading.Event()
    SENTINEL = object()
    todo_iter = iter(todo)
    iter_lock = threading.Lock()
    dl_workers = max(1, min(args.download_workers, len(todo)))

    def prefetcher():
        while not stop_ev.is_set():
            with iter_lock:
                r = next(todo_iter, None)
            if r is None:
                return
            try:
                st = stage_sample(r, scratch_root, args.out, args, env)
                q.put(("OK", st))
            except Exception as e:
                shutil.rmtree(os.path.join(scratch_root, r["sample"]), ignore_errors=True)
                q.put(("DL_ERR", (r, str(e))))

    def dl_coordinator():
        threads = [threading.Thread(target=prefetcher, daemon=True) for _ in range(dl_workers)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        q.put(SENTINEL)

    t_dl = threading.Thread(target=dl_coordinator, daemon=True)
    t_dl.start()

    fail_log = os.path.join(args.out, "failed_samples.tsv")
    if not os.path.isfile(fail_log):
        open(fail_log, "w").write("sample\tstage\terror\n")

    done_n = already
    ok_this_run = 0
    t_start = time.time()

    try:
        while True:
            item = q.get()
            if item is SENTINEL:
                break
            kind, payload = item
            if kind == "DL_ERR":
                r, err = payload
                done_n += 1
                print(f"[{done_n:>3d}/{len(rows)}] {r['sample']} ({r.get('biome_group','')}) DOWNLOAD FAILED: {err}", flush=True)
                with open(fail_log, "a") as fh:
                    fh.write(f"{r['sample']}\tdownload\t{err}\n")
                continue

            sid = payload["row"]["sample"]
            biome = payload["row"].get("biome_group", "")
            country = payload["row"].get("country", "")
            try:
                res = process_staged(payload, args.out, kma_db, kma_ref_fa, shared_genes, args, env)
                done_n += 1
                ok_this_run += 1
                elapsed_min = (time.time() - t_start) / 60.0
                rate = ok_this_run / max(elapsed_min / 60.0, 1e-6)
                print(
                    f"[{done_n:>3d}/{len(rows)}] {sid} ({biome}/{country}) {res['status']} | "
                    f"dl={res['dl_s']:.0f}s argj={res['argj_s']:.0f}s asm={res['asm_s']:.0f}s | "
                    f"genes={res['n_genes']} >=1kb(ARGj={res['argj_1kb']},Asm={res['asm_1kb']}) "
                    f"multi={res['multi_wins']} cliff={res['past_cliff']} | "
                    f"free={free_gb(args.out):.1f}GiB ({rate:.1f} smp/hr)",
                    flush=True,
                )
                if args.checkpoint_every > 0 and ok_this_run % args.checkpoint_every == 0:
                    run_post_analysis(args.manifest, args.out, db_vals, args, env, run_cohort=False, run_annotate=True)
            except Exception as e:
                done_n += 1
                print(f"[{done_n:>3d}/{len(rows)}] {sid} ({biome}) PROCESS FAILED: {e}", flush=True)
                with open(fail_log, "a") as fh:
                    fh.write(f"{sid}\tprocess\t{e}\n")
    except KeyboardInterrupt:
        stop_ev.set()
        print("\nInterrupted by user; cleaning scratch and running checkpoint summary...", flush=True)
    finally:
        stop_ev.set()
        shutil.rmtree(scratch_root, ignore_errors=True)

    print("\n=== Running final summary, per-BioProject cohort, assembly comparison & annotation ===", flush=True)
    run_post_analysis(args.manifest, args.out, db_vals, args, env, run_cohort=True, run_annotate=True)


if __name__ == "__main__":
    main()
