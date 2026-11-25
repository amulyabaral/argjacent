#!/usr/bin/env python3
"""Storage-safe streaming runner and PacBio HiFi validator for frozen ARGjacent v0.5.0
on fresh unseen real paired Illumina + PacBio HiFi metagenomes (Suite 4).

Storage guarantees:
  - Streams Illumina R1/R2 directly from ENA HTTPS, trims with fastp (--detect_adapter_for_pe),
    deletes raw FASTQs immediately after fastp, deletes trimmed FASTQs immediately after ARGjacent.
  - Streams PacBio HiFi FASTQ.gz directly from ENA HTTPS through gzip -dc | minimap2 | samtools fasta
    so ZERO full HiFi runs ever touch disk (keeping only compact hifi_near.fa.gz, ~5-30 MB).
  - Evaluates:
      1. One-sided certified flanks (`certain_bp` and `certain_end_bp`) with the exact frozen
         `benchmarks/heldout/scorer/claims_lr3.py` v3.1 + `real_triage.py`.
      2. Lossless GFA v1.0 compaction & 100% byte-for-byte round-trip check.
      3. Single-sample two-sided pairing (W=500 and W=1,000, using v0.5.0 `depth_bins`) vs
         Cartesian, Main-to-Main, and Unguarded Greedy baselines against real single-molecule HiFi reads.
      4. Multi-sample `argjacent-cohort` (Part A pairing + Part B post-repeat exits in JSON-only mode)
         on longitudinal cohorts, scored directly against real single-molecule HiFi reads.
"""
import argparse
import collections
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "benchmarks", "v05_validation"))
sys.path.insert(0, os.path.join(REPO, "benchmarks", "dev-tools", "experiments", "gpc"))

from argjacent.pair import extract_trunks, index_sample_reads, pair_gene_single_sample
from score_v05_cami import cp_ci, naive_pairs, verify_gfa_file
from truth_2side import edit_sim_banded, extract_two_sided_truth, seq_ident_kmer
sys.path.insert(0, os.path.join(REPO, "benchmarks", "heldout", "scorer"))
import claims_lr3 as c3

RESFINDER_FA = os.environ.get("RESFINDER_FA", "resfinder_db/all.fsa")
KMA = os.environ.get("KMA_BIN") or shutil.which("kma") or "kma"
KMA_INDEX = os.environ.get("KMA_INDEX_BIN") or shutil.which("kma_index") or "kma_index"
FASTP = os.environ.get("FASTP_BIN") or shutil.which("fastp") or "fastp"
MINIMAP2 = os.environ.get("MINIMAP2_BIN") or shutil.which("minimap2") or "minimap2"
SAMTOOLS = os.environ.get("SAMTOOLS_BIN") or shutil.which("samtools") or "samtools"
FASTQ_DUMP = os.environ.get("FASTQ_DUMP_BIN") or shutil.which("fastq-dump") or "fastq-dump"
GET_RUN = os.path.join(REPO, "benchmarks", "real_metagenome", "get_run.sh")
ARGJACENT = os.path.join(REPO, "bin", "argjacent")
ARGJACENT_COHORT = os.path.join(REPO, "bin", "argjacent-cohort")
KMA_GENES = os.path.join(REPO, "benchmarks", "cami", "kma_genes.py")
TILE_LR = os.path.join(REPO, "benchmarks", "real_metagenome", "tile_longreads.py")
CLAIMS_LR3 = os.path.join(REPO, "benchmarks", "heldout", "scorer", "claims_lr3.py")
REAL_TRIAGE = os.path.join(REPO, "benchmarks", "heldout", "scorer", "real_triage.py")

WORK_ROOT = "/tmp/argj_val/real_work"
RESULTS_ROOT = os.path.join(REPO, "benchmarks", "v05_validation", "results", "real_hifi")

# Pre-registered unseen real paired Illumina + PacBio HiFi samples (all never seen by any ARGjacent version!)
SAMPLES = [
    {"id": "INF134_8G1", "group": "infant_cohort", "project": "PRJNA1365828", "sr": "SRR36082382", "lr": "SRR36082377", "lr_spots": 2330356, "desc": "Infant gut INF1340011 t8 (NovaSeq + Revio HiFi 11.5 Gbp)"},
    {"id": "INF133_3G1", "group": "infant_cohort", "project": "PRJNA1365828", "sr": "SRR36082375", "lr": "SRR36082374", "lr_spots": 1613835, "desc": "Infant gut INF1330004 t3 (NovaSeq + Revio HiFi 7.8 Gbp)"},
    {"id": "INF133_8G1", "group": "infant_cohort", "project": "PRJNA1365828", "sr": "SRR36082386", "lr": "SRR36082371", "lr_spots": 1083327, "desc": "Infant gut INF1330004 t8 (NovaSeq + Revio HiFi 7.7 Gbp)"},
    {"id": "INF134_2G1", "group": "infant_cohort", "project": "PRJNA1365828", "sr": "SRR36082367", "lr": "SRR36082381", "lr_spots": 912611,  "desc": "Infant gut INF1340011 t2 (NovaSeq + Revio HiFi 7.2 Gbp)"},
    {"id": "INF134_3G1", "group": "infant_cohort", "project": "PRJNA1365828", "sr": "SRR36082366", "lr": "SRR36082384", "lr_spots": 1236486, "desc": "Infant gut INF1340011 t3 (NovaSeq + Revio HiFi 7.1 Gbp)"},
    {"id": "INF134_4G1", "group": "infant_cohort", "project": "PRJNA1365828", "sr": "SRR36082365", "lr": "SRR36082380", "lr_spots": 1233970, "desc": "Infant gut INF1340011 t4 (NovaSeq + Revio HiFi 6.8 Gbp)"},
    {"id": "US_Adult2T", "group": "adult_gut",     "project": "PRJNA754443",  "sr": "SRR15489030", "lr": "SRR15489017", "lr_spots": 1095633, "desc": "US adult gut 2_treatment (NovaSeq + Sequel II HiFi 7.4 Gbp)"},
    {"id": "UK_Adult",   "group": "adult_gut",     "project": "PRJEB90666",   "sr": "ERR15463326", "lr": "ERR15479330", "lr_spots": 987561,  "desc": "UK adult gut 21BR240s1 (NextSeq 2000 + Revio HiFi 6.7 Gbp)"},
    {"id": "Zymo_D6331", "group": "mock_gut",      "project": "PRJNA1446191", "sr": "SRR37882064", "lr": "SRR13128014", "lr_spots": 1978852, "desc": "ZymoBIOMICS Gut Standard D6331 (5 E. coli strains; NextSeq 2000 + Sequel II HiFi 18 Gbp)"},
]


def ensure_sra_nocache():
    """Dedicated HOME for sra-tools with local caching disabled so streamed HiFi never fills disk."""
    shome = "/tmp/argj_val/sra_nocache"
    ncbi = os.path.join(shome, ".ncbi")
    os.makedirs(ncbi, exist_ok=True)
    cfg = os.path.join(ncbi, "user-settings.mkfg")
    if not os.path.isfile(cfg):
        with open(cfg, "w") as f:
            f.write(
                '/LIBS/GUID = "a347909f-3e09-4323-9c73-015941985b61"\n'
                '/config/default = "false"\n'
                '/repository/remote/disabled = "false"\n'
                '/repository/user/cache-disabled = "true"\n'
            )
    return shome


def free_gb(path="/tmp"):
    return shutil.disk_usage(path).free / (1024 ** 3)


def check_disk(min_gb, stage=""):
    fg = free_gb("/tmp")
    print(f"[disk] {stage}: {fg:.2f} GiB free", flush=True)
    if fg < min_gb:
        raise RuntimeError(f"Low disk ({fg:.2f} GiB free < {min_gb:.1f} GiB) at {stage}")


def ena_fastq_urls(acc):
    u = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes"
    lines = urllib.request.urlopen(u, timeout=30).read().decode().strip().splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"ENA returned no records for {acc}")
    r = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    urls = ["https://" + x for x in r["fastq_ftp"].split(";") if x]
    md5s = [x for x in r["fastq_md5"].split(";") if x]
    return urls, md5s


def download_file(url, dest):
    subprocess.run(["curl", "-sSL", "--fail", "--retry", "5", "--retry-delay", "5", "-o", dest, url], check=True)


def match_side_hifi(pred_seq, hifi_seq, side, W):
    p = pred_seq[-W:] if side == "left" else pred_seq[:W]
    h = hifi_seq[-W:] if side == "left" else hifi_seq[:W]
    if len(p) < W or len(h) < W:
        return False
    if seq_ident_kmer(p, h, k=15) >= 0.78:
        return True
    return edit_sim_banded(p, h) >= 0.94


def score_pairing_hifi(doc_genes, obs_map, preds_by_gene, W):
    fam_reads = collections.defaultdict(dict)
    for g2, rl2 in obs_map.items():
        f2 = g2.split("_")[0]
        for r in rl2:
            rid = r["read"]
            if rid not in fam_reads[f2] or r["id_gene"] > fam_reads[f2][rid]["id_gene"]:
                fam_reads[f2][rid] = r

    tot = chk = ok = bad = unc = 0
    by_tier = collections.defaultdict(lambda: [0, 0])
    post_fork_used = 0
    bad_rows = []
    seen_fam_pairs = {}

    for gname, plist in sorted(preds_by_gene.items()):
        fam = gname.split("_")[0]
        rlist = list(fam_reads.get(fam, {}).values())
        for pr in plist:
            tot += 1
            lp, rp = pr["left_prox"], pr["right_prox"]
            if len(lp) < W or len(rp) < W:
                unc += 1
                continue
            if pr.get("depth_source") == "post_fork_bins":
                post_fork_used += 1
            n_both = n_l = n_r = 0
            for r in rlist:
                if len(r["left"]) < W or len(r["right"]) < W:
                    continue
                ml = match_side_hifi(lp, r["left"], "left", W)
                mr = match_side_hifi(rp, r["right"], "right", W)
                if ml and mr:
                    n_both += 1
                elif ml:
                    n_l += 1
                elif mr:
                    n_r += 1
            if n_both >= 1:
                verdict = "confirmed"
                chk += 1; ok += 1
                by_tier[pr["tier"]][1] += 1
            elif (n_l + n_r) >= 3 and (n_l >= 2 or n_r >= 2):
                verdict = "chimera"
                chk += 1; bad += 1
                by_tier[pr["tier"]][0] += 1
                by_tier[pr["tier"]][1] += 1
                bad_rows.append({
                    "gene": gname, "pair_id": pr["id"], "tier": pr["tier"],
                    "dL": pr.get("left_depth", 0), "dR": pr.get("right_depth", 0),
                    "n_l_only": n_l, "n_r_only": n_r, "note": pr.get("note", ""),
                })
            else:
                verdict = "uncheckable"
                unc += 1
            sig = (fam, lp[-120:], rp[:120])
            if verdict in ("confirmed", "chimera"):
                if sig not in seen_fam_pairs or verdict == "confirmed":
                    seen_fam_pairs[sig] = verdict

    fam_chk = len(seen_fam_pairs)
    fam_bad = sum(1 for v in seen_fam_pairs.values() if v == "chimera")

    # Molecule-cluster recall (>= 2 HiFi reads per cluster)
    true_all = true_hit = true_multi = true_multi_hit = 0
    fam_genes = collections.defaultdict(list)
    for g in doc_genes:
        fam_genes[g.split("_")[0]].append(g)
    for fam, glist in sorted(fam_genes.items()):
        rlist = [r for r in fam_reads.get(fam, {}).values() if len(r["left"]) >= W and len(r["right"]) >= W]
        if not rlist:
            continue
        clusters = []
        for r in rlist:
            lw, rw = r["left"][-W:], r["right"][:W]
            placed = False
            for cl in clusters:
                if match_side_hifi(lw, cl[0], "left", W) and match_side_hifi(rw, cl[1], "right", W):
                    cl[2] += 1
                    placed = True
                    break
            if not placed:
                clusters.append([lw, rw, 1])
        solid = [c for c in clusters if c[2] >= 2]
        if not solid:
            continue
        is_multi = len(solid) >= 2
        fam_preds = [pr for g in glist for pr in preds_by_gene.get(g, []) if len(pr["left_prox"]) >= W and len(pr["right_prox"]) >= W]
        for cl in solid:
            true_all += 1
            if is_multi:
                true_multi += 1
            if any(match_side_hifi(pr["left_prox"], cl[0], "left", W) and match_side_hifi(pr["right_prox"], cl[1], "right", W) for pr in fam_preds):
                true_hit += 1
                if is_multi:
                    true_multi_hit += 1

    lo, hi = cp_ci(bad, chk)
    return {
        "W": W, "emitted": tot, "checkable": chk, "confirmed": ok, "chimeras": bad, "uncheckable": unc,
        "chimera_rate_pct": round(100.0 * bad / chk, 2) if chk else 0.0,
        "ci95_pct": [round(100.0 * lo, 2), round(100.0 * hi, 2)],
        "fam_checkable": fam_chk, "fam_chimeras": fam_bad,
        "true_all": true_all, "true_hit": true_hit,
        "recall_pct": round(100.0 * true_hit / true_all, 1) if true_all else 0.0,
        "true_multi": true_multi, "true_multi_hit": true_multi_hit,
        "multi_recall_pct": round(100.0 * true_multi_hit / true_multi, 1) if true_multi else 0.0,
        "by_tier": {k: {"bad": v[0], "chk": v[1]} for k, v in sorted(by_tier.items())},
        "post_fork_bins_pairs": post_fork_used,
        "bad_details": bad_rows,
    }


def make_scorer_json(flanks_json, out_json, use_end=False):
    """Convert v0.5.0 flanks.json to claims_lr3.py schema (adding trusted_len=0)."""
    doc = json.load(open(flanks_json))
    for g, ge in doc["genes"].items():
        for sd in ("left", "right"):
            for p in ge.get(sd, {}).get("paths", []):
                p["trusted_len"] = 0
                if use_end:
                    p["context_len"] = max(int(p.get("context_len", 0)), int(p.get("end_len", 0)))
    json.dump(doc, open(out_json, "w"))


def triage_claims3_exact(genes_fa, scorer_json, paf_path, claims3_tsv, out_tsv):
    """Exact frozen `real_triage.py` classification (gene_copy / structure / rare / fault) reusing `paf_path`."""
    genes = c3.read_fa(genes_fa)
    doc = json.load(open(scorer_json))
    info = {n: v for n, *v in c3.cands(genes, doc)}
    by_key = {(g, side, str(ci)): n for n, (pad, fl, side, p, g, ci) in info.items()}
    failed_rows = []
    with open(claims3_tsv) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if r["verdict"] in ("detail", "contradicted"):
                failed_rows.append(r)
    if not failed_rows:
        with open(out_tsv, "w") as f:
            f.write("gene\tside\tcand\tclaim\tverdict\tclass\tnote\n")
        return

    want = {by_key[(r["gene"], r["side"], r["cand"])]: r for r in failed_rows if (r["gene"], r["side"], r["cand"]) in by_key}
    alns = c3.parse(paf_path, {n: info[n] for n in want})

    def gene_ok(z, pad):
        if z.o0 > -100:
            return False
        x0 = max(z.o0, -pad)
        return z.bases(x0, 0) <= 0.10 * (0 - x0)

    def holds_from0(z, x, y, fl):
        thr, div, sub = c3.MODES["strict"]
        x = max(x, 0)
        if not (z.o0 <= x):
            return False
        if not (z.o1 >= y or (y >= fl - c3.EDGE and z.o1 >= fl - c3.EDGE)):
            return False
        if z.long_inside(thr, x, y) or z.bases(x, y) > div * max(1, y - x):
            return False
        return not (sub is not None and z.hot_inside(sub, x, y))

    with open(out_tsv, "w") as f:
        f.write("gene\tside\tcand\tclaim\tverdict\tclass\tnote\n")
        for n, r in want.items():
            pad, fl, side, p, g, ci = info[n]
            al = alns.get(n, [])
            C = int(r["claim"])
            anch = [z for z in al if gene_ok(z, pad)]
            W = 3000
            wins = [(0, e) for e in range(400, min(C, W), 100)] + [(0, min(C, W))] + [(x, x + W) for x in range(100, C - W + 1, 100)]
            if wins[-1][1] < C:
                wins.append((max(0, C - W), C))
            first = [z for z in anch if holds_from0(z, 0, min(C, 400), fl)]
            okw = bool(first) and all(any(holds_from0(z, w[0], w[1], fl) for z in (anch if w[0] == 0 else al)) for w in wins)
            if okw:
                cls, note = "gene_copy", f"{len({z.rid for z in first})} gene-carrying reads hold from boundary"
            else:
                st = [z for z in anch if z.o1 >= min(C, fl - c3.EDGE) and not z.long_inside(50, 0, C) and not z.hot_inside(75, 0, C)]
                if st:
                    best = min(st, key=lambda z: z.bases(0, C))
                    cls, note = "structure", f"{len({z.rid for z in st})} reads hold order, min diffs={best.bases(0, C)}/{C}"
                else:
                    stop = int(r.get("tiled_strict") or 0)
                    beyond = {z.rid for z in anch if z.o1 >= stop + 1000 and not z.long_inside(50, 0, min(z.o1, C)) and z.bases(0, min(z.o1, C)) <= 0.03 * min(z.o1, C)}
                    if len(beyond) >= 2:
                        cls, note = "rare", f"{len(beyond)} reads follow 1kb past {stop}"
                    else:
                        cls, note = "fault", f"gene-carrying reads={len({z.rid for z in anch})}"
            f.write(f"{g}\t{side}\t{ci}\t{C}\t{r['verdict']}\t{cls}\t{note}\n")


def parse_claims3(tsv_path, triage_path=None):
    triage_map = {}
    if triage_path and os.path.isfile(triage_path):
        with open(triage_path) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= len(hdr):
                    r = dict(zip(hdr, p))
                    triage_map[(r["gene"], r["side"], r["cand"])] = r.get("class", r.get("triage", ""))
    rows = []
    with open(tsv_path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            r["claim"] = int(r["claim"])
            r["triage_cls"] = triage_map.get((r["gene"], r["side"], r["cand"]), "")
            rows.append(r)

    out = {}
    for lim_name, lo_bp, hi_bp in [(">=1000", 1000, 10**9), (">=2000", 2000, 10**9), (">=5000", 5000, 10**9), ("100-999", 100, 999)]:
        sub = [r for r in rows if lo_bp <= r["claim"] <= hi_bp]
        chk = [r for r in sub if r["verdict"] in ("spanned", "tiled", "detail", "contradicted")]
        unc = [r for r in sub if r["verdict"] in ("unseen", "reads_end")]
        raw_fail = [r for r in chk if r["verdict"] in ("detail", "contradicted")]
        adj_fail = [r for r in raw_fail if r["triage_cls"] != "gene_copy"]
        grps = {r["group"] for r in chk if r["group"] != "-"}
        raw_gfail = {r["group"] for r in raw_fail if r["group"] != "-"}
        adj_gfail = {r["group"] for r in adj_fail if r["group"] != "-"}
        lo, hi = cp_ci(len(adj_fail), len(chk))
        glo, ghi = cp_ci(len(adj_gfail), len(grps))
        out[lim_name] = {
            "checkable": len(chk), "uncheckable": len(unc),
            "raw_failed": len(raw_fail), "adj_failed": len(adj_fail),
            "adj_rate_pct": round(100.0 * len(adj_fail) / len(chk), 3) if chk else 0.0,
            "ci95_pct": [round(100.0 * lo, 2), round(100.0 * hi, 2)],
            "groups": len(grps), "raw_failed_groups": len(raw_gfail), "adj_failed_groups": len(adj_gfail),
            "group_ci95_pct": [round(100.0 * glo, 2), round(100.0 * ghi, 2)],
            "molecule_exact": sum(1 for r in chk if r["mol_verdict"] == "molecule"),
            "failed_rows": [
                {"gene": r["gene"], "side": r["side"], "cand": r["cand"], "claim": r["claim"],
                 "verdict": r["verdict"], "triage": r["triage_cls"]}
                for r in raw_fail
            ],
        }
    return out


def ensure_kmadb():
    os.makedirs(WORK_ROOT, exist_ok=True)
    kmadb = os.path.join("/tmp/argj_val/work", "kmadb")
    os.makedirs(os.path.dirname(kmadb), exist_ok=True)
    if not os.path.isfile(kmadb + ".name"):
        subprocess.run([KMA_INDEX, "-i", RESFINDER_FA, "-o", kmadb], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return kmadb


def run_one_real(meta, kmadb, min_free_gb):
    sid = meta["id"]
    res_dir = os.path.join(RESULTS_ROOT, sid)
    score_json = os.path.join(res_dir, "real_v05_score.json")
    flanks_json = os.path.join(res_dir, "run.flanks.json")
    genes_fa = os.path.join(res_dir, "genes.fa")
    hifi_gz = os.path.join(res_dir, "hifi_near.fa.gz")
    truth_json = os.path.join(res_dir, "truth_2side.json")

    if os.path.isfile(score_json) and os.path.isfile(flanks_json) and os.path.isfile(truth_json):
        print(f"[{sid}] already done, skipping", flush=True)
        return res_dir

    os.makedirs(res_dir, exist_ok=True)
    tmp_s = os.path.join(WORK_ROOT, sid)
    shutil.rmtree(tmp_s, ignore_errors=True)
    os.makedirs(tmp_s, exist_ok=True)

    t0 = time.time()
    print(f"\n================ [REAL {sid}: {meta['desc']}] ================", flush=True)
    check_disk(min_free_gb, f"before {sid}")

    sra_env = dict(os.environ)
    sra_env["HOME"] = ensure_sra_nocache()
    sra_env["PATH"] = f"/opt/homebrew/bin:{os.path.dirname(MINIMAP2)}:{os.path.dirname(SAMTOOLS)}:" + sra_env.get("PATH", "")

    run_pre = os.path.join(res_dir, "run")
    if os.path.isfile(flanks_json) and os.path.isfile(genes_fa):
        n_genes = sum(1 for ln in open(genes_fa) if ln.startswith(">"))
        print(f"  reusing completed ARGjacent outputs ({n_genes} genes)", flush=True)
    else:
        # 1. Fetch verified .sra from NCBI S3 public bucket (multipart ETag + vdb-validate),
        #    stream fastq-dump | fastp -> R1/R2.fq.gz, and delete .sra immediately!
        t_dl = time.time()
        subprocess.run(["sh", GET_RUN, meta["sr"], tmp_s], env=sra_env, check=True)
        sra_file = os.path.join(tmp_s, f"{meta['sr']}.sra")
        print(f"  verified {os.path.getsize(sra_file)/1e9:.2f} GB .sra from NCBI S3 in {time.time()-t_dl:.1f}s", flush=True)

        r1, r2 = os.path.join(tmp_s, "R1.fq.gz"), os.path.join(tmp_s, "R2.fq.gz")
        t_fp = time.time()
        cmd_fp = (
            f"set -o pipefail; '{FASTQ_DUMP}' --split-spot --skip-technical -Z '{sra_file}' 2>'{tmp_s}/fqd.log' "
            f"| '{FASTP}' --stdin --interleaved_in -o '{r1}' -O '{r2}' --detect_adapter_for_pe -w 8 -z 4 "
            f"-j '{res_dir}/fastp.json' -h /dev/null >'{tmp_s}/fastp.log' 2>&1"
        )
        subprocess.run(["bash", "-c", cmd_fp], env=sra_env, check=True)
        os.remove(sra_file)
        print(f"  fastq-dump | fastp -> {(os.path.getsize(r1)+os.path.getsize(r2))/1e9:.2f} GB in {time.time()-t_fp:.1f}s ({free_gb('/tmp'):.2f} GiB free)", flush=True)

        # 2. KMA vs ResFinder -> genes.fa
        kma_pre = os.path.join(tmp_s, "kma")
        subprocess.run(
            [KMA, "-ipe", r1, r2, "-o", kma_pre, "-t_db", kmadb, "-1t1", "-cge", "-ef", "-nf", "-t", "8"],
            stdout=subprocess.DEVNULL, stderr=open(os.path.join(res_dir, "kma.log"), "w"), check=True,
        )
        subprocess.run([sys.executable, KMA_GENES, kma_pre, genes_fa], stdout=open(os.path.join(res_dir, "kma_genes.txt"), "w"), check=True)
        shutil.copy2(kma_pre + ".res", os.path.join(res_dir, "kma.res"))
        n_genes = sum(1 for ln in open(genes_fa) if ln.startswith(">"))
        print(f"  KMA found {n_genes} ResFinder seed genes", flush=True)
        if n_genes == 0:
            shutil.rmtree(tmp_s, ignore_errors=True)
            print(f"  [{sid}] 0 seed genes; skipping", flush=True)
            return None

        # 3. Run frozen ARGjacent v0.5.0, then IMMEDIATELY delete trimmed FASTQs!
        t_aw = time.time()
        rc = subprocess.call([ARGJACENT, "-g", genes_fa, "-1", r1, "-2", r2, "-o", run_pre, "-D", "10000", "-t", "8", "--max-mem", "7", "--keep"])
        os.remove(r1); os.remove(r2)
        print(f"  ARGjacent v0.5.0 finished (rc={rc}) in {time.time()-t_aw:.1f}s ({free_gb('/tmp'):.2f} GiB free)", flush=True)
        if rc != 0:
            shutil.rmtree(tmp_s, ignore_errors=True)
            raise RuntimeError(f"argjacent failed on {sid}: rc={rc}")

    # 4. Stream PacBio HiFi directly from NCBI S3 in 6 parallel slices (zero cache, exact spot-count check!)
    ref_fa = os.path.join(tmp_s, "hifi_ref.fa")
    ref_mmi = os.path.join(tmp_s, "hifi_ref.mmi")
    p4_json = os.path.join(tmp_s, "p4.flanks.json")
    p4e_json = os.path.join(tmp_s, "p4e.flanks.json")
    make_scorer_json(flanks_json, p4_json, use_end=False)
    make_scorer_json(flanks_json, p4e_json, use_end=True)

    cand_ref_fa = os.path.join(tmp_s, "cands_ref.fa")
    subprocess.run([sys.executable, TILE_LR, "refs", genes_fa, p4_json, cand_ref_fa], check=True)
    with open(ref_fa, "w") as out_f:
        out_f.write(open(cand_ref_fa).read())
        for idx, pth in enumerate((RESFINDER_FA, genes_fa)):
            tag = "rf" if idx == 0 else "sg"
            n = 0
            for ln in open(pth):
                if ln.startswith(">"):
                    n += 1
                    out_f.write(f">{tag}{n}\n")
                else:
                    out_f.write(ln)
    subprocess.run([MINIMAP2, "-x", "map-hifi", "-d", ref_mmi, ref_fa], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    t_hf = time.time()
    nspot = int(meta["lr_spots"])
    sl = 6
    per = (nspot + sl - 1) // sl
    procs = []
    for i in range(sl):
        a = i * per + 1
        b = min(nspot, (i + 1) * per)
        if a > b:
            continue
        cmd = (
            f"set -o pipefail; '{FASTQ_DUMP}' --fasta 0 -Z -N {a} -X {b} '{meta['lr']}' 2>'{tmp_s}/dump_{i}.log' "
            f"| '{MINIMAP2}' -a -x map-hifi --secondary=no --sam-hit-only -t 1 -K 40M '{ref_mmi}' - 2>'{tmp_s}/mm2_{i}.log' "
            f"| '{SAMTOOLS}' view -h -e 'rlen >= 200' - 2>/dev/null "
            f"| '{SAMTOOLS}' fasta -F 0x900 - 2>/dev/null > '{tmp_s}/near_{i}.fa'"
        )
        procs.append((i, subprocess.Popen(["bash", "-c", cmd], env=sra_env)))
    for i, pr in procs:
        if pr.wait() != 0:
            raise RuntimeError(f"HiFi slice {i} failed for {meta['lr']}")
    # Verify exact total spot count across all slices
    tot_spots = 0
    for i, _ in procs:
        for ln in open(f"{tmp_s}/dump_{i}.log"):
            if ln.startswith("Written "):
                tot_spots += int(ln.split()[1])
    if tot_spots != nspot:
        raise RuntimeError(f"HiFi spot count mismatch for {meta['lr']}: got {tot_spots}, expected {nspot}")
    hifi_raw_fa = os.path.join(tmp_s, "hifi_near.fa")
    with open(hifi_raw_fa, "wb") as out_hf:
        for i, _ in procs:
            out_hf.write(open(f"{tmp_s}/near_{i}.fa", "rb").read())
            os.remove(f"{tmp_s}/near_{i}.fa")
    n_lr = sum(1 for ln in open(hifi_raw_fa) if ln.startswith(">"))
    subprocess.run(f"gzip -1 -c '{hifi_raw_fa}' > '{hifi_gz}'", shell=True, check=True)
    print(f"  streamed {tot_spots} HiFi spots -> kept {n_lr} near-gene HiFi reads ({os.path.getsize(hifi_gz)/1e6:.1f} MB) in {time.time()-t_hf:.1f}s", flush=True)

    # 5. Exact frozen claims_lr3.py v3.1 + real_triage.py on primary (`certain_bp`) and end-anchored (`certain_end_bp`)
    paf_path = os.path.join(tmp_s, "cands_vs_hifi.paf")
    subprocess.run(
        [MINIMAP2, "-c", "-x", "map-hifi", "--eqx", "-N", "1000", "-p", "0.001", "-t", "8", cand_ref_fa, hifi_raw_fa],
        stdout=open(paf_path, "w"), stderr=subprocess.DEVNULL, check=True,
    )
    p4_pre = os.path.join(res_dir, "lr3_primary")
    p4e_pre = os.path.join(res_dir, "lr3_end")
    subprocess.run([sys.executable, CLAIMS_LR3, genes_fa, p4_json, paf_path, p4_pre], check=True)
    subprocess.run([sys.executable, CLAIMS_LR3, genes_fa, p4e_json, paf_path, p4e_pre], check=True)
    triage_p4 = p4_pre + ".triage.tsv"
    triage_p4e = p4e_pre + ".triage.tsv"
    triage_claims3_exact(genes_fa, p4_json, paf_path, p4_pre + ".claims3.tsv", triage_p4)
    triage_claims3_exact(genes_fa, p4e_json, paf_path, p4e_pre + ".claims3.tsv", triage_p4e)

    one_sided_primary = parse_claims3(p4_pre + ".claims3.tsv", triage_p4)
    one_sided_end = parse_claims3(p4e_pre + ".claims3.tsv", triage_p4e)

    # 6. GFA v1.0 lossless verification
    genes_dict = {}
    nm, buf = None, []
    for ln in open(genes_fa):
        if ln.startswith(">"):
            if nm: genes_dict[nm] = "".join(buf).upper()
            nm, buf = ln[1:].split()[0], []
        else: buf.append(ln.strip())
    if nm: genes_dict[nm] = "".join(buf).upper()
    doc = json.load(open(flanks_json))
    doc_genes = doc["genes"]
    gfa_res = verify_gfa_file(run_pre + ".gfa", genes_dict, doc_genes)

    # 7. Real single-molecule two-sided HiFi truth (`truth_2side.py`) & pairing scoring (`W=500` & `W=1,000`)
    obs_map = extract_two_sided_truth(genes_fa, hifi_gz, mm2_bin=MINIMAP2, threads=8, min_w=500, max_w=2500)
    json.dump(obs_map, open(truth_json, "w"))
    reads_src = os.path.join(run_pre + ".tmp", "aw.reads.tsv")
    if not os.path.isfile(reads_src):
        reads_src = None

    pairing_res = {}
    for W in (500, 1000):
        phased = index_sample_reads(genes_dict, doc_genes, reads_src, W=W) if reads_src else {}
        cert_preds = {}
        for gname, gentry in doc_genes.items():
            gseq_s = genes_dict.get(gname)
            if not gseq_s: continue
            l_tr = {t["trunk_id"]: t for t in extract_trunks(gentry.get("left", {}).get("paths", []), "left", W=W)}
            r_tr = {t["trunk_id"]: t for t in extract_trunks(gentry.get("right", {}).get("paths", []), "right", W=W)}
            if W == 500 and "pairs" in gentry:
                plist = []
                for pr in gentry["pairs"]:
                    lt, rt = l_tr.get(pr["left_trunk"]), r_tr.get(pr["right_trunk"])
                    if lt and rt:
                        item = dict(pr)
                        item["left_prox"], item["right_prox"] = lt["prox"], rt["prox"]
                        item["depth_source"] = "post_fork_bins" if (lt["depth_source"] == "post_fork_bins" or rt["depth_source"] == "post_fork_bins") else "boundary"
                        plist.append(item)
                if plist: cert_preds[gname] = plist
            else:
                prs = pair_gene_single_sample(gname, gseq_s, gentry, W=W, pre_phased=phased.get(gname, {}), include_prox=True)
                for pr in prs:
                    lt, rt = l_tr.get(pr["left_trunk"]), r_tr.get(pr["right_trunk"])
                    pr["depth_source"] = "post_fork_bins" if ((lt and lt["depth_source"] == "post_fork_bins") or (rt and rt["depth_source"] == "post_fork_bins")) else "boundary"
                if prs: cert_preds[gname] = prs

        w_out = {"argjacent_v05": score_pairing_hifi(doc_genes, obs_map, cert_preds, W)}
        for bname in ("cartesian", "main", "greedy"):
            bpreds = {g: naive_pairs(ge, W, bname) for g, ge in doc_genes.items()}
            w_out[bname] = score_pairing_hifi(doc_genes, obs_map, bpreds, W)
        pairing_res[f"W{W}"] = w_out

    # Clean up temporary directories immediately!
    shutil.rmtree(run_pre + ".tmp", ignore_errors=True)
    shutil.rmtree(tmp_s, ignore_errors=True)
    if os.path.exists(run_pre + ".flanks.fa") and os.path.getsize(run_pre + ".flanks.fa") > 5_000_000:
        os.remove(run_pre + ".flanks.fa")

    summary = {
        "id": sid, "group": meta["group"], "project": meta["project"], "desc": meta["desc"],
        "sr": meta["sr"], "lr": meta["lr"], "n_genes": n_genes, "n_hifi_near": n_lr,
        "one_sided_primary": one_sided_primary,
        "one_sided_end": one_sided_end,
        "gfa": gfa_res,
        "pairing": pairing_res,
    }
    json.dump(summary, open(score_json, "w"), indent=2)

    p1 = one_sided_primary[">=1000"]
    e1 = one_sided_end[">=1000"]
    pr5 = pairing_res["W500"]["argjacent_v05"]
    pc5 = pairing_res["W500"]["cartesian"]
    print(f"=== REAL {sid} DONE in {time.time()-t0:.1f}s ({free_gb('/tmp'):.2f} GiB free) ===")
    print(f"  1-sided primary >=1kb: adj {p1['adj_failed']}/{p1['checkable']} (raw {p1['raw_failed']}), groups {p1['adj_failed_groups']}/{p1['groups']}, unc={p1['uncheckable']} | end >=1kb: {e1['adj_failed']}/{e1['checkable']}")
    print(f"  GFA v1.0 lossless    : {gfa_res.get('paths',0)} paths, {gfa_res.get('compression_ratio',0)}x compact, mismatches={gfa_res.get('roundtrip_mismatches','?')}")
    print(f"  2-sided HiFi (W=500) : ARGjacent {pr5['chimeras']}/{pr5['checkable']} chimeras ({pr5['recall_pct']}% recall, {pr5['multi_recall_pct']}% multi, post_fork_bins={pr5['post_fork_bins_pairs']}) vs Cartesian {pc5['chimeras']}/{pc5['checkable']}")
    return res_dir


def score_real_cohort(sample_metas, label, out_prefix):
    """Run shipped `argjacent-cohort` (JSON-only mode) on a real cohort and score against real HiFi!"""
    valid = []
    for m in sample_metas:
        rd = os.path.join(RESULTS_ROOT, m["id"])
        if os.path.isfile(os.path.join(rd, "run.flanks.json")) and os.path.isfile(os.path.join(rd, "truth_2side.json")):
            valid.append((m["id"], rd))
    if len(valid) < 2:
        return
    sheet = out_prefix + ".sample_sheet.tsv"
    with open(sheet, "w") as f:
        f.write("sample\tflanks_json\tgenes_fa\n")
        for sid, rd in valid:
            f.write(f"{sid}\t{rd}/run.flanks.json\t{rd}/genes.fa\n")

    for W in (500, 1000):
        pre = f"{out_prefix}_W{W}"
        subprocess.run([ARGJACENT_COHORT, "--sample-sheet", sheet, "-o", pre, "--pair-bp", str(W)], check=True)
        cdoc = json.load(open(pre + ".cohort.json"))
        tot_chk = tot_bad = tot_ok = tot_unc = tot_true = tot_hit = tot_multi = tot_multi_hit = 0
        cohort_chk = cohort_bad = 0
        by_tier = collections.defaultdict(lambda: [0, 0])
        bad_pairs = []

        for sid, rd in valid:
            doc_genes = json.load(open(os.path.join(rd, "run.flanks.json")))["genes"]
            obs_map = json.load(open(os.path.join(rd, "truth_2side.json")))
            c_preds = {}
            for gname, plist in cdoc.get("pairs", {}).get(sid, {}).items():
                ge = doc_genes.get(gname, {})
                l_tr = {t["trunk_id"]: t for t in extract_trunks(ge.get("left", {}).get("paths", []), "left", W=W)}
                r_tr = {t["trunk_id"]: t for t in extract_trunks(ge.get("right", {}).get("paths", []), "right", W=W)}
                out_l = []
                for pr in plist:
                    lt, rt = l_tr.get(pr["left_trunk"]), r_tr.get(pr["right_trunk"])
                    if lt and rt:
                        item = dict(pr)
                        item["left_prox"], item["right_prox"] = lt["prox"], rt["prox"]
                        out_l.append(item)
                if out_l: c_preds[gname] = out_l
            ps = score_pairing_hifi(doc_genes, obs_map, c_preds, W)
            tot_chk += ps["checkable"]; tot_bad += ps["chimeras"]; tot_ok += ps["confirmed"]; tot_unc += ps["uncheckable"]
            tot_true += ps["true_all"]; tot_hit += ps["true_hit"]; tot_multi += ps["true_multi"]; tot_multi_hit += ps["true_multi_hit"]
            for tn, td in ps["by_tier"].items():
                by_tier[tn][0] += td["bad"]; by_tier[tn][1] += td["chk"]
                if tn.startswith("cohort_"):
                    cohort_bad += td["bad"]; cohort_chk += td["chk"]
            for bd in ps["bad_details"]:
                bad_pairs.append(dict(bd, sample=sid))

        res = {
            "label": f"{label}_W{W}", "n_samples": len(valid), "W": W,
            "checkable": tot_chk, "confirmed": tot_ok, "chimeras": tot_bad, "uncheckable": tot_unc,
            "cohort_only_chimeras": cohort_bad, "cohort_only_checkable": cohort_chk,
            "true_all": tot_true, "true_hit": tot_hit, "recall_pct": round(100.0 * tot_hit / max(1, tot_true), 1),
            "true_multi": tot_multi, "true_multi_hit": tot_multi_hit, "multi_recall_pct": round(100.0 * tot_multi_hit / max(1, tot_multi), 1),
            "by_tier": {k: {"bad": v[0], "chk": v[1]} for k, v in sorted(by_tier.items())},
            "bad_pairs": bad_pairs,
            "n_partB_repeats": sum(len(exts) for smap in cdoc.get("repeats", {}).values() for gmap in smap.values() for exts in gmap.values()),
        }
        json.dump(res, open(pre + ".real_cohort_score.json", "w"), indent=2)
        print(f"=== REAL COHORT {label} (S={len(valid)}, W={W}) ===")
        print(f"  Part A pairing : {tot_bad}/{tot_chk} chimeras (cohort-added: {cohort_bad}/{cohort_chk}) | recall {tot_hit}/{tot_true} ({res['recall_pct']}%), multi {tot_multi_hit}/{tot_multi} ({res['multi_recall_pct']}%) | by_tier={res['by_tier']}")
        print(f"  Part B repeats : {res['n_partB_repeats']} exits emitted")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="", help="comma-separated sample IDs to run (default: all)")
    ap.add_argument("--min-free-gb", type=float, default=7.0)
    a = ap.parse_args()
    want = {x.strip() for x in a.only.split(",") if x.strip()}
    sel = [m for m in SAMPLES if (not want or m["id"] in want)]
    kmadb = ensure_kmadb()
    for m in sel:
        run_one_real(m, kmadb, a.min_free_gb)
    inf_metas = [m for m in sel if m["group"] == "infant_cohort"]
    if len(inf_metas) >= 2:
        score_real_cohort(inf_metas, "infant_PRJNA1365828", os.path.join(RESULTS_ROOT, "cohort_infant"))
    shutil.rmtree(WORK_ROOT, ignore_errors=True)
    check_disk(0.0, "end real_hifi suite")


if __name__ == "__main__":
    main()
