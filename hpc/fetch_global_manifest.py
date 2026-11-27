#!/usr/bin/env python3
"""Query ENA (and optionally MGnify + Logan) to build a curated, 1-to-1 paired-read + pre-made-assembly
manifest for global ARGjacent HPC benchmarking.

Each row links ONE paired-end WGS sequencing run (ERR/SRR) to ONE primary metagenome assembly (ERZ)
from ENA/MGnify (and optionally its Logan S3 contig archive), together with harmonized One Health biome
labels, geography, dates, FASTQ sizes, and direct HTTPS download URLs.

Examples:
    # Preview counts & build a stratified 2023 manifest (200 samples/biome, >= 3 Gbp Illumina PE)
    python3 hpc/fetch_global_manifest.py --year 2023 --per-biome 200 --max-per-study 50 \\
        --check-assembler --check-logan -o global_2023.tsv

    # Build a manifest across 2023-2024 for wastewater + human_gut + livestock + soil_plant
    python3 hpc/fetch_global_manifest.py --year 2023-2024 \\
        --biomes human_gut,wastewater,livestock,soil_plant,aquatic \\
        --min-gbp 3.0 --max-gbp 25.0 --min-study-samples 4 \\
        -o global_onehealth_2023_2024.tsv

    # Query specific curated BioProjects (e.g. Global Sewage)
    python3 hpc/fetch_global_manifest.py --studies PRJEB47609,PRJEB60351 -o curated_studies.tsv
"""
import argparse
import collections
import concurrent.futures
import json
import os
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml"
MGNIFY_API = "https://www.ebi.ac.uk/metagenomics/api/v1"
LOGAN_BASE = "https://s3.amazonaws.com/logan-pub/c"
UA = {"User-Agent": "argjacent-global-manifest/0.5.1"}

# Ordered One Health biome classification rules applied to combined metadata text
BIOME_RULES = [
    ("wastewater", re.compile(
        r"\b(wastewater|sewage|sludge|effluent|influent|wwtp|treatment plant|hospital.*water|drainage|sewer)\b", re.I
    )),
    ("human_gut", re.compile(
        r"\b(human gut|human feces|human stool|human fecal|infant gut|gut metagenome|feces metagenome|fecal metagenome|stool metagenome)\b", re.I
    )),
    ("livestock", re.compile(
        r"\b(pig|porcine|swine|chicken|poultry|gallus|bovine|cattle|cow|rumen|sheep|ovine|goat|caprine|manure|livestock|broiler|turkey|dairy|calf|feedlot|horse|equine)\b", re.I
    )),
    ("soil_plant", re.compile(
        r"\b(soil|rhizosphere|compost|peat|permafrost|plant|root|phyllosphere|forest|grassland|agricultural|farmland|tundra|sediment.*terrestrial|moss|lichen)\b", re.I
    )),
    ("aquatic", re.compile(
        r"\b(marine|seawater|ocean|freshwater|river|lake|estuar|aquatic|sediment|coastal|coral|reef|plankton|groundwater|pond|stream|hot spring|hydrothermal|glacier|ice|brackish|fish|aquaculture)\b", re.I
    )),
    ("human_other", re.compile(
        r"\b(human oral|human skin|human vaginal|human lung|human nasal|human nasopharyngeal|human respiratory|human urinary|urine metagenome|saliva|plaque|sputum|homo sapiens|human metagenome|airway|throat|milk metagenome)\b", re.I
    )),
]

MISSING_VALS = {"", "missing", "not applicable", "not collected", "not provided", "unknown", "na", "n/a", "none", "null"}


def http_request(url, data=None, method="GET", headers=None, timeout=60, retries=3):
    hdrs = dict(UA)
    if headers:
        hdrs.update(headers)
    body = urllib.parse.urlencode(data).encode("utf-8") if isinstance(data, dict) else data
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=body, headers=hdrs, method=method)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.status, dict(resp.headers), resp.read()
        except urllib.error.HTTPError as e:
            if e.code in (400, 404, 410):
                raise
            last_err = e
        except Exception as e:
            last_err = e
        time.sleep(1.5 * (attempt + 1))
    raise last_err


def ena_post_tsv(params, timeout=120):
    _, _, raw = http_request(ENA_SEARCH, data=params, method="POST", timeout=timeout)
    text = raw.decode("utf-8", errors="replace").strip()
    if not text:
        return []
    lines = text.splitlines()
    if len(lines) < 2:
        return []
    hdr = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < len(hdr):
            parts.extend([""] * (len(hdr) - len(parts)))
        rows.append(dict(zip(hdr, parts)))
    return rows


def to_https(ftp_path):
    p = (ftp_path or "").strip()
    if not p:
        return ""
    if p.startswith("https://") or p.startswith("http://"):
        return p
    if p.startswith("ftp://"):
        return "https://" + p[6:]
    return "https://" + p.lstrip("/")


def pick_assembly_url(row):
    """Select the best FASTA contig URL from ENA generated_ftp / submitted_ftp."""
    cands = []
    for col in ("generated_ftp", "submitted_ftp"):
        for item in (row.get(col) or "").split(";"):
            item = item.strip()
            if item:
                cands.append(item)
    # Prefer standard compressed FASTA contig archives
    for rx in (
        re.compile(r"contig.*\.f(a|na|asta)\.gz$", re.I),
        re.compile(r"\.f(a|na|asta)\.gz$", re.I),
        re.compile(r"\.f(a|na|asta)$", re.I),
    ):
        for c in cands:
            if rx.search(c):
                return to_https(c)
    return to_https(cands[0]) if cands else ""


def pick_paired_fastqs(fastq_ftp, fastq_bytes):
    """Extract strictly matched paired-end (_1, _2) FASTQ URLs and total compressed bytes."""
    urls = [u.strip() for u in (fastq_ftp or "").split(";") if u.strip()]
    bs_raw = [b.strip() for b in (fastq_bytes or "").split(";") if b.strip()]
    bmap = {}
    for i, u in enumerate(urls):
        try:
            bmap[u] = int(bs_raw[i]) if i < len(bs_raw) and bs_raw[i].isdigit() else 0
        except ValueError:
            bmap[u] = 0

    r1 = r2 = ""
    for u in urls:
        if re.search(r"_1\.f(ast)?q\.gz$", u):
            r1 = u
        elif re.search(r"_2\.f(ast)?q\.gz$", u):
            r2 = u
    if not (r1 and r2) and len(urls) == 2:
        r1, r2 = urls[0], urls[1]
    if not (r1 and r2):
        return "", "", 0
    return to_https(r1), to_https(r2), bmap.get(r1, 0) + bmap.get(r2, 0)


def classify_biome(row):
    """Map ENA/MGnify sample metadata into clean One Health biome categories."""
    # Check specific host/livestock signals before generic 'gut metagenome'
    sci = (row.get("scientific_name") or "").strip()
    env = (row.get("environment_biome") or "").strip()
    mat = (row.get("environment_material") or "").strip()
    iso = (row.get("isolation_source") or "").strip()
    host = (row.get("host") or "").strip()
    mgn = (row.get("mgnify_biome") or "").strip()
    combined = " | ".join(x for x in (mgn, sci, env, mat, iso, host) if x)

    # Disambiguate animal vs human gut first
    if re.search(r"\b(wastewater|sewage|sludge|effluent|influent|wwtp|sewer)\b", combined, re.I):
        return "wastewater"
    if re.search(r"\b(pig|porcine|swine|chicken|poultry|gallus|bovine|cattle|cow|rumen|sheep|ovine|goat|caprine|manure|livestock|broiler|turkey|dairy|calf|feedlot|sus scrofa|bos taurus)\b", combined, re.I):
        return "livestock"
    if re.search(r"\b(human gut|human feces|human stool|human fecal|infant gut)\b", combined, re.I):
        return "human_gut"
    if re.search(r"\b(gut metagenome|feces metagenome|fecal metagenome|stool metagenome)\b", combined, re.I):
        if re.search(r"\b(mouse|murine|mus musculus|rat|rattus|dog|canine|cat|feline|fish|insect|bee|termite|bird|avian|macaque|primate|wildlife|panda|bear|bat)\b", combined, re.I):
            return "other"
        if re.search(r"\b(human|homo sapiens|patient|infant|adult|child|cohort|clinical)\b", combined, re.I) or not host:
            return "human_gut"

    for name, rx in BIOME_RULES:
        if rx.search(combined):
            return name
    return "other"


def clean_country(raw):
    c = (raw or "").split(":")[0].strip()
    return "" if c.lower() in MISSING_VALS else c


def parse_year_spec(spec):
    spec = spec.strip()
    if "-" in spec:
        a, b = spec.split("-", 1)
        return int(a), int(b)
    y = int(spec)
    return y, y


def fetch_ena_analyses(year_spec=None, date_field="first_public", studies=None, limit=0):
    """Fetch primary metagenome assemblies (ERZ) from ENA Portal API."""
    clauses = ['analysis_type="SEQUENCE_ASSEMBLY"', 'assembly_type="primary metagenome"']
    if year_spec:
        y0, y1 = parse_year_spec(year_spec)
        clauses.append(f"{date_field}>={y0}-01-01")
        clauses.append(f"{date_field}<={y1}-12-31")
    if studies:
        s_or = " OR ".join(f'study_accession="{s}"' for s in studies)
        clauses.append(f"({s_or})")

    query = " AND ".join(clauses)
    fields = ",".join([
        "analysis_accession",
        "run_accession",
        "study_accession",
        "sample_accession",
        "first_public",
        "generated_ftp",
        "submitted_ftp",
        "assembly_software",
        "country",
        "environment_biome",
        "environment_material",
        "scientific_name",
    ])
    sys.stderr.write(f"[1/4] Querying ENA primary metagenome assemblies: {query}\n")
    rows = ena_post_tsv({
        "result": "analysis",
        "query": query,
        "fields": fields,
        "limit": str(limit),
        "format": "tsv",
    })
    sys.stderr.write(f"      Retrieved {len(rows):,} ENA primary metagenome assembly rows.\n")
    return rows


def fetch_mgnify_biome_erz(lineage, max_pages=50):
    """Fetch ERZ assembly accessions from MGnify API v1 for a specific biome lineage."""
    sys.stderr.write(f"[1/4] Querying MGnify API v1 assemblies for biome lineage '{lineage}'...\n")
    erz_to_biome = {}
    url = f"{MGNIFY_API}/analyses?experiment_type=assembly&lineage={urllib.parse.quote(lineage, safe='')}&page_size=100"
    page = 0
    while url and page < max_pages:
        page += 1
        _, _, raw = http_request(url, headers={"Accept": "application/json"}, timeout=60)
        doc = json.loads(raw.decode("utf-8"))
        for item in doc.get("data", []):
            assem = (item.get("relationships", {}).get("assembly", {}).get("data") or {}).get("id")
            if assem:
                erz_to_biome[assem] = lineage
        url = (doc.get("links") or {}).get("next")
    sys.stderr.write(f"      Found {len(erz_to_biome):,} unique MGnify ERZ assemblies in '{lineage}'.\n")
    return erz_to_biome


def fetch_ena_analyses_by_erz(erz_list, chunk_size=250):
    """Fetch ENA analysis metadata for an explicit list of ERZ accessions."""
    fields = ",".join([
        "analysis_accession", "run_accession", "study_accession", "sample_accession",
        "first_public", "generated_ftp", "submitted_ftp", "assembly_software",
        "country", "environment_biome", "environment_material", "scientific_name",
    ])
    out = []
    for i in range(0, len(erz_list), chunk_size):
        chunk = erz_list[i:i + chunk_size]
        q = " OR ".join(f'analysis_accession="{a}"' for a in chunk)
        out.extend(ena_post_tsv({"result": "analysis", "query": q, "fields": fields, "limit": "0", "format": "tsv"}))
    return out


def fetch_read_runs(run_accs, chunk_size=300, workers=4):
    """Batch-fetch ENA read_run metadata for a list of run accessions."""
    fields = ",".join([
        "run_accession", "study_accession", "sample_accession",
        "instrument_platform", "instrument_model",
        "library_strategy", "library_layout",
        "read_count", "base_count",
        "fastq_ftp", "fastq_bytes",
        "collection_date", "first_public", "country",
        "environment_biome", "environment_material",
        "isolation_source", "host", "scientific_name",
    ])
    chunks = [run_accs[i:i + chunk_size] for i in range(0, len(run_accs), chunk_size)]
    sys.stderr.write(f"[2/4] Batch-fetching ENA read_run metadata for {len(run_accs):,} single-run accessions ({len(chunks)} chunks)...\n")
    by_run = {}

    def _one(ch):
        q = " OR ".join(f'run_accession="{r}"' for r in ch)
        return ena_post_tsv({"result": "read_run", "query": q, "fields": fields, "limit": "0", "format": "tsv"})

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        for res in ex.map(_one, chunks):
            for r in res:
                ra = r.get("run_accession", "").strip()
                if ra:
                    by_run[ra] = r
    sys.stderr.write(f"      Matched {len(by_run):,} / {len(run_accs):,} runs in ENA read_run.\n")
    return by_run


def lookup_assembler_xml(erz):
    """Parse assembler <PROGRAM> and <COVERAGE> from ENA analysis XML."""
    try:
        _, _, raw = http_request(f"{ENA_XML}/{erz}", timeout=20, retries=2)
        root = ET.fromstring(raw)
        prog = ""
        for el in root.iter():
            tag = el.tag.split("}")[-1].upper()
            if tag == "PROGRAM" and el.text and el.text.strip():
                prog = el.text.strip()
                break
        return erz, prog
    except Exception:
        return erz, ""


def lookup_logan_head(run_acc):
    """Check whether Logan S3 hosts a pre-built contig archive for run_acc."""
    url = f"{LOGAN_BASE}/{run_acc}/{run_acc}.contigs.fa.zst"
    try:
        status, hdrs, _ = http_request(url, method="HEAD", timeout=15, retries=2)
        if status == 200:
            cl = int(hdrs.get("Content-Length") or hdrs.get("content-length") or 0)
            return run_acc, url, cl
    except Exception:
        pass
    return run_acc, "", 0


def stratify_records(records, per_biome=0, max_per_study=0, min_study_samples=1, seed=42):
    """Apply study-level cohort thresholds and deterministic stratified subsampling across biomes."""
    rng = random.Random(seed)

    # Filter out studies with fewer than min_study_samples if requested
    if min_study_samples > 1:
        scnt = collections.Counter(r["study_accession"] for r in records)
        records = [r for r in records if scnt[r["study_accession"]] >= min_study_samples]

    # Group by biome, then cap per study and sample evenly across studies inside each biome
    by_biome = collections.defaultdict(list)
    for r in records:
        by_biome[r["biome_group"]].append(r)

    selected = []
    for biome in sorted(by_biome):
        items = list(by_biome[biome])
        # Group by study so we draw evenly across studies (preventing 1 huge BioProject dominating)
        by_study = collections.defaultdict(list)
        for r in items:
            by_study[r["study_accession"]].append(r)
        for st in by_study:
            by_study[st].sort(key=lambda x: x["run_accession"])
            rng.shuffle(by_study[st])
            if max_per_study > 0:
                by_study[st] = by_study[st][:max_per_study]

        if per_biome > 0:
            # Round-robin across studies (shuffled deterministically) until per_biome reached
            study_keys = sorted(by_study)
            rng.shuffle(study_keys)
            # Prioritize studies with >= 4 samples so argjacent-cohort works well
            study_keys.sort(key=lambda k: -(len(by_study[k]) >= 4))
            picked = []
            idx = 0
            while len(picked) < per_biome:
                added_any = False
                for st in study_keys:
                    if idx < len(by_study[st]):
                        picked.append(by_study[st][idx])
                        added_any = True
                        if len(picked) >= per_biome:
                            break
                if not added_any:
                    break
                idx += 1
            selected.extend(picked)
        else:
            for st in sorted(by_study):
                selected.extend(by_study[st])

    selected.sort(key=lambda r: (r["biome_group"], r["study_accession"], r["run_accession"]))
    return selected


def main():
    ap = argparse.ArgumentParser(
        description="Build a curated global metagenome manifest linking 1-to-1 ENA raw paired reads with pre-made ERZ (and optional Logan) assemblies."
    )
    src = ap.add_argument_group("Query source (provide --year, --studies, and/or --mgnify-biome)")
    src.add_argument("--year", help="Release/collection year or range, e.g. '2023' or '2023-2024'")
    src.add_argument("--date-field", choices=("first_public", "collection_date"), default="first_public",
                     help="ENA date field to filter --year against (default: first_public)")
    src.add_argument("--studies", help="Comma-separated ENA/SRA BioProject accessions (or path to file with one per line)")
    src.add_argument("--mgnify-biome", help="MGnify biome lineage (e.g. 'root:Engineered:Wastewater') to query assemblies from")
    src.add_argument("--mgnify-max-pages", type=int, default=30, help="Max MGnify pages of 100 analyses to fetch (default: 30)")
    src.add_argument("--ena-limit", type=int, default=0, help="Optional raw ENA analysis limit for fast testing (0 = all)")

    flt = ap.add_argument_group("Quality & metadata filters")
    flt.add_argument("--platforms", default="ILLUMINA",
                     help="Comma-separated allowed sequencing platforms (default: ILLUMINA; e.g. ILLUMINA,DNBSEQ)")
    flt.add_argument("--min-gbp", type=float, default=3.0,
                     help="Minimum sequenced bases in Gbp (~10M 2x150 pairs = 3.0 Gbp; default: 3.0)")
    flt.add_argument("--max-gbp", type=float, default=30.0,
                     help="Maximum sequenced bases in Gbp to avoid extreme outliers (default: 30.0; 0 = unlimited)")
    flt.add_argument("--min-bp-per-record", type=float, default=135.0,
                     help="Minimum base_count/read_count ratio to exclude ultra-short reads (default: 135.0)")
    flt.add_argument("--biomes", help="Comma-separated One Health biome buckets to keep (e.g. human_gut,wastewater,livestock,soil_plant,aquatic)")
    flt.add_argument("--require-country", action="store_true", help="Keep only samples with non-empty country metadata")
    flt.add_argument("--per-biome", type=int, default=0, help="Stratified sample target per biome bucket (0 = keep all)")
    flt.add_argument("--max-per-study", type=int, default=0, help="Cap samples per BioProject within each biome (0 = unlimited)")
    flt.add_argument("--min-study-samples", type=int, default=1, help="Require at least K passing samples per BioProject (e.g. 4 for cohort mode)")
    flt.add_argument("--seed", type=int, default=42, help="Deterministic random seed for stratified sampling (default: 42)")

    enr = ap.add_argument_group("Optional enrichment (run on final selected rows)")
    enr.add_argument("--check-assembler", action="store_true", help="Fetch ENA ERZ XML <PROGRAM> (MEGAHIT / metaSPAdes) concurrently")
    enr.add_argument("--check-logan", action="store_true", help="HTTP HEAD check AWS Open Data Logan contig availability per run")
    enr.add_argument("--workers", type=int, default=8, help="Concurrent HTTP workers for batch lookups (default: 8)")

    ap.add_argument("-o", "--out", required=True, help="Output manifest TSV path")
    args = ap.parse_args()

    if not (args.year or args.studies or args.mgnify_biome):
        ap.error("Provide at least one of --year, --studies, or --mgnify-biome.")

    studies = None
    if args.studies:
        if os.path.isfile(args.studies):
            studies = [l.strip() for l in open(args.studies) if l.strip() and not l.startswith("#")]
        else:
            studies = [s.strip() for s in args.studies.split(",") if s.strip()]

    mgn_biomes = {}
    if args.mgnify_biome:
        mgn_biomes = fetch_mgnify_biome_erz(args.mgnify_biome, max_pages=args.mgnify_max_pages)
        erz_list = sorted(mgn_biomes)
        if args.ena_limit > 0:
            erz_list = erz_list[:args.ena_limit]
        arows = fetch_ena_analyses_by_erz(erz_list)
    else:
        arows = fetch_ena_analyses(
            year_spec=args.year, date_field=args.date_field, studies=studies, limit=args.ena_limit
        )

    # Filter strictly to 1-to-1 single-run assemblies with valid FASTA URLs
    single_by_run = {}
    for r in arows:
        ra = (r.get("run_accession") or "").strip()
        erz = (r.get("analysis_accession") or "").strip()
        if not ra or not erz or ";" in ra or "," in ra or " " in ra:
            continue
        asm_url = pick_assembly_url(r)
        if not asm_url:
            continue
        # Keep the newest analysis per run_accession
        prev = single_by_run.get(ra)
        if prev is None or (r.get("first_public") or "") > (prev.get("first_public") or ""):
            r = dict(r)
            r["assembly_url"] = asm_url
            if erz in mgn_biomes:
                r["mgnify_biome"] = mgn_biomes[erz]
            single_by_run[ra] = r

    sys.stderr.write(f"      1-to-1 single-run assemblies with FASTA URLs: {len(single_by_run):,}\n")
    if not single_by_run:
        sys.exit("No 1-to-1 single-run assemblies matched the query.")

    run_meta = fetch_read_runs(sorted(single_by_run), workers=min(args.workers, 6))

    allowed_platforms = {p.strip().upper() for p in args.platforms.split(",") if p.strip()}
    allowed_biomes = {b.strip().lower() for b in args.biomes.split(",") if b.strip()} if args.biomes else None

    sys.stderr.write("[3/4] Filtering for paired-end WGS reads & classifying One Health biomes...\n")
    passing = []
    drop_reasons = collections.Counter()
    for ra, ar in single_by_run.items():
        rr = run_meta.get(ra)
        if not rr:
            drop_reasons["missing_read_run"] += 1
            continue
        plat = (rr.get("instrument_platform") or "").strip().upper()
        if allowed_platforms and plat not in allowed_platforms:
            drop_reasons[f"platform:{plat or 'blank'}"] += 1
            continue
        if (rr.get("library_strategy") or "").strip().upper() != "WGS":
            drop_reasons["non_WGS"] += 1
            continue
        if (rr.get("library_layout") or "").strip().upper() != "PAIRED":
            drop_reasons["non_PAIRED"] += 1
            continue
        r1_url, r2_url, fq_bytes = pick_paired_fastqs(rr.get("fastq_ftp"), rr.get("fastq_bytes"))
        if not (r1_url and r2_url):
            drop_reasons["missing_paired_fastq_urls"] += 1
            continue
        try:
            rc = int(rr.get("read_count") or 0)
            bc = int(rr.get("base_count") or 0)
        except ValueError:
            drop_reasons["bad_counts"] += 1
            continue
        gbp = bc / 1e9
        if gbp < args.min_gbp:
            drop_reasons["below_min_gbp"] += 1
            continue
        if args.max_gbp > 0 and gbp > args.max_gbp:
            drop_reasons["above_max_gbp"] += 1
            continue
        bp_per_rec = (bc / rc) if rc > 0 else 0.0
        if bp_per_rec < args.min_bp_per_record:
            drop_reasons["short_reads"] += 1
            continue

        merged = {
            "sample": ra,
            "run_accession": ra,
            "erz_accession": ar["analysis_accession"].strip(),
            "study_accession": (rr.get("study_accession") or ar.get("study_accession") or "").strip(),
            "sample_accession": (rr.get("sample_accession") or ar.get("sample_accession") or "").strip(),
            "scientific_name": (rr.get("scientific_name") or ar.get("scientific_name") or "").strip(),
            "environment_biome": (rr.get("environment_biome") or ar.get("environment_biome") or "").strip(),
            "environment_material": (rr.get("environment_material") or ar.get("environment_material") or "").strip(),
            "isolation_source": (rr.get("isolation_source") or "").strip(),
            "host": (rr.get("host") or "").strip(),
            "mgnify_biome": ar.get("mgnify_biome", ""),
            "country": clean_country(rr.get("country") or ar.get("country")),
            "collection_date": (rr.get("collection_date") or "").strip(),
            "first_public": (ar.get("first_public") or rr.get("first_public") or "").strip(),
            "instrument_platform": plat,
            "instrument_model": (rr.get("instrument_model") or "").strip(),
            "read_count": rc,
            "base_count": bc,
            "est_gbp": f"{gbp:.2f}",
            "fastq_gb": f"{fq_bytes / (1024 ** 3):.2f}" if fq_bytes > 0 else "",
            "fastq_bytes": fq_bytes,
            "r1_url": r1_url,
            "r2_url": r2_url,
            "assembly_url": ar["assembly_url"],
            "assembler": (ar.get("assembly_software") or "").strip(),
            "logan_url": "",
            "logan_mb": "",
        }
        merged["biome_group"] = classify_biome(merged)
        if allowed_biomes and merged["biome_group"] not in allowed_biomes:
            drop_reasons[f"biome_excluded:{merged['biome_group']}"] += 1
            continue
        if args.require_country and not merged["country"]:
            drop_reasons["missing_country"] += 1
            continue
        passing.append(merged)

    sys.stderr.write(f"      Passed all read/QC filters: {len(passing):,} samples.\n")
    if drop_reasons:
        top_drops = ", ".join(f"{k}={v}" for k, v in drop_reasons.most_common(6))
        sys.stderr.write(f"      Top filter drops: {top_drops}\n")

    if not passing:
        sys.exit("No samples passed all filters. Try relaxing --min-gbp, --platforms, or --biomes.")

    selected = stratify_records(
        passing,
        per_biome=args.per_biome,
        max_per_study=args.max_per_study,
        min_study_samples=args.min_study_samples,
        seed=args.seed,
    )
    sys.stderr.write(f"[4/4] Selected {len(selected):,} samples after stratification.\n")

    if args.check_assembler and selected:
        need_xml = [r["erz_accession"] for r in selected if not r["assembler"]]
        if need_xml:
            sys.stderr.write(f"      Fetching ENA ERZ XML assembler provenance for {len(need_xml):,} assemblies...\n")
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
                amap = dict(ex.map(lookup_assembler_xml, need_xml))
            for r in selected:
                if not r["assembler"]:
                    r["assembler"] = amap.get(r["erz_accession"], "")

    if args.check_logan and selected:
        sys.stderr.write(f"      Checking AWS Open Data Logan contig availability for {len(selected):,} runs...\n")
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
            lres = list(ex.map(lookup_logan_head, [r["run_accession"] for r in selected]))
        lmap = {ra: (u, b) for ra, u, b in lres}
        for r in selected:
            u, b = lmap.get(r["run_accession"], ("", 0))
            r["logan_url"] = u
            r["logan_mb"] = f"{b / (1024 ** 2):.1f}" if b > 0 else ""

    cols = [
        "sample", "run_accession", "erz_accession", "study_accession", "sample_accession",
        "biome_group", "scientific_name", "environment_biome", "country",
        "collection_date", "first_public", "instrument_platform", "instrument_model",
        "read_count", "base_count", "est_gbp", "fastq_gb",
        "r1_url", "r2_url", "assembly_url", "assembler", "logan_url", "logan_mb",
    ]
    out_dir = os.path.dirname(os.path.abspath(args.out))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in selected:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # Print rich summary statistics to stdout
    biomes = collections.Counter(r["biome_group"] for r in selected)
    studies_c = collections.Counter(r["study_accession"] for r in selected)
    countries = collections.Counter(r["country"] for r in selected if r["country"])
    asms = collections.Counter(r["assembler"] or "(not_checked)" for r in selected)
    tot_tb = sum(r["fastq_bytes"] for r in selected) / (1024 ** 4)
    tot_gbp = sum(r["base_count"] for r in selected) / 1e9
    n_cohort_studies = sum(1 for _, c in studies_c.items() if c >= 4)
    n_logan = sum(1 for r in selected if r["logan_url"])

    print(f"\n=== Global Metagenome Manifest Written: {args.out} ===")
    print(f"  Total samples      : {len(selected):,} (across {len(studies_c):,} BioProjects, {n_cohort_studies} with >=4 samples for cohort mode)")
    print(f"  Unique countries   : {len(countries):,}" + (f" (top: {', '.join(f'{k}={v}' for k, v in countries.most_common(6))})" if countries else ""))
    print(f"  Total sequencing   : {tot_gbp:,.1f} Gbp | Transient FASTQ download: {tot_tb:.2f} TB")
    print(f"  Est. HPC compute   : ~{len(selected) * 0.09:.0f}-{len(selected) * 0.14:.0f} node-hrs @ 8 CPUs (~{max(1.0, tot_tb * 1024 / 200 / 3.6):.1f}-{max(1.0, tot_tb * 1024 / 100 / 3.6):.1f} hrs download @ 100-200 MB/s)")
    print("  Biome breakdown    :")
    for b, c in biomes.most_common():
        print(f"    - {b:<14s}: {c:,}")
    if args.check_assembler:
        print("  Assemblers (XML)   : " + ", ".join(f"{k}={v}" for k, v in asms.most_common(6)))
    if args.check_logan:
        print(f"  Logan S3 matches   : {n_logan:,} / {len(selected):,} ({100.0 * n_logan / len(selected):.1f}%)")


if __name__ == "__main__":
    main()
