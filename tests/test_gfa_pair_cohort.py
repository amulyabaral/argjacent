"""Tests for 0.5.0 GFA v1.0 export, depth_bins, single-sample pairing, and argjacent-cohort deconvolution."""
import json
import os
import random
import subprocess
import sys

from argjacent import cohort, gfa, pair

COMP = str.maketrans("ACGT", "TGCA")


def _rc(s):
    return s.translate(COMP)[::-1]


def _rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _sim_reads(rng, molecules, prefix, rlen=150, ins_mean=420, ins_sd=40, err=0.002):
    """Simulate paired-end reads from circular molecules `[(seq, cov), ...]`."""
    q = "I" * rlen
    with open(prefix + "_1.fq", "w") as f1, open(prefix + "_2.fq", "w") as f2:
        for mi, (seq, cov) in enumerate(molecules):
            if cov <= 0:
                continue
            L = len(seq)
            cc = seq + seq[:2000]
            n_pairs = int(cov * L / (2 * rlen))
            for k in range(n_pairs):
                ins = max(rlen + 20, int(rng.gauss(ins_mean, ins_sd)))
                pos = rng.randrange(0, L)
                frag = cc[pos:pos + ins]
                if rng.random() < 0.5:
                    frag = _rc(frag)
                r1 = list(frag[:rlen])
                r2 = list(_rc(frag)[:rlen])
                for arr in (r1, r2):
                    for idx in range(rlen):
                        if rng.random() < err:
                            arr[idx] = rng.choice([b for b in "ACGT" if b != arr[idx]])
                f1.write("@m%d_%d/1\n%s\n+\n%s\n" % (mi, k, "".join(r1), q))
                f2.write("@m%d_%d/2\n%s\n+\n%s\n" % (mi, k, "".join(r2), q))


def _run_cohort(args, expect=0):
    r = subprocess.run([sys.executable, "-m", "argjacent.cohort"] + list(args), capture_output=True, text=True)
    if expect is not None:
        assert r.returncode == expect, r.stdout + r.stderr
    return r.returncode, r.stdout, r.stderr


def _parse_gfa(path):
    segs, links, paths = {}, [], {}
    header = None
    for raw in open(path):
        line = raw.rstrip("\n")
        if not line:
            continue
        toks = line.split("\t")
        rec = toks[0]
        if rec == "H":
            header = line
        elif rec == "S":
            tags = {t.split(":", 2)[0]: t.split(":", 2)[2] for t in toks[3:]}
            segs[toks[1]] = {"seq": toks[2], "tags": tags}
        elif rec == "L":
            links.append(tuple(toks[1:6]))
        elif rec == "P":
            paths[toks[1]] = [x[:-1] for x in toks[2].split(",")]
    return header, segs, links, paths


def test_gfa_outputs_depth_bins_and_flags(data, cli, tmp_path):
    out = str(tmp_path / "rep")
    cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "repeat_1.fq"),
         "-2", os.path.join(data, "repeat_2.fq"), "-o", out, "-t", "2"])
    assert os.path.getsize(out + ".gfa") > 0
    assert os.path.getsize(out + ".pairs.tsv") > 0

    doc = json.load(open(out + ".flanks.json"))
    g_ent = doc["genes"]["demo_gene"]
    gene_seq = "".join(l.strip() for l in open(os.path.join(data, "gene.fa")) if not l.startswith(">"))

    # depth_bins & gene_depth reflect the 5-copy 1.3 kb repeat at 3,000-4,300 bp on the right side
    r_c1 = g_ent["right"]["paths"][0]
    assert r_c1["gene_depth"] >= 20.0
    bins = r_c1["depth_bins"]
    assert len(bins) >= 45
    uniq_med = sorted(bins[5:25])[10]
    rep_med = sorted(bins[32:41])[4]
    assert rep_med >= 3.0 * uniq_med, (uniq_med, rep_med)

    # GFA v1.0 lossless roundtrip and local DP:f: jump on repeat segment
    hdr, segs, links, paths = _parse_gfa(out + ".gfa")
    assert hdr == "H\tVN:Z:1.0"
    assert all(cig == "0M" for _, _, _, _, cig in links)
    assert segs["demo_gene|GENE"]["seq"] == gene_seq
    assert segs["demo_gene|GENE"]["tags"]["TP"] == "gene"
    assert segs["demo_gene|GENE"]["tags"]["CL"] == "#d95f02"
    for sd in ("left", "right"):
        for p in g_ent[sd]["paths"]:
            pname = "demo_gene|%s|%s" % (sd, p["id"])
            assert pname in paths
            recon = "".join(segs[s]["seq"] for s in paths[pname])
            expect_seq = (p["flank"] + gene_seq) if sd == "left" else (gene_seq + p["flank"])
            assert recon == expect_seq
    # Right side must have both certain (teal) and candidate (grey) segments, and elevated DP at the repeat
    r_dps = [float(v["tags"]["DP"]) for k, v in segs.items() if "|R" in k]
    assert max(r_dps) >= 2.0 * min(r_dps)

    # --no-gfa and --no-pair suppress those files; failed rerun removes stale .gfa and .pairs.tsv
    out2 = str(tmp_path / "clean_min")
    cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "clean_1.fq"),
         "-2", os.path.join(data, "clean_2.fq"), "-o", out2, "-t", "2", "--no-gfa", "--no-pair"])
    assert not os.path.exists(out2 + ".gfa")
    assert not os.path.exists(out2 + ".pairs.tsv")
    assert "pairs" not in json.load(open(out2 + ".flanks.json"))["genes"]["demo_gene"]

    cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "clean_1.fq"),
         "-2", os.path.join(data, "clean_2.fq"), "-o", out, "-t", "2", "--map", os.path.join(data, "gene.fa")], expect=1)
    assert not os.path.exists(out + ".gfa") and not os.path.exists(out + ".pairs.tsv")


def test_single_sample_pairing_clean_mixed_and_shared(data, cli, tmp_path):
    truth = json.load(open(os.path.join(data, "truth.json")))

    # clean -> 1x1 unbranched pair, also emitted as a P line in .gfa
    out_c = str(tmp_path / "clean")
    cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "clean_1.fq"),
         "-2", os.path.join(data, "clean_2.fq"), "-o", out_c, "-t", "2"])
    doc_c = json.load(open(out_c + ".flanks.json"))
    prs_c = doc_c["genes"]["demo_gene"].get("pairs", [])
    assert len(prs_c) == 1 and prs_c[0]["tier"] == "unbranched"
    rows_c = [l.rstrip("\n").split("\t") for l in open(out_c + ".pairs.tsv")][1:]
    assert len(rows_c) == 1 and rows_c[0][8] == "unbranched"
    _, _, _, gpaths = _parse_gfa(out_c + ".gfa")
    assert "demo_gene|pair|p1|unbranched" in gpaths

    # mixed (80x vs 50x) -> 0 chimeras (any emitted pair must match same genome A-A or B-B)
    out_m = str(tmp_path / "mixed")
    cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "mixed_1.fq"),
         "-2", os.path.join(data, "mixed_2.fq"), "-o", out_m, "-t", "2"])
    doc_m = json.load(open(out_m + ".flanks.json"))["genes"]["demo_gene"]
    l_by_id = {p["id"]: p["flank"] for p in doc_m["left"]["paths"]}
    r_by_id = {p["id"]: p["flank"] for p in doc_m["right"]["paths"]}
    for pr in doc_m.get("pairs", []):
        lf, rf = l_by_id[pr["left_id"]][-500:], r_by_id[pr["right_id"]][:500]
        ok_a = (lf == truth["mixed"]["A"]["left"][-500:] and rf == truth["mixed"]["A"]["right"][:500])
        ok_b = (lf == truth["mixed"]["B"]["left"][-500:] and rf == truth["mixed"]["B"]["right"][:500])
        assert ok_a or ok_b, "Chimeric pair emitted on mixed!"

    # shared (foreign locus shares last 160 bp of gene at 60x vs true 40x) -> never pairs foreign locus
    out_s = str(tmp_path / "shared")
    cli(["-g", os.path.join(data, "gene.fa"), "-1", os.path.join(data, "shared_1.fq"),
         "-2", os.path.join(data, "shared_2.fq"), "-o", out_s, "-t", "2"])
    doc_s = json.load(open(out_s + ".flanks.json"))["genes"]["demo_gene"]
    prs_s = doc_s.get("pairs", [])
    r_map = {p["id"]: p["flank"] for p in doc_s["right"]["paths"]}
    for pr in prs_s:
        assert r_map[pr["right_id"]][:500] == truth["shared"]["A"]["right"][:500]


def test_post_fork_depth_bins_unlocks_shared_promoter_pairing(cli, tmp_path):
    """Two plasmids (180x vs 50x) share an identical 250-bp cassette on BOTH flanks of a 1,100-bp gene.
    Boundary depth (0-200 bp) is summed (~230x) on both trunks, so legacy 0.4.1 JSON must abstain;
    0.5.0 depth_bins measures exclusive post-fork depth (~180x vs ~50x) and resolves both pairs with 0 chimeras.
    Also verifies that at 120x vs 30x (where one rarer side stops before W=500), pruned-fork + distal-depth
    guards prevent false 2x1 Kirchhoff chimeras.
    """
    rng = random.Random(101)
    gene = _rand_seq(rng, 1100)
    prom_l = _rand_seq(rng, 250)
    prom_r = _rand_seq(rng, 250)
    L1, R1 = _rand_seq(rng, 4500), _rand_seq(rng, 4500)
    L2, R2 = _rand_seq(rng, 4500), _rand_seq(rng, 4500)
    P1 = L1 + prom_l + gene + prom_r + R1
    P2 = L2 + prom_l + gene + prom_r + R2

    gfa_in = str(tmp_path / "g.fa")
    open(gfa_in, "w").write(">argP\n%s\n" % gene)
    pre = str(tmp_path / "reads")
    _sim_reads(rng, [(P1, 180), (P2, 50)], pre)

    out = str(tmp_path / "out")
    cli(["-g", gfa_in, "-1", pre + "_1.fq", "-2", pre + "_2.fq", "-o", out, "-D", "1500", "-t", "2"])
    doc = json.load(open(out + ".flanks.json"))
    g_ent = doc["genes"]["argP"]
    prs = g_ent.get("pairs", [])
    assert len(prs) == 2, prs
    assert all(pr["tier"] == "depth_match" for pr in prs), prs

    l_by_id = {p["id"]: p["flank"] for p in g_ent["left"]["paths"]}
    r_by_id = {p["id"]: p["flank"] for p in g_ent["right"]["paths"]}
    matched = set()
    for pr in prs:
        lf, rf = l_by_id[pr["left_id"]][-500:], r_by_id[pr["right_id"]][:500]
        if lf == (L1 + prom_l)[-500:] and rf == (prom_r + R1)[:500]:
            matched.add("P1")
        elif lf == (L2 + prom_l)[-500:] and rf == (prom_r + R2)[:500]:
            matched.add("P2")
        else:
            raise AssertionError("Chimeric pair emitted: %r" % (pr,))
    assert matched == {"P1", "P2"}

    # Verify legacy JSON without depth_bins safely abstains instead of guessing from pre-fork summed depth
    legacy = json.loads(json.dumps(g_ent))
    for sd in ("left", "right"):
        for p in legacy[sd]["paths"]:
            p.pop("depth_bins", None)
    assert pair.pair_gene_single_sample("argP", gene, legacy, W=500) == []

    # Partially verified 120x vs 30x case: must emit 0 chimeras even if one rarer side stops at the fork
    pre_low = str(tmp_path / "reads_low")
    _sim_reads(random.Random(101), [(P1, 120), (P2, 30)], pre_low)
    out_low = str(tmp_path / "out_low")
    cli(["-g", gfa_in, "-1", pre_low + "_1.fq", "-2", pre_low + "_2.fq", "-o", out_low, "-D", "1500", "-t", "2"])
    ent_low = json.load(open(out_low + ".flanks.json"))["genes"]["argP"]
    lm_l = {p["id"]: p["flank"] for p in ent_low["left"]["paths"]}
    rm_l = {p["id"]: p["flank"] for p in ent_low["right"]["paths"]}
    for pr in ent_low.get("pairs", []):
        lf, rf = lm_l[pr["left_id"]][-500:], rm_l[pr["right_id"]][:500]
        ok1 = (lf == (L1 + prom_l)[-500:] and rf == (prom_r + R1)[:500])
        ok2 = (lf == (L2 + prom_l)[-500:] and rf == (prom_r + R2)[:500])
        assert ok1 or ok2, "Chimeric pair on partially verified shared promoter: %r" % (pr,)


def test_span_pair_and_equal_depth_abstention(cli, tmp_path):
    rng = random.Random(202)
    # 1) Short 420-bp gene at equal 55x/55x depth -> resolved by direct span_pair
    short_g = _rand_seq(rng, 420)
    L1, R1 = _rand_seq(rng, 4000), _rand_seq(rng, 4000)
    L2, R2 = _rand_seq(rng, 4000), _rand_seq(rng, 4000)
    g1 = str(tmp_path / "short.fa")
    open(g1, "w").write(">shortG\n%s\n" % short_g)
    rpre1 = str(tmp_path / "short_r")
    _sim_reads(rng, [(L1 + short_g + R1, 55), (L2 + short_g + R2, 55)], rpre1, ins_mean=480, ins_sd=45)
    out1 = str(tmp_path / "short_o")
    cli(["-g", g1, "-1", rpre1 + "_1.fq", "-2", rpre1 + "_2.fq", "-o", out1, "-D", "1000", "-t", "2"])
    ent1 = json.load(open(out1 + ".flanks.json"))["genes"]["shortG"]
    prs1 = ent1.get("pairs", [])
    assert len(prs1) == 2 and all(pr["tier"] == "span_pair" for pr in prs1)
    lm1 = {p["id"]: p["flank"] for p in ent1["left"]["paths"]}
    rm1 = {p["id"]: p["flank"] for p in ent1["right"]["paths"]}
    got1 = set()
    for pr in prs1:
        lf, rf = lm1[pr["left_id"]][-500:], rm1[pr["right_id"]][:500]
        if lf == L1[-500:] and rf == R1[:500]:
            got1.add("P1")
        elif lf == L2[-500:] and rf == R2[:500]:
            got1.add("P2")
        else:
            raise AssertionError("Chimeric span_pair!")
    assert got1 == {"P1", "P2"}

    # 2) Long 1,400-bp identical gene at equal 55x/55x depth -> single-sample pairing must abstain
    long_g = _rand_seq(rng, 1400)
    g2 = str(tmp_path / "long.fa")
    open(g2, "w").write(">longG\n%s\n" % long_g)
    rpre2 = str(tmp_path / "long_r")
    _sim_reads(rng, [(L1 + long_g + R1, 55), (L2 + long_g + R2, 55)], rpre2)
    out2 = str(tmp_path / "long_o")
    cli(["-g", g2, "-1", rpre2 + "_1.fq", "-2", rpre2 + "_2.fq", "-o", out2, "-D", "1000", "-t", "2"])
    ent2 = json.load(open(out2 + ".flanks.json"))["genes"]["longG"]
    assert ent2.get("pairs", []) == []


def test_cohort_pairing_resolves_tied_and_shared_arm_plasmids(cli, tmp_path):
    """Three plasmids across S=6 samples carrying an identical 1,400-bp gene:
    P1 = L1-G-R1, P2 = L2-G-R2, P3 = L1-G-R3 (1-to-2 shared left arm + tied depths).
    """
    rng = random.Random(303)
    gene = _rand_seq(rng, 1400)
    L1, L2 = _rand_seq(rng, 4000), _rand_seq(rng, 4000)
    R1, R2, R3 = _rand_seq(rng, 4000), _rand_seq(rng, 4000), _rand_seq(rng, 4000)
    P = {"P1": (L1, R1), "P2": (L2, R2), "P3": (L1, R3)}

    gfa_in = str(tmp_path / "gene.fa")
    open(gfa_in, "w").write(">blaC\n%s\n" % gene)

    # S=6 samples including single-plasmid anchors, tied 50/50 pairs, and 3-plasmid mixture
    design = [
        ("S1", {"P1": 60, "P2": 0,  "P3": 0}),
        ("S2", {"P1": 0,  "P2": 60, "P3": 0}),
        ("S3", {"P1": 50, "P2": 50, "P3": 0}),   # equal-depth tie!
        ("S4", {"P1": 45, "P2": 0,  "P3": 45}),  # shared left arm L1 + equal depth!
        ("S5", {"P1": 0,  "P2": 50, "P3": 50}),  # equal-depth tie!
        ("S6", {"P1": 40, "P2": 60, "P3": 35}),  # all three plasmids!
    ]
    json_args = []
    docs = {}
    for sid, covs in design:
        rpre = str(tmp_path / (sid + "_r"))
        mols = [(P[k][0] + gene + P[k][1], c) for k, c in covs.items() if c > 0]
        _sim_reads(rng, mols, rpre)
        opre = str(tmp_path / sid)
        cli(["-g", gfa_in, "-1", rpre + "_1.fq", "-2", rpre + "_2.fq", "-o", opre, "-D", "1000", "-t", "2"])
        json_args.append("%s=%s.flanks.json" % (sid, opre))
        docs[sid] = json.load(open(opre + ".flanks.json"))

    cout = str(tmp_path / "coh")
    gdir = str(tmp_path / "gfa_out")
    # Invoke via `argjacent cohort` subcommand dispatch
    cli(["cohort", "-g", gfa_in, "-i"] + json_args + ["-o", cout, "--gfa-dir", gdir])

    cdoc = json.load(open(cout + ".cohort.json"))
    truth_total = sum(sum(1 for c in covs.values() if c > 0) for _, covs in design)
    found_true = 0
    chimeras = 0
    for sid, covs in design:
        active = {k for k, c in covs.items() if c > 0}
        g_ent = docs[sid]["genes"]["blaC"]
        lm = {p["id"]: p["flank"] for p in g_ent["left"]["paths"]}
        rm = {p["id"]: p["flank"] for p in g_ent["right"]["paths"]}
        matched = set()
        for pr in cdoc["pairs"].get(sid, {}).get("blaC", []):
            lf, rf = lm[pr["left_id"]][-500:], rm[pr["right_id"]][:500]
            hit = None
            for k in active:
                if lf == P[k][0][-500:] and rf == P[k][1][:500]:
                    hit = k
                    break
            if hit is not None:
                matched.add(hit)
            else:
                chimeras += 1
        found_true += len(matched)
        # Verify updated GFA written in --gfa-dir
        assert os.path.getsize(os.path.join(gdir, sid + ".gfa")) > 0

    assert chimeras == 0
    assert found_true == truth_total, (found_true, truth_total)


def test_cohort_json_only_post_repeat_deconvolution_and_collinear_trap(cli, tmp_path):
    """Post-repeat IS26 exit deconvolution directly from 0.5.0 flanks.json (NO --keep / NO aw.reads.tsv!),
    plus adversarial 1:1 collinear trap refusal.
    """
    rng = random.Random(404)
    gene = _rand_seq(rng, 1000)
    is26 = _rand_seq(rng, 900)
    up_g = _rand_seq(rng, 5000)
    up_d1, up_d2 = _rand_seq(rng, 5000), _rand_seq(rng, 5000)
    ex_true = _rand_seq(rng, 5000)
    ex_d1, ex_d2 = _rand_seq(rng, 5000), _rand_seq(rng, 5000)

    mol_true = up_g + gene + is26 + ex_true
    mol_d1 = up_d1 + is26 + ex_d1
    mol_d2 = up_d2 + is26 + ex_d2

    gfa_in = str(tmp_path / "g.fa")
    open(gfa_in, "w").write(">argIS\n%s\n" % gene)

    # 1) Independent varying abundances across S=7 samples where decoys always outweigh gene
    #    inside IS26 (c1 + c2 > cg), so single-sample argjacent cannot cross IS26 (NO --keep!)
    cov_table = [
        ("R1", 85, 50, 50),
        ("R2", 25, 85, 30),
        ("R3", 25, 30, 85),
        ("R4", 75, 85, 30),
        ("R5", 70, 30, 85),
        ("R6", 35, 75, 75),
        ("R7", 55, 55, 45),
    ]
    j_args = []
    docs = {}
    for sid, cg, c1, c2 in cov_table:
        rpre = str(tmp_path / (sid + "_r"))
        _sim_reads(rng, [(mol_true, cg), (mol_d1, c1), (mol_d2, c2)], rpre)
        opre = str(tmp_path / sid)
        cli(["-g", gfa_in, "-1", rpre + "_1.fq", "-2", rpre + "_2.fq", "-o", opre, "-D", "3200", "-t", "2"])
        assert not os.path.exists(opre + ".tmp")  # confirm JSON-only mode without --keep!
        j_args.append("%s=%s.flanks.json" % (sid, opre))
        docs[sid] = json.load(open(opre + ".flanks.json"))
        # Single-sample argjacent must stop at/before IS26 far end (< 1500 bp)
        assert max(p["context_len"] for p in docs[sid]["genes"]["argIS"]["right"]["paths"]) < 1500

    cout = str(tmp_path / "rep_coh")
    _run_cohort(["-g", gfa_in, "-i"] + j_args + ["-o", cout, "--target-bp", "3000"])
    cdoc = json.load(open(cout + ".cohort.json"))
    certified = 0
    false_exits = 0
    n_proposed = sum(
        1 for sid, _, _, _ in cov_table
        if any(p["flank"][:2500] == (is26 + ex_true)[:2500] for p in docs[sid]["genes"]["argIS"]["right"]["paths"])
    )
    for sid, _, _, _ in cov_table:
        exts = cdoc["repeats"].get(sid, {}).get("argIS", {}).get("right", [])
        for ex in exts:
            assert ex["source"] == "depth_bins"
            assert ex["cohort_certified_bp"] >= 2500
            c_seq = ex["flank"][:2500]
            if c_seq == (is26 + ex_true)[:2500]:
                certified += 1
            else:
                false_exits += 1
    assert false_exits == 0
    assert n_proposed >= 5 and certified == n_proposed, (certified, n_proposed)

    # 2) Adversarial 1:1 collinear trap: decoy D1 always has identical coverage to gene
    trap_table = [
        ("T1", 25, 25, 60),
        ("T2", 50, 50, 30),
        ("T3", 75, 75, 80),
        ("T4", 35, 35, 90),
        ("T5", 65, 65, 20),
        ("T6", 40, 40, 70),
    ]
    t_args = []
    for sid, cg, c1, c2 in trap_table:
        rpre = str(tmp_path / (sid + "_r"))
        _sim_reads(rng, [(mol_true, cg), (mol_d1, c1), (mol_d2, c2)], rpre)
        opre = str(tmp_path / sid)
        cli(["-g", gfa_in, "-1", rpre + "_1.fq", "-2", rpre + "_2.fq", "-o", opre, "-D", "3200", "-t", "2"])
        t_args.append("%s=%s.flanks.json" % (sid, opre))

    tout = str(tmp_path / "trap_coh")
    _run_cohort(["-g", gfa_in, "-i"] + t_args + ["-o", tout, "--target-bp", "3000"])
    tdoc = json.load(open(tout + ".cohort.json"))
    total_trap_claims = sum(
        len(exts)
        for smap in tdoc["repeats"].values()
        for gmap in smap.values()
        for exts in gmap.values()
    )
    assert total_trap_claims == 0, tdoc["repeats"]
