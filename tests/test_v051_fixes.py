"""Regression unit tests for v0.5.1 pairing, cohort, GFA, and CLI fixes."""
import json
import os
import random
import subprocess
import sys

from argjacent import __version__ as VERSION
from argjacent import cli, cohort, gfa, pair


COMP = str.maketrans("ACGT", "TGCA")


def _rc(s):
    return s.translate(COMP)[::-1]


def _rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _make_path(cid, flank, depth, cert_bp=800, bins=None, gene_depth=None):
    nb = max(1, len(flank) // 100)
    dbins = list(bins) if bins is not None else [float(depth)] * nb
    return {
        "id": cid,
        "flank": flank,
        "flank_len": len(flank),
        "context_len": cert_bp,
        "end_len": cert_bp,
        "shared_end": 0,
        "depth": float(depth),
        "gene_depth": float(gene_depth if gene_depth is not None else depth),
        "depth_bins": dbins,
        "main": cid == "c1",
    }


def test_span_pair_disabled_above_650bp_deletion_variant(tmp_path):
    """Reproduces the Strain Madness s30 aadA2 (792 bp) cassette-deletion failure mode:
    a >650 bp multi-trunk gene must never emit span_pair from reads spanning a co-existing
    deletion variant where L1 and R1 sit directly adjacent without the gene, whereas a
    short (<= 650 bp) gene with real spanning pairs must still emit span_pair.
    """
    rng = random.Random(511)
    L1, L2 = _rand_seq(rng, 800), _rand_seq(rng, 800)
    R1, R2 = _rand_seq(rng, 800), _rand_seq(rng, 800)

    # 1) 792-bp gene (> 650 bp): true contexts are (L1, R2) and (L2, R1), PLUS a deletion
    #    variant L1 + R1 that generates 20 read pairs bridging L1 directly to R1.
    long_g = _rand_seq(rng, 792)
    doc_long = {
        "aadA2": {
            "left": {"read_depth": 100.0, "paths": [_make_path("c1", L1, 50.0), _make_path("c2", L2, 50.0)]},
            "right": {"read_depth": 100.0, "paths": [_make_path("c1", R1, 50.0), _make_path("c2", R2, 50.0)]},
        }
    }
    reads_long = str(tmp_path / "long_reads.tsv")
    del_mol = L1 + R1
    with open(reads_long, "w") as fh:
        for pid in range(25):
            # Mate 1 lies inside L1, mate 2 lies inside R1 (simulating deletion cassette)
            s1 = del_mol[len(L1) - 180 + (pid % 20):len(L1) - 30 + (pid % 20)]
            s2 = _rc(del_mol[len(L1) + 20 + (pid % 20):len(L1) + 170 + (pid % 20)])
            fh.write("%d\t1\t%s\t\t%s\n" % (pid, s1, "I" * len(s1)))
            fh.write("%d\t2\t%s\t\t%s\n" % (pid, s2, "I" * len(s2)))

    phased_long = pair.index_sample_reads({"aadA2": long_g}, doc_long, reads_long, W=500)
    assert phased_long == {}, phased_long
    pairs_long = pair.pair_sample({"aadA2": long_g}, doc_long, W=500, reads_source=reads_long)
    assert pairs_long == {}, pairs_long

    # 2) 450-bp gene (<= 650 bp): real spanning pairs across (L1, R1) and (L2, R2) still fire.
    short_g = _rand_seq(rng, 450)
    doc_short = {
        "fosB": {
            "left": {"read_depth": 100.0, "paths": [_make_path("c1", L1, 50.0), _make_path("c2", L2, 50.0)]},
            "right": {"read_depth": 100.0, "paths": [_make_path("c1", R1, 50.0), _make_path("c2", R2, 50.0)]},
        }
    }
    reads_short = str(tmp_path / "short_reads.tsv")
    m1 = L1 + short_g + R1
    m2 = L2 + short_g + R2
    with open(reads_short, "w") as fh:
        pid = 0
        for mol in (m1, m2):
            for k in range(12):
                s1 = mol[len(L1) - 160 + k:len(L1) - 10 + k]
                s2 = _rc(mol[len(L1) + len(short_g) + 10 + k:len(L1) + len(short_g) + 160 + k])
                fh.write("%d\t1\t%s\t\t%s\n" % (pid, s1, "I" * len(s1)))
                fh.write("%d\t2\t%s\t\t%s\n" % (pid, s2, "I" * len(s2)))
                pid += 1

    pairs_short = pair.pair_sample({"fosB": short_g}, doc_short, W=500, reads_source=reads_short)
    assert len(pairs_short.get("fosB", [])) == 2
    assert all(p["tier"] == "span_pair" for p in pairs_short["fosB"])
    got = {(p["left_id"], p["right_id"]) for p in pairs_short["fosB"]}
    assert got == {("c1", "c1"), ("c2", "c2")}


def test_cohort_preserves_stored_span_pair_and_snv_phase_and_anchors(tmp_path):
    """argjacent-cohort must keep single-sample span_pair and snv_phase records stored in
    flanks.json, propagate them as cohort_anchor seeds, emit them into --gfa-dir GFAs,
    refuse to extrapolate saved pairs from smaller pair_bp < W, and drop stale >650 bp span_pairs.
    """
    rng = random.Random(512)
    short_g = _rand_seq(rng, 450)
    long_g = _rand_seq(rng, 1200)
    L1, L2 = _rand_seq(rng, 800), _rand_seq(rng, 800)
    R1, R2 = _rand_seq(rng, 800), _rand_seq(rng, 800)

    def make_sample(sid, gname, gseq, stored_pairs, pair_bp=500, d1=50.0, d2=50.0):
        return {
            "sid": sid,
            "genes": {gname: gseq},
            "doc": {
                "version": VERSION,
                "parameters": {"pair_bp": pair_bp, "read_len": 150},
                "genes": {
                    gname: {
                        "left": {"read_depth": d1 + d2, "paths": [_make_path("c1", L1, d1), _make_path("c2", L2, d2)]},
                        "right": {"read_depth": d1 + d2, "paths": [_make_path("c1", R1, d1), _make_path("c2", R2, d2)]},
                        "pairs": stored_pairs,
                    }
                },
            },
            "reads_tsv": None,
        }

    stored_span = [
        {"id": "p1", "left_trunk": "L1", "right_trunk": "R1", "left_id": "c1", "right_id": "c1",
         "left_certain_bp": 800, "right_certain_bp": 800, "left_depth": 60.0, "right_depth": 60.0,
         "tier": "span_pair", "score": 99.0, "note": "pairs=15"},
        {"id": "p2", "left_trunk": "L2", "right_trunk": "R2", "left_id": "c2", "right_id": "c2",
         "left_certain_bp": 800, "right_certain_bp": 800, "left_depth": 40.0, "right_depth": 40.0,
         "tier": "span_pair", "score": 99.0, "note": "pairs=14"},
    ]
    # S1 has (60x, 40x) — too close for 3-sigma depth_match alone, so it relies on stored span_pair,
    # while differing in ratio from S2's tied (50x, 50x) so both links transfer via cohort_anchor.
    s1 = make_sample("S1", "fosB_1_X89875", short_g, stored_span, pair_bp=500, d1=60.0, d2=40.0)
    s2 = make_sample("S2", "fosB_1_X89875", short_g, [], pair_bp=500, d1=50.0, d2=50.0)

    res1 = cohort.cohort_pair_all([s1, s2], W=500, threads=1)
    res2 = cohort.cohort_pair_all([s1, s2], W=500, threads=2)
    assert res1 == res2

    s1_pairs = res1["S1"]["fosB_1_X89875"]
    s2_pairs = res1["S2"]["fosB_1_X89875"]
    assert {(p["left_id"], p["right_id"], p["tier"]) for p in s1_pairs} == {
        ("c1", "c1", "span_pair"), ("c2", "c2", "span_pair")
    }
    assert {(p["left_id"], p["right_id"], p["tier"]) for p in s2_pairs} == {
        ("c1", "c1", "cohort_anchor"), ("c2", "c2", "cohort_anchor")
    }

    # Even in a rank-1 identical (50x, 50x) vs (50x, 50x) tie, S2 still resolves both pairs
    # (one via cohort_anchor and the other via 2x2 cohort_complement).
    s1_eq = make_sample("E1", "fosB_1_X89875", short_g, stored_span, pair_bp=500, d1=50.0, d2=50.0)
    s2_eq = make_sample("E2", "fosB_1_X89875", short_g, [], pair_bp=500, d1=50.0, d2=50.0)
    res_eq = cohort.cohort_pair_all([s1_eq, s2_eq], W=500, threads=1)
    assert {(p["left_id"], p["right_id"]) for p in res_eq["E2"]["fosB_1_X89875"]} == {("c1", "c1"), ("c2", "c2")}

    # Verify CLI + --gfa-dir preserves span_pair inside emitted GFA
    gfa_fa = str(tmp_path / "genes.fa")
    open(gfa_fa, "w").write(">fosB_1_X89875\n%s\n" % short_g)
    j1, j2 = str(tmp_path / "S1.flanks.json"), str(tmp_path / "S2.flanks.json")
    json.dump(s1["doc"], open(j1, "w"))
    json.dump(s2["doc"], open(j2, "w"))
    out_pre = str(tmp_path / "coh")
    gdir = str(tmp_path / "gfas")
    r = subprocess.run(
        [sys.executable, "-m", "argjacent.cohort", "-g", gfa_fa, "-i", j1, j2, "-o", out_pre, "--gfa-dir", gdir, "-t", "2"],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    gfa_s1 = open(os.path.join(gdir, "S1.gfa")).read()
    assert "|pair|p1|span_pair" in gfa_s1
    assert "|pair|p2|span_pair" in gfa_s1

    # snv_phase on a 1,200-bp gene is preserved and seeds cohort_anchor
    stored_snv = [
        {"id": "p1", "left_trunk": "L1", "right_trunk": "R1", "left_id": "c1", "right_id": "c1",
         "left_certain_bp": 800, "right_certain_bp": 800, "left_depth": 50.0, "right_depth": 50.0,
         "tier": "snv_phase", "score": 99.0, "note": "snvs=4,r_links=9"},
    ]
    snv1 = make_sample("A1", "tetM_1_X90939", long_g, stored_snv, pair_bp=500)
    snv2 = make_sample("A2", "tetM_1_X90939", long_g, [], pair_bp=500)
    res_snv = cohort.cohort_pair_all([snv1, snv2], W=500, threads=1)
    assert res_snv["A1"]["tetM_1_X90939"][0]["tier"] == "snv_phase"
    assert {(p["left_id"], p["right_id"]) for p in res_snv["A2"]["tetM_1_X90939"]} == {("c1", "c1"), ("c2", "c2")}

    # Radius guard: saved pair_bp=500 must NOT be extrapolated when W=700
    res_wide = cohort.cohort_pair_all([s1, s2], W=700, threads=1)
    assert res_wide["S1"] == {} and res_wide["S2"] == {}

    # Legacy stale >650 bp span_pair must NOT be carried over
    stale1 = make_sample("B1", "aadA2_1_JQ364967", long_g, stored_span, pair_bp=500)
    stale2 = make_sample("B2", "aadA2_1_JQ364967", long_g, [], pair_bp=500)
    res_stale = cohort.cohort_pair_all([stale1, stale2], W=500, threads=1)
    assert res_stale["B1"] == {} and res_stale["B2"] == {}


def test_gene_family_modes_safe_for_refseq_card_and_resfinder():
    # ResFinder allele copies collapse cleanly in auto mode
    assert cohort.gene_family("blaTEM-1B_1_JF910132", "auto") == "blaTEM-1B"
    assert cohort.gene_family("blaTEM-1B_2_AF123456", "auto") == "blaTEM-1B"
    assert cohort.gene_family("blaZ_54_NZ_JVRB01000003", "auto") == "blaZ"
    assert cohort.gene_family("tet(W)_1_NC_012345", "auto") == "tet(W)"
    assert cohort.gene_family("aac(6')-Ib-cr_1_DQ303918", "auto") == "aac(6')-Ib-cr"

    # Non-ResFinder headers stay distinct in auto mode (no silent prefix collisions!)
    for name in (
        "NC_000913_blaTEM", "NC_000913_tetA", "NZ_CP012345.1_1", "WP_000123456.1",
        "ARO_3000194", "CARD_tetM", "PROKKA_00001", "contig_1_orf1", "contig_1_orf2",
        "gene_1", "gene_2", "my_target_A",
    ):
        assert cohort.gene_family(name, "auto") == name
        assert cohort.gene_family(name, "exact") == name

    # Legacy resfinder mode splits on first underscore; exact never merges
    assert cohort.gene_family("blaTEM-1B_1_JF910132", "exact") == "blaTEM-1B_1_JF910132"
    assert cohort.gene_family("NC_000913_blaTEM", "resfinder") == "NC"

    # End-to-end: two unrelated RefSeq-prefixed genes in the same samples must NOT collide under auto
    rng = random.Random(513)
    gA, gB = _rand_seq(rng, 1100), _rand_seq(rng, 1100)
    LA1, LA2, RA1, RA2 = (_rand_seq(rng, 700) for _ in range(4))
    LB1, LB2, RB1, RB2 = (_rand_seq(rng, 700) for _ in range(4))

    def mk_sd(sid, a_depths, b_depths):
        def mk_gene(L1s, L2s, R1s, R2s, d1, d2):
            lp, rp = [], []
            if d1 > 0:
                lp.append(_make_path("c1", L1s, d1))
                rp.append(_make_path("c1", R1s, d1))
            if d2 > 0:
                lp.append(_make_path("c%d" % (len(lp) + 1), L2s, d2))
                rp.append(_make_path("c%d" % (len(rp) + 1), R2s, d2))
            return {"left": {"read_depth": d1 + d2, "paths": lp}, "right": {"read_depth": d1 + d2, "paths": rp}}
        return {
            "sid": sid,
            "genes": {"NC_000913_geneA": gA, "NC_000913_geneB": gB},
            "doc": {
                "genes": {
                    "NC_000913_geneA": mk_gene(LA1, LA2, RA1, RA2, *a_depths),
                    "NC_000913_geneB": mk_gene(LB1, LB2, RB1, RB2, *b_depths),
                }
            },
            "reads_tsv": None,
        }

    # In every sample geneA is deeper than geneB; under legacy split('_')[0], geneA would shadow geneB!
    samples = [
        mk_sd("S1", (80.0, 0.0), (40.0, 0.0)),
        mk_sd("S2", (0.0, 80.0), (0.0, 40.0)),
        mk_sd("S3", (60.0, 60.0), (30.0, 30.0)),
    ]
    res_auto = cohort.cohort_pair_all(samples, W=500, threads=2, family_mode="auto")
    # S3 is an equal-depth tie for BOTH geneA and geneB -> both must resolve via cohort_anchor!
    assert len(res_auto["S3"].get("NC_000913_geneA", [])) == 2
    assert len(res_auto["S3"].get("NC_000913_geneB", [])) == 2


def test_part_b_modal_fork_bp_resists_single_outlier_sample():
    # Direct unit checks on _select_modal_fork_bp
    sets = [{820} for _ in range(20)] + [{410}]
    assert cohort._select_modal_fork_bp(sets, tol_bp=50) == 820
    assert cohort._select_modal_fork_bp([{500}, {500}, {900}, {900}], tol_bp=50) == 500
    assert cohort._select_modal_fork_bp([], tol_bp=50) is None

    # End-to-end Part B test: 7 samples fork at 850 bp (IS26), plus 1 outlier sample forking at 420 bp
    rng = random.Random(514)
    gene = _rand_seq(rng, 1000)
    is26 = _rand_seq(rng, 850)
    ex_true = _rand_seq(rng, 2200)
    ex_d1 = _rand_seq(rng, 2200)
    ex_d2 = _rand_seq(rng, 2200)
    weird_early = _rand_seq(rng, 2630)

    covs = [
        ("R1", 80.0, 45.0, 45.0),
        ("R2", 25.0, 85.0, 30.0),
        ("R3", 25.0, 30.0, 85.0),
        ("R4", 75.0, 85.0, 30.0),
        ("R5", 70.0, 30.0, 85.0),
        ("R6", 35.0, 75.0, 75.0),
        ("R7", 55.0, 55.0, 45.0),
    ]
    samples = []
    for sid, cg, c1, c2 in covs:
        # Bins 0..8 lie in IS26 (summed depth), bins 9..29 lie in exclusive post-IS26 exits
        tot = cg + c1 + c2
        p_t = _make_path("c1", is26 + ex_true, cg, cert_bp=400, bins=[tot] * 9 + [cg] * 21, gene_depth=cg)
        p_1 = _make_path("c2", is26 + ex_d1, c1, cert_bp=400, bins=[tot] * 9 + [c1] * 21, gene_depth=cg)
        p_2 = _make_path("c3", is26 + ex_d2, c2, cert_bp=400, bins=[tot] * 9 + [c2] * 21, gene_depth=cg)
        samples.append({
            "sid": sid,
            "genes": {"argIS": gene},
            "doc": {"genes": {"argIS": {"right": {"read_depth": cg, "paths": [p_t, p_1, p_2]}}}},
            "reads_tsv": None,
        })

    # Outlier sample with a 420-bp fork (would have broken min(fork_votes)!)
    out_p1 = _make_path("c1", is26[:420] + weird_early, 40.0, cert_bp=300, bins=[80.0] * 4 + [40.0] * 26, gene_depth=40.0)
    out_p2 = _make_path("c2", is26 + ex_d1, 40.0, cert_bp=300, bins=[80.0] * 4 + [40.0] * 26, gene_depth=40.0)
    samples.append({
        "sid": "S_OUTLIER",
        "genes": {"argIS": gene},
        "doc": {"genes": {"argIS": {"right": {"read_depth": 40.0, "paths": [out_p1, out_p2]}}}},
        "reads_tsv": None,
    })

    resolved = cohort.cohort_deconvolve_repeats(samples, "argIS", "right", min_repeat_bp=400, target_bp=2800)
    assert len(resolved) == 7, resolved
    for sid, _, _, _ in covs:
        exts = resolved[sid]
        assert len(exts) == 1
        assert 840 <= exts[0]["fork_bp"] <= 860
        assert exts[0]["flank"] == is26 + ex_true
        assert exts[0]["cohort_certified_bp"] >= 2500


def test_gfa_sanitized_id_collision_and_read_len_rc(tmp_path):
    rng = random.Random(515)
    g1, g2 = _rand_seq(rng, 600), _rand_seq(rng, 600)
    f1, f2 = _rand_seq(rng, 500), _rand_seq(rng, 500)
    genes = {"bla-1": g1, "bla_1": g2}
    doc_genes = {
        "bla-1": {"left": {"read_depth": 30.0, "paths": []},
                  "right": {"read_depth": 30.0, "paths": [_make_path("c1", f1, 30.0, cert_bp=500)]}},
        "bla_1": {"left": {"read_depth": 40.0, "paths": []},
                  "right": {"read_depth": 40.0, "paths": [_make_path("c1", f2, 40.0, cert_bp=500)]}},
    }

    gfa150 = str(tmp_path / "r150.gfa")
    gfa100 = str(tmp_path / "r100.gfa")
    gfa.write_gfa(genes, doc_genes, gfa150, verify=True, read_len=150)
    gfa.write_gfa(genes, doc_genes, gfa100, verify=True, read_len=100)

    def read_s_records(path):
        s_ids = []
        rc_by_id = {}
        for ln in open(path):
            if ln.startswith("S\t"):
                toks = ln.rstrip("\n").split("\t")
                s_ids.append(toks[1])
                for t in toks[3:]:
                    if t.startswith("RC:i:"):
                        rc_by_id[toks[1]] = int(t.split(":")[2])
        return s_ids, rc_by_id

    ids150, rc150 = read_s_records(gfa150)
    ids100, rc100 = read_s_records(gfa100)
    # Zero duplicate segment IDs even though 'bla-1' and 'bla_1' both sanitize to 'bla_1'
    assert len(ids150) == len(set(ids150)) == 4, ids150
    assert "bla_1|GENE" in ids150 and "bla_1__2|GENE" in ids150
    # RC:i: scales inversely with read_len (100 bp reads need 1.5x as many reads as 150 bp)
    assert rc100["bla_1|GENE"] == int(round(30.0 * 600 / 100.0))
    assert rc150["bla_1|GENE"] == int(round(30.0 * 600 / 150.0))


def test_load_cohort_inputs_caches_genes_fa_and_binary_versions(tmp_path):
    gfa_fa = str(tmp_path / "shared.fa")
    open(gfa_fa, "w").write(">g1\n" + "ACGT" * 200 + "\n")
    sheet = str(tmp_path / "sheet.tsv")
    with open(sheet, "w") as fh:
        fh.write("sample\tflanks_json\tgenes_fa\n")
        for i in range(4):
            jp = str(tmp_path / ("s%d.flanks.json" % i))
            json.dump({"version": VERSION, "genes": {}}, open(jp, "w"))
            fh.write("s%d\t%s\t%s\n" % (i, jp, gfa_fa))

    class _Args:
        genes = None
        inputs = None
        sample_sheet = sheet

    loaded = cohort._load_cohort_inputs(_Args())
    assert len(loaded) == 4
    # All four sample records must share the single cached FASTA dict instance
    assert all(sd["genes"] is loaded[0]["genes"] for sd in loaded[1:])

    # All three C binaries must match argjacent.__version__
    for name in ("aw", "awx", "awv"):
        b = cli.find_bin(name)
        assert b is not None
        r = subprocess.run([b, "-v"], capture_output=True, text=True)
        assert r.returncode == 0
        assert VERSION in r.stdout, (name, r.stdout, VERSION)


def test_terminal_crowded_allele_hijack_rejected_while_uniform_ortholog_kept(tmp_path):
    """Reproduces the real-world CPH_RA sewage cfxA6 vs cfxA5 C-terminal hijack:
    1) If r300 walks onto a co-occurring relative with 20 mismatches crowded into the last
       30 bp (total d300=20 <= 30, but d100=20 > 15 and d30=20 > 9), pipeline.merge() must
       reject that relative-terminus hijack and fall back to r0 (the database allele's true flank).
    2) An r300 candidate with 17 mismatches crowded into the last 100 bp (83 % local identity,
       d100=17 > 15) is likewise rejected and falls back to r0.
    3) A genuinely ~90.7 %-identical ortholog with 28 mismatches distributed uniformly across
       the 300-bp REDO window (9 in last 100 bp, 3 in last 30 bp, e.g. Klebsiella OqxA) is kept.
    """
    from argjacent import pipeline

    rng = random.Random(516)
    gene_cfx = _rand_seq(rng, 996)
    gene_ecoli = _rand_seq(rng, 900)
    gene_oqxa = _rand_seq(rng, 1100)
    genes = [("cfxA6", gene_cfx), ("rel83", gene_ecoli), ("OqxA", gene_oqxa)]

    def mut_at(seq, positions):
        arr = list(seq)
        for p in positions:
            arr[p] = rng.choice([b for b in "ACGT" if b != arr[p]])
        return "".join(arr)

    fl_true_cfx = _rand_seq(rng, 718)
    fl_hijack_cfx5 = _rand_seq(rng, 927)
    fl_true_ecoli = _rand_seq(rng, 600)
    fl_hijack_ecoli = _rand_seq(rng, 1800)
    fl_oqxa_sample = _rand_seq(rng, 1500)

    # 1) cfxA6: 20 mismatches crowded into the terminal 30 bp [370..400)
    end_cfx = gene_cfx[-400:]
    hijack_cfx = mut_at(end_cfx, list(range(370, 390)))

    # 2) rel83: 17 mismatches in [300..370) (d300=17 <= 30, d30=0 <= 9, but d100=17 > 15)
    end_ecoli = gene_ecoli[-400:]
    hijack_ecoli = mut_at(end_ecoli, list(range(300, 368, 4)))  # 17 sites in [300..364]

    # 3) OqxA: 28 mismatches uniformly spaced across [100..400) -> d300=28, d100=9, d30=3
    end_oqxa = gene_oqxa[-400:]
    unif_oqxa = mut_at(end_oqxa, list(range(105, 396, 10))[:28])

    pre = str(tmp_path / "x")
    cdir = str(tmp_path / "c")
    os.makedirs(cdir)

    def write_run(redo, rows):
        with open("%s.r%d.fa" % (pre, redo), "w") as fa, open("%s.r%d.tsv" % (pre, redo), "w") as tsv:
            tsv.write("gene\tside\tcand\tflank\tstop\tborn_at\tpenalty\tvariants\tforks\tone_read\tjoin_first\tmass_end\n")
            for gi, side, cid, gl, seq, stop in rows:
                fa.write(">%d|%s|%d|%d|0|0\n%s\n" % (gi, side, cid, gl, seq))
                tsv.write("%d\t%s\t%d\t%d\t%s\t0\t0\t0\t0\t0\t-99999\t1.0000\n" % (gi, side, cid, len(seq) - gl, stop))

    write_run(0, [
        (0, "right", 1, 400, end_cfx + fl_true_cfx, "reached target"),
        (1, "right", 1, 400, end_ecoli + fl_true_ecoli, "reached target"),
        (2, "right", 1, 400, end_oqxa, "ran out of reads"),
    ])
    write_run(300, [
        (0, "right", 1, 400, hijack_cfx + fl_hijack_cfx5, "reached target"),
        (1, "right", 1, 400, hijack_ecoli + fl_hijack_ecoli, "reached target"),
        (2, "right", 1, 400, unif_oqxa + fl_oqxa_sample, "reached target"),
    ])

    cands, st = pipeline.merge(genes, pre, cdir)
    # cfxA6 and rel83 both reject the crowded relative end and recover their true r0 flank!
    assert len(cands[0]) == 1 and cands[0][0]["gene_end"] == "database"
    assert cands[0][0]["seq"][400:] == fl_true_cfx
    assert len(cands[1]) == 1 and cands[1][0]["gene_end"] == "database"
    assert cands[1][0]["seq"][400:] == fl_true_ecoli
    # OqxA keeps its uniform ~90.7 %-identical sample allele!
    assert len(cands[2]) == 1 and cands[2][0]["gene_end"] == "sample"
    assert cands[2][0]["gene_end_diffs"] == 28
    assert cands[2][0]["gene_end_diffs_near"] <= pipeline.MAXNEAR
    assert cands[2][0]["gene_end_diffs_tip"] <= pipeline.MAXTIP
    assert cands[2][0]["seq"][400:] == fl_oqxa_sample
    assert st["sides_database_end"] == 2 and st["sides_sample_end"] == 1
    assert st["sides_sample_end_too_far"] == 2 and st["dropped_too_far"] == 0


def test_end_to_end_cooccurring_cterminal_alleles_do_not_bleed(tmp_path):
    """End-to-end simulation of cfxA6 (12x) co-occurring with cfxA5 (45x), sharing 966 bp identically
    except for a 20-bp mismatch in the terminal 30 bp of cfxA6: ARGjacent must certify cfxA6's own
    true right flank rather than hijacking onto cfxA5's higher-depth right flank.
    """
    rng = random.Random(517)
    shared_core = _rand_seq(rng, 966)
    tail_cfx6 = _rand_seq(rng, 30)
    # Mutate 20 of the 30 terminal bases to form cfxA5's 3' terminus
    tail_cfx5 = list(tail_cfx6)
    for idx in range(0, 30, 3):
        tail_cfx5[idx] = rng.choice([b for b in "ACGT" if b != tail_cfx5[idx]])
    for idx in range(1, 30, 3):
        tail_cfx5[idx] = rng.choice([b for b in "ACGT" if b != tail_cfx5[idx]])
    tail_cfx5 = "".join(tail_cfx5)  # exactly 20 diffs in last 30 bp

    gene_cfx6 = shared_core + tail_cfx6
    gene_cfx5 = shared_core + tail_cfx5
    left_6, right_6 = _rand_seq(rng, 1500), _rand_seq(rng, 1500)
    left_5, right_5 = _rand_seq(rng, 1500), _rand_seq(rng, 1500)

    mol6 = left_6 + gene_cfx6 + right_6
    mol5 = left_5 + gene_cfx5 + right_5

    genes_fa = str(tmp_path / "genes.fa")
    with open(genes_fa, "w") as fh:
        fh.write(">cfxA6\n%s\n>cfxA5\n%s\n" % (gene_cfx6, gene_cfx5))

    r1_fq = str(tmp_path / "R1.fq")
    r2_fq = str(tmp_path / "R2.fq")
    pid = 0
    with open(r1_fq, "w") as f1, open(r2_fq, "w") as f2:
        for mol, cov in ((mol6, 18), (mol5, 60)):
            n_pairs = int(cov * len(mol) / 300)
            for _ in range(n_pairs):
                ins = max(200, min(450, int(rng.gauss(320, 35))))
                pos = rng.randrange(0, len(mol) - ins)
                s1 = mol[pos:pos + 150]
                s2 = _rc(mol[pos + ins - 150:pos + ins])
                f1.write("@r%d/1\n%s\n+\n%s\n" % (pid, s1, "I" * 150))
                f2.write("@r%d/2\n%s\n+\n%s\n" % (pid, s2, "I" * 150))
                pid += 1

    # 1) Default --max-paths 16: r300 drops the 60x cfxA5 hijack (diffs=20) and keeps cfxA6's true 0-diff path!
    out_pre = str(tmp_path / "out16")
    res = subprocess.run(
        [sys.executable, "-m", "argjacent", "-g", genes_fa, "-1", r1_fq, "-2", r2_fq, "-o", out_pre, "-D", "1200", "-t", "2"],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, res.stdout + res.stderr
    doc = json.load(open(out_pre + ".flanks.json"))["genes"]

    p6 = doc["cfxA6"]["right"]["paths"]
    assert len(p6) == 1, p6
    best6 = p6[0]
    assert best6["gene_end_diffs"] == 0 and best6["gene_end_diffs_near"] == 0, best6
    assert best6["context_len"] >= 900, best6
    assert best6["flank"][:best6["context_len"]] == right_6[:best6["context_len"]]
    assert not best6["flank"].startswith(right_5[:100])

    best5 = doc["cfxA5"]["right"]["paths"][0]
    assert best5["context_len"] >= 900, best5
    assert best5["flank"][:best5["context_len"]] == right_5[:best5["context_len"]]

    # 2) Crowded beam (--max-paths 1): even when cfxA5 (60x) crowds cfxA6 (18x) out of awx's
    #    single beam slot, pipeline.merge() rejects the hijacked r300 candidate and verify.py
    #    refuses to certify any chimeric right_5 extension on cfxA6 (context_len == 0)!
    out_crowded = str(tmp_path / "out1")
    res1 = subprocess.run(
        [sys.executable, "-m", "argjacent", "-g", genes_fa, "-1", r1_fq, "-2", r2_fq, "-o", out_crowded,
         "-D", "1200", "--max-paths", "1", "-t", "2"],
        capture_output=True, text=True,
    )
    assert res1.returncode == 0, res1.stdout + res1.stderr
    doc1 = json.load(open(out_crowded + ".flanks.json"))["genes"]
    b6_1 = doc1["cfxA6"]["right"]["paths"][0]
    assert b6_1["gene_end"] == "database" and b6_1["gene_end_diffs"] == 0, b6_1
    assert b6_1["context_len"] == 0, b6_1
