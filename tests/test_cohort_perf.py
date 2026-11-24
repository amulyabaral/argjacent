"""Unit and performance tests for exact optimizations in argjacent.pair and argjacent.cohort."""
import random
import time

import numpy as np

from argjacent import cohort, pair


def _orig_seq_sim(a, b):
    L = min(len(a), len(b))
    if L == 0:
        return 1.0
    direct = sum(x == y for x, y in zip(a[:L], b[:L])) / L
    if direct >= 0.97 or L < 40:
        return direct
    best = direct
    for sh in (-12, -8, -4, -2, -1, 1, 2, 4, 8, 12):
        if sh > 0 and L - sh > 20:
            m = sum(x == y for x, y in zip(a[:L - sh], b[sh:L])) / (L - sh)
            if m > best:
                best = m
        elif sh < 0 and L + sh > 20:
            s = -sh
            m = sum(x == y for x, y in zip(a[s:L], b[:L - s])) / (L - s)
            if m > best:
                best = m
    return best


def _orig_div_bp(out_a, out_b):
    L = min(len(out_a), len(out_b))
    diffs = [1 if out_a[i] != out_b[i] else 0 for i in range(L)]
    win = sum(diffs[:15])
    if win >= 4:
        for i in range(min(15, L)):
            if diffs[i]:
                return i
    for i in range(15, L):
        win += diffs[i] - diffs[i - 15]
        if win >= 4:
            for j in range(i - 14, i + 1):
                if diffs[j]:
                    return j
    return L


def _orig_nnls(A, b, w_sqrt):
    Aw = A * w_sqrt[:, None]
    bw = b * w_sqrt
    n = Aw.shape[1]
    if n == 0:
        return np.zeros(0), float(np.sum(bw ** 2))
    P_set = np.zeros(n, dtype=bool)
    x = np.zeros(n, dtype=np.float64)
    w_grad = Aw.T @ (bw - Aw @ x)
    it = 0
    while np.any((~P_set) & (w_grad > 1e-9)) and it < 5 * n + 20:
        it += 1
        j = int(np.argmax(np.where(~P_set, w_grad, -1e30)))
        P_set[j] = True
        while True:
            idx = np.nonzero(P_set)[0]
            if len(idx) == 0:
                break
            Ap = Aw[:, idx]
            sp, _, _, _ = np.linalg.lstsq(Ap, bw, rcond=1e-10)
            if np.all(sp > 1e-10):
                x[:] = 0.0
                x[idx] = sp
                break
            neg = sp <= 1e-10
            alphas = x[idx][neg] / (x[idx][neg] - sp[neg] + 1e-30)
            alpha = float(np.min(alphas))
            x_new = np.zeros(n)
            x_new[idx] = x[idx] + alpha * (sp - x[idx])
            x = np.maximum(0.0, x_new)
            P_set[idx[x[idx] <= 1e-10]] = False
        w_grad = Aw.T @ (bw - Aw @ x)
    res = float(np.sum((bw - Aw @ x) ** 2))
    return x, res


def test_seq_sim_and_divergence_bp_exact():
    rng = random.Random(17)
    base = "".join(rng.choice("ACGT") for _ in range(500))
    cases = [
        ("", ""),
        ("", "ACGT"),
        ("ACGTACGT", "ACGTACGT"),
        ("ACGTACGT", "TTTTAAAA"),
        (base[:10], base[:10]),
        (base[:35], "T" * 35),
        (base, base),
        (base, base[:250] + "A" * 250),
        (base, "AA" + base[:-2]),       # +2 shift
        (base, base[4:] + "GGGG"),      # -4 shift
        (base, base[:300]),             # unequal lengths
        (base[:300], base),             # unequal lengths reversed
        (base[:300], base[8:308]),      # shift + unequal lengths
        (base, base[:120] + "G" + base[121:]),  # isolated SNP
    ]
    for a, b in cases:
        assert abs(pair.seq_sim_fast(a, b) - _orig_seq_sim(a, b)) < 1e-12
        assert pair.divergence_bp(a, b) == _orig_div_bp(a, b)


def test_exact_seq_clusters_matches_naive():
    rng = random.Random(29)
    roots = ["".join(rng.choice("ACGT") for _ in range(rng.choice([120, 300, 500]))) for _ in range(25)]
    seqs = []
    for r in roots:
        seqs.append(r)
        # Shifted variants (+/- 1, 2, 4, 8 bp) and truncated variants
        if len(r) >= 100:
            seqs.append(r[2:] + "AT")
            seqs.append("CG" + r[:-2])
            seqs.append(r[:80])
            seqs.append(r[4:160])
    rng.shuffle(seqs)

    naive_clusters = []
    naive_ids = []
    for s in seqs:
        hit = None
        for idx, c in enumerate(naive_clusters):
            if _orig_seq_sim(s, c) >= 0.96:
                hit = idx
                break
        if hit is None:
            hit = len(naive_clusters)
            naive_clusters.append(s)
        naive_ids.append(hit)

    fast_cl = cohort._ExactSeqClusters(max_len=500, thr=0.96)
    fast_ids = [fast_cl.match_or_add(s, s) for s in seqs]
    assert fast_ids == naive_ids
    assert fast_cl.cnt == len(naive_clusters)
    # Lookup-only mode must also match
    for s, exp in zip(seqs, naive_ids):
        assert fast_cl.match_or_add(s, None, add=False) == exp


def test_nnls_weighted_exact_on_sparse_dense_and_collinear():
    rng = np.random.RandomState(123)
    for shape, density in [((20, 10), 1.0), ((60, 120), 0.12), ((40, 0), 1.0)]:
        A = (rng.rand(*shape) < density).astype(float) * rng.uniform(1.0, 50.0, shape)
        b = rng.uniform(0.0, 100.0, shape[0])
        b[rng.rand(shape[0]) < 0.3] = 0.0
        w = 1.0 / np.sqrt(np.maximum(16.0, b + (0.2 * b) ** 2))
        x_ref, r_ref = _orig_nnls(A, b, w)
        x_opt, r_opt = cohort._nnls_weighted(A, b, w)
        assert abs(r_opt - r_ref) < 1e-8
        if shape[1] > 0:
            assert np.max(np.abs(x_opt - x_ref)) < 1e-8

    # Collinear columns and zero target
    A = (rng.rand(50, 80) < 0.15).astype(float) * rng.uniform(5.0, 80.0, (50, 80))
    A[:, 5] = A[:, 0]
    A[:, 6] = A[:, 0] + A[:, 1]
    b = rng.uniform(10.0, 200.0, 50)
    w = 1.0 / np.sqrt(16.0 + (0.18 * b) ** 2)
    _, r_ref = _orig_nnls(A, b, w)
    _, r_opt = cohort._nnls_weighted(A, b, w)
    assert abs(r_opt - r_ref) < 1e-8

    x_z, r_z = cohort._nnls_weighted(A, np.zeros(50), w)
    assert r_z == 0.0 and np.all(x_z == 0.0)


def test_high_diversity_no_copresent_fast_and_thread_invariant():
    """Hyper-variable gene family with 80 samples and ~240 unique flanks on each side
    (co_present < 2 everywhere) must skip wasted NNLS and finish quickly, with
    identical results for threads=1 and threads=2."""
    rng = random.Random(77)
    gene = "".join(rng.choice("ACGT") for _ in range(1000))
    samples = []
    for i in range(80):
        l_paths, r_paths = [], []
        for k in range(3):
            fl = "".join(rng.choice("ACGT") for _ in range(550))
            fr = "".join(rng.choice("ACGT") for _ in range(550))
            l_paths.append({"id": "c%d" % (k + 1), "flank": fl, "flank_len": 550, "context_len": 150,
                            "depth": 30.0 + k * 5, "shared_end": 0, "gene_depth": 100.0, "depth_bins": [30.0] * 5})
            r_paths.append({"id": "c%d" % (k + 1), "flank": fr, "flank_len": 550, "context_len": 150,
                            "depth": 30.0 + k * 5, "shared_end": 0, "gene_depth": 100.0, "depth_bins": [30.0] * 5})
        samples.append({
            "sid": "S%03d" % i,
            "genes": {"lnuC_1": gene},
            "doc": {"genes": {"lnuC_1": {
                "left": {"read_depth": 100.0, "paths": l_paths},
                "right": {"read_depth": 100.0, "paths": r_paths},
            }}},
            "reads_tsv": None,
        })

    t0 = time.time()
    p1 = cohort.cohort_pair_all(samples, W=500, threads=1)
    dt = time.time() - t0
    assert dt < 5.0, "High-diversity Part A took too long: %.2fs" % dt
    p2 = cohort.cohort_pair_all(samples, W=500, threads=2)
    assert p1 == p2

    r1 = cohort.cohort_deconvolve_all_repeats(samples, threads=1)
    r2 = cohort.cohort_deconvolve_all_repeats(samples, threads=2)
    assert r1 == r2
