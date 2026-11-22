"""How far a candidate context is certain: the rules that read the records of awv.

A candidate was built by following reads; that says what CAN lie next to the gene, not what does. Where reads of another
locus come onto the path (a repeat, a shared element) the path may go on along that locus. The rules below vouch for a
position x of a candidate only on evidence that molecules which carry the gene reach x:

  G   two or more matching templates (a read, or a read pair with both mates matching) begin inside the gene and cover x
  J   another locus is on the path, and 8 or more templates that begin before it came on cover x
  B   budget: the matching reads next to the gene, minus every read seen leaving the path since (charged in the measure
      in which the leaving reads' mates lie on the path towards the gene), still number 8 or more, and stand three
      standard deviations clear of the chance in what was charged (0.4)
  F   flow: more matching reads cover x than all loci that came on can supply, allowing for chance in both numbers and
      for the supply being known only to about a tenth. What came on is what the reads seen coming on stand for, and at
      least what the depth has risen above its lowest level since the gene (beyond chance): a locus that comes on softly
      is seen by few of its reads (0.4; the faults of 0.3.0 on block E). And it never finds more of the gene's own reads
      than the budget allows at most: what was next to the gene, minus the reads charged as leaving, plus all that came on

and in any case only where two reads carry the 31 bases that end at x, where reads carry the 61 bases that end at x if
the depth on both sides says they should, and where the place is not one that reads do not cross (a hole in the
coverage with reads piled up at its edges). The verified length is the distance from the gene to the first position
that nothing vouches for.

The rules are run on two sets of records of the same candidates (awv.c): the strict measure, in which only reads without
any difference count, and the tolerant measure (patches=True), in which reads with a few single differences, small
patches of differing bases or one small insertion or deletion count as well, a small patch is no place where a read
leaves, and the records of reads that leave or come on carry a class: 1 = the read only passes a small patch (left out
here), 2 = too near the read's end to tell (counted in the measure of class 0 against class 1 at that place).
Under the tolerant measure a read that goes on along a version a few per cent away is not seen leaving anywhere, so
there the budget (B) is also charged with what the DEPTH says has left (all that was on the path next to the gene, plus
all that came on, minus what is on the path here, beyond chance), and it does not vouch among a crowd, where other loci
on the path outnumber it more than twice and the depth could not show its reads going.
A position is certain if it is vouched for under either measure (cli.py takes the longer of the two lengths).

How many reads a cluster of leaving or coming-on reads stands for: all reads seen, divided by the share of reads over a
place that can show it; and, where that says more, the reads with 70 bases or more beyond (before) the place divided
by their share: a version a few per cent away shows itself only in such reads (0.4).

History: the rules are those of the fifth-round prototype (benchmarks/dev-tools/experiments/frozen6_draft/vcall15.py with
--nomate --tclean --res3 --bbridge --sumfa --jmate --lmate --ball --aguard --gap 0.2 --rag --fq 0.1, policy r8), which
0.3.0 carried unchanged; 0.4 adds the rise term of the flow rule and the tolerant measure
(benchmarks/dev-tools/experiments/tol/NOTE).
"""
import numpy as np

NARR = 20                 # arrays per candidate in the .i32 file of awv
A_EXACT, A_S31, A_S61, A_ENTER, A_LEAVE, A_SOFT, A_RAGR, A_RAGL, A_PATH = 1, 2, 3, 6, 7, 14, 17, 18, 19
BIG = 1 << 29
# upper and lower 95 % limits of a Poisson mean, for counts 0, 1, 2, ...
UCB95 = (3.0, 4.75, 6.3, 7.76, 9.16, 10.52, 11.85, 13.15, 14.44, 15.71, 16.97)
LCB95 = (0.0, 0.05, 0.36, 0.82, 1.37, 1.97, 2.61, 3.29, 3.98, 4.70, 5.43, 6.17, 6.92, 7.69, 8.46, 9.25, 10.04, 10.83, 11.63, 12.44, 13.25)

DELTA = 15                # a template "begins before" a place if it begins this many bp before it
Z = 2.5                   # flow: own reads must exceed Z standard deviations of the depth
ZR = 4.0                  # flow: a rise of the depth above its lowest level, beyond what came on seen and ZR standard deviations, came on unseen
ZC = 2.0                  # tolerant measure, budget: what the depth says has left, less ZC standard deviations (of the depth and of what came on)
KAPPA = 2.0               # tolerant measure, budget: vouches only where what other loci may have on the path is at most KAPPA times the budget
FCAP = True               # flow: held to what the budget allows at most (see where it is used)
BSIG = 3.0                # budget: must stand BSIG standard deviations clear of the chance in what was charged to it
FARB = 70                 # a cluster of leaving (coming-on) reads is also counted from its reads with FARB bases or more beyond (before) the
                          # place (for reads of 150 bp; in proportion for shorter ones)
SOFTW = 300               # strict measure, budget: reads that leave the path in ONE base up to SOFTW bp before a place where reads are
                          # seen leaving are charged with them (for reads of 150 bp; in proportion for shorter ones; 0: not at all),
SOFTF = 0.25              # where another locus may be on the path with at least SOFTF of the depth next to the gene
GINSIDE = True            # pairs "from the gene" begin before any place in the gene part where another locus was seen coming on,
GMARGIN = 60              # by this many bases (the stretch over which the record program judges what a read is like before it comes on)
FQ = 0.1                  # flow: the supply of loci that came on is known to this share of itself
TJ = 8                    # templates from before a join that vouch beyond it
COHD = 15.0               # depth from which two reads must carry every 61-bp window
BMIN = 8.0                # budget: reads of the gene's own that must be left
ZB = 1.0
LMIN = 2                  # reads that make a cluster of leaving reads
UNEV = 0.20               # uneven depth: the reads of ONE locus are not spread evenly along it. Between two places their depth
                          # differs by more than chance in counts: by 9 to 19 % on long real contexts (20 Nov 2026), taken as UNEV.
                          # The budget compares a count at one place (the gene's reads next to the gene) with counts at others
                          # (the reads seen leaving): what is left of it must stand clear of UNEV times what was charged, not only
                          # of chance in counts. And under the tolerant measure the budget is in effect the depth here less
                          # what other loci supply, a supply counted where they came on: it is first taken down by UPOW * UNEV
                          # of all that other loci brought on (the gene's own swing cancels out of that difference). Made-up
                          # cases with depth swinging by 20 % along every molecule: 6 wrong claims in 25,740 without this, 1
                          # with it. (The same allowance in the flow rule and its cap changed no case at 20 % and cost four
                          # times the length; it is not made.)
UPOW = 2.0
FCAPU = False             # ... the same allowance on the flow rule's cap under the tolerant measure: OFF. It closes a hole the made-up
                          # cases show (see where it is used: 2 wrong claims in 25,740 at a 20 % swing), but takes the rarer of
                          # two contexts of a gene from certain for 10 kb to 6 kb in the self-test's "mixed" case (80x and 50x):
                          # the price was judged too high. 20 Nov 2026.
UBLOCK = 300              # what was charged at places UBLOCK bases or more apart errs independently (see where it is used)
BRIDGE = 100              # a stretch that nothing vouches for is bridged when it is at most BRIDGE bases long, reads carry the path
SUSTAIN = 300             # all through it, and beyond it the flow rule, pairs from the gene or the budget vouch again for SUSTAIN
                          # bases without a break (see where it is used)
BRQ = 0.2                 # ... and only if what was charged to the budget in the stretch is at most BRQ of what the flow rule took
                          # for the gene's own reads before it: where most of the gene's reads are seen leaving, it is a fork
AFTERG = True             # only reads that leave AFTER the gene end are charged to the budget. The prototypes also took those that
                          # leave in the gene's last 20 bases; harmless then, since the budget starts at the gene end and what
                          # was charged before cancels out of it. But three later terms read the SUM of all charges: the
                          # budget's bar (chance in what was charged), the flow rule's cap (what was next to the gene less
                          # what was charged) and the tolerant measure's charge by the depth (the larger of what the depth
                          # says has left since the gene and what was charged). Where relatives of the gene leave in its last
                          # bases (E-G2, erm(F): 1,994 reads leave there and 393 go on) the cap was negative, the bar a
                          # multiple of the budget, and the tolerant budget was never charged by the depth. 20 Nov 2026.
NEAR = 100                # ... in the last NEAR bases every group of reads coming on is a join; deeper in, see where it is used
INSIDE = 300              # reads coming on are joins from the last 300 bases of the gene part on: the stretch that the proposer reads
                          # off the sample (0.3.0 looked at the last 100 only; a candidate that went over to another locus 100 to
                          # 300 bases before the gene end was then judged as if that locus's reads were the gene's: known-answer
                          # cases from real genomes, 19 Nov 2026)
REACH = 400               # how far beyond a join a read pair that begins before it reaches
ATTMIN = 8
RESOPEN = 0.1
LAMFRAC = 0.1
TFRAC = 0.1
SOFTK, SOFTADD = 2.0, 10.0
BEXP, BT0 = 6.0, 8
RHO = 0.8                 # share of the gene's own leaving reads whose mate lies on the path, at a true fork
MROOM = 700               # a join this near the candidate's end: its reads' mates would lie beyond the end
SDET = 0.5
GAP, GAPQ = 0.2, 0.1      # a place reads do not cross: under GAP of the depth carry it, under GAPQ of the depth leave or come on
CUT = 30                  # bases taken off the end of every verified part before it is called certain
DEPTH_BIN = 100           # bin width (bp) for compact per-candidate flank depth profiles exported to flanks.json / GFA


def lcb(n):
    n = int(n)
    return LCB95[n] if n < len(LCB95) else n - 1.645 * n ** 0.5


def ucb(n):
    n = int(round(n))
    return UCB95[n] if n < len(UCB95) else n + 1.645 * n ** 0.5 + 1.5


def ucbf(v):                                              # the same upper limit for a count that is not whole
    v = max(0.0, float(v))
    return float(np.interp(v, np.arange(len(UCB95)), UCB95)) if v < len(UCB95) - 1 else v + 1.645 * v ** 0.5 + 1.5


def _ucb_array(v):
    return np.where(v <= 10, np.interp(v, np.arange(len(UCB95)), UCB95), v + 1.645 * np.sqrt(v) + 1.5)


def read_length_constants(read_len):
    """What depends on the read length. A read shows that it comes onto the path, or leaves it, only if 40 of its bases
    lie on the path and 12 beyond: (L - 52) / L of the reads over a place can show it. The rules were shaped and held
    out with reads of 150 bp and longer; for those the constants are those of 150 bp (for longer reads that is the
    cautious side). Shorter reads get their own."""
    L = max(60, min(150, int(read_len)))
    return {"L": L, "eff": round((L - 52) / float(L), 2), "show": float(L - 52), "tail": float(L - 12), "stretch_max": L - 25,
            "far": min(FARB, int(round(FARB * L / 150.0)))}


def _med41(a):
    """Running median of width 41 (w = 20, edge-padded) on a 1-D int32 array -> float64."""
    v = np.lib.stride_tricks.sliding_window_view(np.pad(a, (20, 20), mode="edge"), 41)
    return np.partition(v, 20, axis=1)[:, 20].astype(np.float64)


def _group_by_cand(rows):
    if len(rows) == 0:
        return rows, np.zeros(1, dtype=np.int64)
    rows = rows[np.argsort(rows[:, 0], kind="stable")]
    cuts = np.searchsorted(rows[:, 0], np.arange(0, int(rows[-1, 0]) + 2))
    return rows, cuts


def _load_records(pre):
    with open(pre + ".idx") as f:
        idx_rows = []
        for line in f:
            name, ln, off = line.rstrip("\n").split("\t")
            idx_rows.append((name, int(ln), int(off), int(name.split("|")[3])))
    with open(pre + ".i32", "rb") as f:
        blob = f.read()
    ent, ecuts = _group_by_cand(np.fromfile(pre + ".ent", dtype=np.int32).reshape(-1, 8))
    lev, lcuts = _group_by_cand(np.fromfile(pre + ".lev", dtype=np.int32).reshape(-1, 8))
    tpl, tcuts = _group_by_cand(np.fromfile(pre + ".tpl", dtype=np.int32).reshape(-1, 3))
    return idx_rows, blob, ent, ecuts, lev, lcuts, tpl, tcuts


def verify_gene(pre, read_len=150, detail=None, patches=False):
    """The rules on one gene's records. Per candidate: "verified", how far it is certain as a context of the gene as given
    (reads of other loci coming on are looked for over the gene's last INSIDE bases); and "verified_end", how far it is
    certain as a context of a gene that ENDS like the one given (they are looked for over its last NEAR bases only): where
    another sequence comes onto the gene's last NEAR to INSIDE bases (a relative that shares them, a fragment of the gene
    elsewhere, an insertion in the gene), the second can be longer; "shared_end" is how many of the gene's last bases lie
    after the deepest such place (0: none)."""
    loaded = _load_records(pre)
    out = _verify(loaded, read_len, detail, patches, INSIDE, None)
    again = {r["name"] for r in out if r["shared_end"]}
    if again:
        end = {r["name"]: r for r in _verify(loaded, read_len, None, patches, NEAR, again)}
        for r in out:
            e_ = end.get(r["name"])
            if e_ is not None and e_["verified"] > r["verified"]:
                r["verified_end"], r["stop_end"] = e_["verified"], e_["stop"]
    return out


def _verify(loaded, read_len, detail, patches, inside, only):
    """loaded: pre-read records of one gene -> one dict per candidate, in the order of the candidate file:
    name (as in the file), flank_len, verified (bp from the gene, before the cut), joins, open_joins, depth, stop.
    patches: the records were made under the tolerant measure (pipeline.TOLERANT)"""
    idx_rows, blob, ent, ecuts, lev, lcuts, tpl, cuts = loaded
    BUB = bool(patches)
    rc = read_length_constants(read_len)
    eff, rl, farb = rc["eff"], rc["L"], rc["far"]
    # Tolerant measure: a template is "from before" a place where a locus came on only if it begins a read length before it.
    # The locus's reads can lie on the path unseen that far before the place where they are SEEN coming on (they pass the
    # first small difference as a patch).
    dlj = rl if patches else DELTA
    softw = max(30, int(round(SOFTW * rl / 150.0)))
    out = []
    for ci, (name, ln, off, g) in enumerate(idx_rows):
        if only is not None and name not in only:
            continue
        A = np.frombuffer(blob, dtype=np.int32, count=NARR * (ln + 1), offset=off).reshape(NARR, ln + 1)
        dex, s31 = A[A_EXACT].astype(np.float64), A[A_S31]
        x = np.arange(ln + 1)
        ds = _med41(A[A_EXACT])
        dp = _med41(A[A_PATH])
        t = tpl[cuts[ci]:cuts[ci + 1]] if ci + 1 < len(cuts) else tpl[:0]
        st, en = t[:, 1], t[:, 2]
        st_all, en_all = st, en
        if len(st):
            # a pair whose unread middle holds a whole variable stretch (reads leave the path and reads come on again, each
            # cluster at least TFRAC of the depth) may belong to the other version of that stretch: it says nothing here
            thr = np.maximum(2.0, TFRAC * eff * ds)

            def pts(arr):
                clu = np.convolve(arr[:ln + 1], np.ones(4), mode="same")
                return np.nonzero(clu >= thr)[0]
            LP, EP = pts(A[A_LEAVE]), pts(A[A_ENTER])
            if len(LP) and len(EP):
                g0 = st + rl
                g1 = en - rl
                i = np.searchsorted(LP, g0, side="left")
                pl = LP[np.minimum(i, len(LP) - 1)]
                c1 = (i < len(LP)) & (pl < g1) & (g1 > g0)
                j = np.searchsorted(EP, pl, side="right")
                pe = EP[np.minimum(j, len(EP) - 1)]
                c2 = (j < len(EP)) & (pe < g1)
                keep_t = ~(c1 & c2)
                st, en = st[keep_t], en[keep_t]

        def thru(p, xs, dl=DELTA):                            # templates that begin before p - DELTA and cover xs
            ends = np.sort(en[st <= p - dl])
            return len(ends) - np.searchsorted(ends, xs, side="right")

        def thru_all(p, xs):                                  # the same, of all templates
            ends = np.sort(en_all[st_all <= p - DELTA])
            return len(ends) - np.searchsorted(ends, xs, side="right")

        # ---- reads that come on: groups of two or more within 3 bp are a locus that joins the path
        ec0 = ent[ecuts[ci]:ecuts[ci + 1]] if ci + 1 < len(ecuts) else ent[:0]
        e = ec0[ec0[:, 1] >= max(0, g - inside)]
        e = e[np.argsort(e[:, 1], kind="stable")].copy()
        egl = e[:, 4] // 10000                                # what each read is like before the place: 1 like the path, 2 foreign, 0 not told
        e[:, 4] = e[:, 4] % 10000
        ebub = np.zeros(0, dtype=np.int64)
        if BUB:
            ecls = e[:, 4] // 1000
            e = e.copy(); e[:, 4] = e[:, 4] % 1000
            ebub = np.sort(e[ecls == 1, 1])                   # where reads come back on after a small patch
            eamb = (ecls == 2)[ecls != 1]
            egl = egl[ecls != 1]
            e = e[ecls != 1]

        def bshare(pos_lo, pos_hi, n0, where):                # share of the reads at a place that are not such returns (add-one estimate)
            n1 = int(np.searchsorted(where, pos_hi + 3, side="right") - np.searchsorted(where, pos_lo - 3, side="left"))
            return min(1.0, (n0 + 1.0) / (n0 + n1 + 1.0))
        groups = []
        i = 0
        while i < len(e):
            j = i + 1
            while j < len(e) and e[j, 1] - e[j - 1, 1] <= 3:
                j += 1
            groups.append((i, j))
            i = j
        # what may be left of a locus beyond a pair's reach is the share of its reads whose mate lies ON the path further on
        # (matching, or with single differences only); groups with fewer than 4 such reads share one estimate
        jpool = 1.0
        sE = sN = 0
        for i, j in groups:
            if j - i >= 2 and int(e[i, 1]) + MROOM <= ln:
                fw_ = (e[i:j, 5] == -1) & (e[i:j, 6] != -3)
                if int(fw_.sum()) < 4:
                    sN += int(fw_.sum())
                    sE += int((fw_ & ((e[i:j, 6] >= 0) | (e[i:j, 6] == -2))).sum())
        if sN >= 8:
            jpool = min(1.0, ucb(sE) / float(sN))
        joins = []         # (from, reads, until, reads seen, what it can supply, may still be on the path, where it came on)
        came = []          # (where, how many reads it stands for): every locus that came on, once
        cameu = []         # the same with the upper limit of what it can supply
        islands = []       # (where a locus came on, where its own reads show it leaving, reads seen coming on)
        deep = []          # where another sequence comes onto the gene part deeper than its last NEAR bases
        for i, j in groups:
            n = j - i
            if n < 2:
                continue
            nf = float(n)
            p0 = int(e[i, 1])
            if p0 < g - NEAR:
                # Deeper in the gene part than its last NEAR bases: reads that are LIKE THE GENE before the place where they come
                # on (the record program marks them: at most 15 % of the up to 60 bases before it differ) carry the sample's
                # other allele of the gene; they are the gene's own, as they were for 0.3.0, which did not look this far in.
                # Reads that are foreign before the place are another locus that shares the gene's last bases. A group is a
                # join unless the reads like the gene are at least as many as the foreign ones.
                like, foreign = int((egl[i:j] == 1).sum()), int((egl[i:j] == 2).sum())
                if like >= 1 and like >= foreign:
                    continue
            if BUB:
                n2 = int(eamb[i:j].sum())
                if n2:
                    nf = (n - n2) + n2 * bshare(p0, int(e[j - 1, 1]), n - n2, ebub)
                if nf < 2:
                    continue
            if p0 < g - NEAR:
                deep.append(p0)
            ee, bb, a40, mt, mA = e[i:j, 2], e[i:j, 3], e[i:j, 4], e[i:j, 5], e[i:j, 6]
            gone = (bb - ee >= 12) & (a40 >= 2)
            u = BIG
            frr = 1.0
            if gone.any():
                # only reads that reach 15 bp beyond the place can show whether the locus stays; one that is not seen leaving
                # (it matches on, or differs in a single base) stays. Of the reads whose fate is seen, the share that stays:
                q = int(np.percentile(ee[gone], 90))
                stay = int(((bb >= q + 15) & ~(gone & (ee < q + 15))).sum())
                left_ = int((gone & (ee < q + 15)).sum())
                frr = (stay + 0.5) / (left_ + stay + 1.0)
                u = q + 5
            fa = ucb(nf) / eff
            if True:
                # of a locus that comes on at one place, (rl - 52) / rl of the reads over the place show it, and (rl - 110) / rl
                # of them with 70 bases or more before the place; one that comes on softly shows itself only in the latter
                nfar = int(((e[i:j, 1] - (bb - rl)) >= farb).sum())
                fa = max(fa, nfar * rl / float(max(10, rl - 40 - farb)))
            if BUB and nf < n:
                n = max(2, int(round(nf)))
            fwm = (mt == -1) & (mA != -3)
            nfm = int(fwm.sum())
            Em = int((fwm & ((mA >= 0) | (mA == -2))).sum())
            if p0 + MROOM > ln:
                frp = 1.0                                     # the mates would lie beyond the candidate's end
            elif nfm >= 4:
                frp = (Em + 0.5) / (nfm + 1.0)
            else:
                frp = jpool
            pr = min(p0 + REACH, ln)
            came.append((p0, nf / eff))
            cameu.append((p0, fa))

            def seg_(s0, s1, f_):
                if s1 <= s0:
                    return
                r_ = fa * f_
                joins.append((s0, nf, s1, n, r_, True if f_ >= 1.0 else r_ > RESOPEN * max(1.0, ds[min(ln, s0)]), p0))
            # the locus is on the path in full up to where its own reads show it leaving (u) or, failing that, up to a pair's
            # reach; beyond, what may be left of it is bounded by its reads that stay and by the mates of its reads
            if BUB:
                # Tolerant measure: what came on is supply for the rest of the candidate. Under this measure a relative's reads
                # flicker in and out of the count unseen (they pass the next cluster of differences as a patch, or begin
                # beyond it), so neither its reads seen leaving nor where their mates lie shows that the locus has left
                # (made-up cases with a diverged repeat, 19 Nov 2026: 49 wrong claims in 2,700 before this).
                seg_(p0, BIG, 1.0)
            elif u < pr:
                joins.append((p0, nf, u, n, fa, True, p0))
                islands.append((p0, u, n))
                seg_(min(u, ln), pr, frr)
                seg_(pr, BIG, min(frr, frp))
            else:
                seg_(p0, pr, 1.0)
                seg_(pr, BIG, min(frr, frp))
        if ln - g > 200:
            # a locus that comes on softly (its reads differ by single bases before, so no cluster of reads is seen coming on):
            # the depth rises far above anything since the gene, more than the joins seen in the last 200 bp can supply
            rmax = np.maximum.accumulate(np.where(x >= g, ds, 0.0))
            rlag = np.concatenate([np.full(min(150, ln + 1), rmax[g] if g <= ln else 0.0), rmax[:max(0, ln + 1 - 150)]])
            rlag = np.maximum(rlag, ds[g] if g <= ln else 0.0)
            rec = np.zeros(ln + 2)
            for p0, n, u, nraw, fa, isopen, anc in joins:
                if p0 <= ln:
                    rec[max(0, p0)] += fa
                    rec[min(ln + 1, max(0, p0) + 200)] -= fa
            rec = np.cumsum(rec)[:ln + 1]
            soft = (ds > SOFTK * (rlag + rec) + SOFTADD) & (x > g + 30)
            last = -10 ** 9
            for xs in np.nonzero(soft)[0]:
                if xs - last < 300:
                    continue
                last = int(xs)
                fa_s = float(ds[min(ln, xs + 150)] - rlag[xs])
                if fa_s > 0:
                    joins.append((max(g, int(xs) - 20), fa_s * eff, BIG, -1, fa_s, True, max(g, int(xs) - 20)))     # nraw -1 marks a soft join
                    came.append((max(g, int(xs) - 20), fa_s))
                    cameu.append((max(g, int(xs) - 20), fa_s))
        joins.sort(key=lambda z: z[0])
        # one upper limit for all the reads seen coming on that may still be on the path
        Nsum = np.zeros(ln + 2)
        Fsoft = np.zeros(ln + 2)
        pstar = np.full(ln + 1, -1, dtype=np.int64)   # the earliest place where a locus that may still be on the path came on
        for p0, n, u, nraw, fa, isopen, anc in joins:
            hi = min(u, ln + 1)
            if hi > p0:
                if nraw < 0:
                    Fsoft[p0] += fa
                    Fsoft[hi] -= fa
                else:
                    part = n * fa / (ucb(n) / eff)            # (with what its far reads say it can supply, where that is more)
                    Nsum[p0] += part
                    Nsum[hi] -= part
                if isopen:
                    seg = pstar[p0:hi]
                    seg[(seg < 0) | (seg > anc)] = anc
        Ns = np.maximum(0.0, np.cumsum(Nsum)[:ln + 1])
        Fa = np.where(Ns > 1e-9, np.where(Ns < 10, np.interp(Ns, np.arange(len(UCB95)), UCB95), Ns + 1.645 * np.sqrt(Ns) + 1.5), 0.0) / eff \
            + np.cumsum(Fsoft)[:ln + 1]
        # ---- what holds at every position whatever vouches: reads carry the path
        s61 = A[A_S61]
        dref = np.minimum(np.concatenate([np.full(100, 1e9), ds[:-100]]) if ln >= 100 else ds,
                          np.concatenate([ds[100:], np.full(100, 1e9)]) if ln >= 100 else ds)
        thin = (s61 < 2) & (dref >= COHD) & (dref < 1e8)      # the depth 100 bp to either side says reads should carry the 61-bp window
        # a place that reads do not cross: the reads with the path's 61-bp window are fewer than GAP times the depth on both sides
        # (upper 95 % limit of their count); it is not a fork (the reads seen leaving or coming on within the window, 15 bp to
        # spare, are fewer than GAPQ times that depth); and reads that stop or start at ONE position there with a few differing
        # bases at their very end (the pile within 1 bp, the largest inside the window) outnumber the reads that carry the window
        s61f = s61.astype(np.float64)
        evc = np.concatenate([[0.0], np.cumsum((A[A_ENTER] + A[A_LEAVE]).astype(np.float64))])
        ev = evc[np.minimum(ln + 1, x + 16)] - evc[np.maximum(0, x - 75)]
        rg = (A[A_RAGR] + A[A_RAGL]).astype(np.float64)
        pile = rg + np.concatenate([[0.0], rg[:-1]]) + np.concatenate([rg[1:], [0.0]])
        m67_ = np.pad(pile, (63, 3))
        for sh67_ in (1, 2, 4, 8, 16, 32):
            m67_ = np.maximum(m67_[:-sh67_], m67_[sh67_:])
        mx = np.maximum(m67_[:-3], m67_[3:])
        hole = (_ucb_array(s61f) < GAP * dref) & (ev < GAPQ * dref) & (dref < 1e8) & (x >= g + 61) & (mx > np.maximum(3.0, s61f))
        base = (s31 >= 2) & ~thin & ~hole
        # ---- F: flow
        clow = np.minimum.accumulate(np.where(x >= g, ds, 1e9))
        # What has come onto the path is at least what the depth has risen above its lowest level since the gene. Where that
        # rise is more than the loci seen coming on can supply, and more than chance (ZR standard deviations of the depth),
        # the rest came on unseen; it is counted as supply from there on (it cannot be seen leaving either).
        exc = np.where(x >= g, ds - clow - Fa - ZR * np.sqrt(ds), 0.0)
        FaF = Fa + np.maximum.accumulate(np.maximum(exc, 0.0))
        own = np.minimum(ds - FaF, clow)
        barF = np.maximum(3.0, Z * np.sqrt(ds + (FQ * FaF) ** 2))
        okF = (own >= barF) & base
        # ---- J and G: templates from before a join, or from the gene
        tau = np.zeros(ln + 1, dtype=np.int64)
        for p in np.unique(pstar[pstar >= 0]):
            xs = np.nonzero(pstar == p)[0]
            tau[xs] = thru(p, xs, dlj)
        # "From the gene" means from before any place in the gene part where another locus was seen coming on (they are looked
        # for in its last INSIDE bases): the reads of a locus that shares the gene's last bases begin in the gene part too, and
        # their mates lie along ITS flank (E. coli benchmark with the proposer that keeps reads through variant places,
        # 19 Nov 2026: a context that leaves the gene's locus 75 bp before the gene end was vouched for by 33 such pairs).
        gin = [p0_ for p0_, n_ in came if p0_ <= g - 10] if GINSIDE else []
        # ... and from GMARGIN bases before it: in the gene part a read counts as matching with up to 4 % of its bases differing,
        # so a relative's reads that begin a few bases before the place where they are seen coming on pass for matching
        # (E. coli x026 with the proposer that keeps reads through variant places: four such pairs, 255 bp).
        tg = thru(min(g - 10 + DELTA, min(gin) + DELTA - max(dlj, GMARGIN)) if gin else g - 10 + DELTA, x)
        # ---- B: budget
        lc0 = lev[lcuts[ci]:lcuts[ci + 1]] if ci + 1 < len(lcuts) else lev[:0]
        lv = lc0[(lc0[:, 3] >= g - 20) & (lc0[:, 1] == lc0[:, 2])]     # not those that had come on themselves
        lbub = np.zeros(0, dtype=np.int64)
        lamb = None
        if BUB:
            lcls = lv[:, 4] // 1000
            lbub = np.sort(lv[lcls == 1, 3])                  # where reads leave into a small patch and come back
            lv = lv[lcls != 1]
        if AFTERG:
            lv = lv[lv[:, 3] > g]
        bcaps = []
        if len(st):
            Lt_ = np.sort((en - st).astype(np.float64))
            cs_ = np.concatenate([[0.0], np.cumsum(Lt_)])
            tot_ = cs_[-1]

            def Rshare(d):                                    # share of the templates over a place that begin d bp or more before it
                k_ = np.searchsorted(Lt_, d, side="right")
                return ((tot_ - cs_[k_]) - d * (len(Lt_) - k_)) / max(tot_, 1.0)
            xg = min(ln, g + 50)
            tratio = float(ds[xg]) / max(1.0, float(np.sum((st <= xg) & (en > xg))))     # read depth per template over a place, at the gene
        else:
            def Rshare(d):
                return 0.0
            tratio = 1.0
        Lc = np.zeros(ln + 2)
        scale = np.ones(ln + 2)
        if len(lv):
            # a short stretch shared with a locus that came on and left again: the reads that leave at its end and begin inside
            # it may be that locus's (as many as its reads coming on lead one to expect) and are not charged
            keep = np.ones(len(lv), dtype=bool)
            for p0, u, nraw in islands:
                q = u - 5
                m = q - p0
                if m < 40 or m > rc["stretch_max"]:
                    continue
                att = (lv[:, 1] >= p0 - 2) & (np.abs(lv[:, 3] - q) <= 8)
                k = int(att.sum())
                if k and k <= 1.5 * nraw * max(0, m - 40) / rc["show"] + 3:
                    keep &= ~att
                    # ... but they may be the gene's own. What is known to stay are the templates from before the stretch that go
                    # on beyond it: the budget after the stretch is at most what they stand for (lower 95 % limit, scaled by
                    # the share of templates that can show it)
                    Sb = int(thru(p0, np.array([min(ln, q + 6)]))[0])
                    d_ = (q + 6) - (p0 - DELTA)
                    bcaps.append((min(ln, q), lcb(Sb) / max(Rshare(d_), 0.02) * tratio))
                    lo_, hi_ = max(0, q - 8), min(ln + 1, q + 9)
                    scale[lo_:hi_] = np.maximum(scale[lo_:hi_], rc["show"] / max(20.0, rc["tail"] - m))   # reads from before the stretch that leave here: fewer can show it
            lv = lv[keep]
        sp_ = np.sort((en_all - st_all)[(en_all - st_all) > rl + 10]) if len(st_all) else np.zeros(0)       # spans of the pairs on this candidate

        def wlaid(p):                                         # share of pairs short enough for the mate of a read at p to lie on the candidate
            return 1.0 if len(sp_) == 0 else max(0.05, float(np.searchsorted(sp_, p, side="right")) / len(sp_))
        clusters = []
        Mown = Uprev = cumnl = 0.0
        lv = lv[np.argsort(lv[:, 3], kind="stable")]
        softv = A[A_SOFT].astype(np.float64).copy() if (SOFTW > 0 and not BUB) else None
        softbg = 0.0                                          # reads that differ in one base by a read error: per read and base, the
        if softv is not None and ln - g >= 400:               # candidate's own rate (median over windows of 200 bp)
            softbg = float(np.median([float(softv[g + 200 * k_:g + 200 * (k_ + 1)].sum()) / max(1.0, float(dp[g + 200 * k_:g + 200 * (k_ + 1)].sum()))
                                      for k_ in range((ln - g) // 200)]))
        i = 0
        while i < len(lv):
            j = i + 1
            while j < len(lv) and lv[j, 3] - lv[j - 1, 3] <= 3:
                j += 1
            nl = float(j - i)
            if BUB:
                n2 = int((lv[i:j, 4] // 1000 == 2).sum())
                if n2:
                    nl = (j - i - n2) + n2 * bshare(int(lv[i, 3]), int(lv[j - 1, 3]), j - i - n2, lbub)
            if nl >= LMIN:
                mt, mAl = lv[i:j, 5], lv[i:j, 6]
                q = min(ln, int(lv[i, 3]))
                frac = None
                ps = int(pstar[q])
                if ps >= 0 and q - ps <= REACH and q + 6 <= ln:
                    lam = int(np.sum((mt == -1) & (mAl >= 0) & (mAl < ps - DELTA)))    # they leave, and their mates begin before the other locus came on
                    tq = int(thru(ps, np.array([q + 5]), dlj)[0])                    # templates from before it that go on
                    # only where the gene's own reads are not seen leaving: then the reads that leave are another locus's
                    if tq + lam >= ATTMIN and lam <= max(1.0, LAMFRAC * tq):
                        frac = ucb(lam) / (ucb(lam) + tq)
                # a read that leaves is charged as one of the gene's own in the measure in which leaving reads have their mate
                # ON the path towards the gene, given the share of pairs whose mate can lie on the candidate at that place
                rvm = (mt == -1) & (mAl != -3)
                Rm = int(rvm.sum())
                Xm = int((rvm & ((mAl >= 0) | (mAl == -2))).sum())
                own_ = nl if Rm == 0 else min(nl, nl * Xm / (float(Rm) * wlaid(q) * RHO))
                Mown += own_
                cumnl += nl
                Unew = min(cumnl, ucbf(Mown))                 # one upper limit for the running sum
                chg = max(0.0, Unew - Uprev) / eff * float(scale[q])
                if nl > 0:                                # ... and of reads that leave: the lower limit of those far beyond, where it says more
                    nfar = int(((lv[i:j, 1] + rl - lv[i:j, 3]) >= farb).sum())
                    chg *= max(1.0, (lcb(nfar) * rl / float(max(10, rl - 40 - farb))) / (nl / eff))
                if softv is not None and nl > 0 and float(Fa[q]) >= SOFTF * max(0.0, float(ds[g] - Fa[g])):
                    # Strict measure: reads that leave the path in one base (the first of their differences) before a place
                    # where reads are seen leaving are the same reads going. Where the gene's molecules go off along a relative
                    # 7 to 12 % away, most of their reads leave like that, up to 200 bp before the first is SEEN leaving, and
                    # what was left of the budget passed every bar (made-up cases, 19 Nov 2026: 25 wrong claims in 23,040, most
                    # at 240x and more). They are counted less what read errors alone give in the window, each once. Only where
                    # another locus is on the path: without one the depth falls with the gene's reads, and the budget is
                    # never more than the depth. This does not close the hole (2 of 8,640 at 240x and 480x were left); what
                    # does, a check of the budget by the depth, takes a fifth to a half of all strict claims on gut samples.
                    lo_, hi_ = max(g, int(lv[i, 3]) - softw), min(ln, int(lv[j - 1, 3]) + 10)
                    if hi_ >= lo_:
                        seg = softv[lo_:hi_ + 1]
                        bg_ = softbg * float(dp[lo_:hi_ + 1].sum())
                        tot_ = float(seg.sum())
                        use_ = max(0.0, tot_ - 1.5 * bg_ - 2.0 * bg_ ** 0.5)
                        if tot_ > 0 and use_ > 0:
                            seg *= 1.0 - use_ / tot_
                            chg *= 1.0 + use_ / nl
                Uprev = max(Uprev, Unew)
                clusters.append((q, chg, frac, ps))
            i = j
        c0 = max(0.0, float(ds[g] - Fa[g]))
        if gin:
            # A locus that came on inside the gene part may all be gone again at the gene end (its reads end where the two
            # part), yet it is taken off the depth there as if it stayed. The gene's own reads at the gene end are at least
            # the pairs that come from before that place and cover the gene end: three reads in four of such templates.
            c0 = max(c0, 0.75 * float(tg[g]))
        # never more than the depth at any place passed, and what left since: O(x) = min over x' <= x of depth(x') - leaves(x'..x]
        cap = np.where(x > g, ds, 1e9)
        cap[g] = c0
        for q_, v_ in bcaps:
            cap[q_] = min(cap[q_], v_)
        for q, chg, frac, ps_ in clusters:
            if frac is not None:
                Lcum = np.cumsum(Lc)[:ln + 1]
                Oq = float((np.minimum.accumulate(cap + Lcum) - Lcum)[max(g, q - 1)])
                # the reads that leave are written off to the locus that came on at ps_ only if no OTHER locus of at least SDET
                # times the budget was on the path before it: among such a crowd the gene's own reads are too few to tell
                older = 0.0
                for p0, n, u, nraw, fa, isopen, anc in joins:
                    if p0 <= q < min(u, ln + 1) and anc < ps_ - DELTA:
                        older = max(older, (fa * eff) if nraw < 0 else (nraw / eff) * min(1.0, fa / (ucb(n) / eff)))
                if not older >= SDET * max(Oq, 1e-9):
                    chg = min(chg, max(0.0, frac * Oq))
            Lc[q] += chg
        Lcum = np.cumsum(Lc)[:ln + 1]
        # Uneven depth: what was charged at places UBLOCK bases or more apart errs independently (the depth's swing at one
        # place says nothing about the next), so the allowance below is UNEV times the root of the sum of squares of the
        # charges, block by block
        Lseen = Lcum.copy()
        bid = np.maximum(0, np.arange(ln + 1) - g) // UBLOCK
        nb_ = int(bid[-1]) + 1
        bstart = np.array([0] + [g + b_ * UBLOCK for b_ in range(1, nb_)])
        bend = np.array([min(ln, g + (b_ + 1) * UBLOCK - 1) for b_ in range(nb_)])
        before_ = np.where(bstart > 0, Lseen[np.maximum(bstart - 1, 0)], 0.0)
        Sb_ = Lseen[bend] - before_
        Uq = np.concatenate([[0.0], np.cumsum(Sb_ ** 2)])[:nb_][bid] + (Lseen - before_[bid]) ** 2
        if BUB:
            # Tolerant measure. A read that goes on along a version a few per cent away has no place where it is SEEN leaving
            # (its clusters of differences pass for patches); it only stops matching. So the budget is charged at least with
            # what the depth says has left since the gene: all that was on the path next to the gene, plus all that came on
            # (seen, or shown by a rise of the depth), minus what is on the path here; less ZC standard deviations of the
            # depth and of the count of what came on. The depth is that of reads on the path for as long as each is on it.
            jc = np.zeros(ln + 2)
            for p0_, fa_ in came:
                if p0_ > g:
                    jc[min(ln + 1, p0_)] += fa_
            Jc = np.cumsum(jc)[:ln + 1] + (FaF - Fa)
            lcons = np.where(x > g, dp[g] + Jc - dp - ZC * np.sqrt(np.maximum(dp, 1.0) + Jc / eff), 0.0)
            lcons = lcons * min(1.0, float(ds[g]) / max(1.0, float(dp[g])))
            Lcum = np.maximum.accumulate(np.maximum(Lcum, lcons))
        O = np.minimum.accumulate(cap + Lcum) - Lcum
        Oall = O
        if BUB:
            # Uneven depth, tolerant measure. What the budget keeps here is the depth less the other loci's supply, and that
            # supply was counted where they came on; where their depth swings up along the path, the surplus is theirs and
            # passed for the gene's own (made-up case, 20 % swing: the gene's 30x gone at a fork, another locus at 30x counted
            # in a trough, 20 reads "left"). The budget is taken down by UPOW * UNEV of all that other loci have brought on
            # before its bars are asked. (The gene's own depth is on both sides of the difference and needs no allowance.)
            O = O - UPOW * UNEV * (float(Fa[g]) + Jc)
        okB = (O >= max(BMIN, ZB * c0 ** 0.5)) & base & (x >= g)
        # The charges are counts of reads seen leaving. Where (nearly) all of the gene's reads leave at a fork while another
        # locus keeps the depth up, what is left of the budget is chance around zero; 8 reads were no bar against that
        # (made-up cases, 19 Nov 2026: 8 wrong claims in 1,800). The budget must stand clear of the chance in its charges.
        # ... and of the unevenness of the depth between the place where the budget starts and the places where it was charged
        # (under the tolerant measure also in what the depth says has left beyond the reads seen)
        okB &= O >= BSIG * np.sqrt(np.maximum(Lcum, 0.0) / eff + UNEV ** 2 * (Uq + np.maximum(0.0, np.maximum(Lcum, 0.0) - Lseen) ** 2))
        if BUB:
            # ... and among a crowd (other loci on the path more than KAPPA times the budget) the depth cannot show the gene's
            # own reads going: their number is within chance of the crowd's. There the tolerant budget does not vouch.
            okB &= FaF <= KAPPA * np.maximum(O, 0.0)
        O = Oall
        # ---- a break in the path: templates from before a place where a locus came on should reach beyond it
        injoin = pstar >= 0
        broken = ln + 1
        if len(st):
            Lt = np.sort((en_all - st_all).astype(np.float64))
            csum = np.concatenate([[0.0], np.cumsum(Lt)])
            tot = csum[-1]

            def R(d):                                         # share of the templates over a place that reach d bp beyond it
                k = np.searchsorted(Lt, d, side="right")
                return ((tot - csum[k]) - d * (len(Lt) - k)) / max(tot, 1.0)
            for p0, n, u, nraw, fa, isopen, anc in joins:
                if p0 < g or (u < BIG and u != min(p0 + REACH, ln + 1)) or p0 + 45 > ln or p0 >= broken:
                    continue
                T0 = int(thru_all(p0, np.array([p0 + 5]))[0])
                if T0 < BT0:
                    continue
                dd = np.arange(40, 420, 20)
                ex = T0 * np.array([R(d + 20) for d in dd]) / max(R(20), 1e-9)
                okd = dd[(ex >= BEXP) & (p0 + dd <= ln)]
                if not len(okd):
                    continue
                dstar = int(okd[-1])
                if int(thru_all(p0, np.array([p0 + dstar]))[0]) <= 1:
                    xs_ = np.arange(p0, min(ln, p0 + dstar) + 1)
                    tt = thru_all(p0, xs_)
                    bad_ = np.nonzero(tt < 2)[0]
                    broken = min(broken, int(xs_[bad_[0]]) if len(bad_) else p0 + dstar)
        okG = (tg >= 2) & base
        okJ = injoin & (tau >= TJ) & base
        if FCAP:
            # The flow rule cannot find more of the gene's own reads than the budget allows at most: what was next to the gene,
            # minus the reads charged as leaving, plus what came on AND IS GONE (the charged reads may have been those; what
            # may still be on the path is taken off the depth by this rule, so it cannot explain the charges as well).
            # Without this the two rules could contradict each other and the flow rule win: all the gene's reads leave at a
            # fork and are charged, the other locus was not SEEN coming on (it came on softly, or through small patches the
            # tolerant measure lets pass), and its depth passes for the gene's own (E. coli benchmark, 19 Nov 2026: one context
            # wrong for 660 bp at 30x and for 3.5 kb at 60x; two strains at 15x and 60x: four contexts wrong for 2.9 kb).
            jcu = np.zeros(ln + 2)
            for p0_, fa_ in cameu:
                if p0_ > g:
                    jcu[min(ln + 1, p0_)] += fa_
            # ("what was next to the gene": every read on the path there, also those that leave within their own length and so
            # are not among the matching reads the budget starts from; they are among the reads charged)
            capF = max(c0, float(dp[g] - Fa[g])) + np.maximum(0.0, float(Fa[g]) + np.cumsum(jcu)[:ln + 1] - Fa) - Lcum
            if BUB and FCAPU:
                # (Off: FCAPU.) Uneven depth, tolerant measure: what the cap leaves is a count next to the gene less counts
                # further on, and is off by UNEV of what was charged (block by block, as in the budget's bar). Under this
                # measure a locus can come on unseen, so the cap is all that holds the flow rule there: where the gene's
                # molecules leave at a repeat and another locus comes on through it, a remainder within the swing passes for
                # the gene's own (made-up cases, 20 % swing, a relative at the gene end: 2 wrong claims in 25,740; none with
                # this; 12 against 9 at a 40 % swing). KNOWN WEAKNESS while this is off.
                capF = capF - Z * UNEV * np.sqrt(Uq + np.maximum(0.0, np.maximum(Lcum, 0.0) - Lseen) ** 2)
            okF = okF & (np.minimum(own, capF) >= barF)
        ok = okF | okG | okJ | okB
        if broken <= ln:
            ok &= x < broken
        bad = np.nonzero(~ok[g:])[0]
        V = min(int(bad[0]) if len(bad) else ln - g, ln - g)
        # Bridging. The flow rule, pairs from the gene and the budget are each evidence FOR THEIR POSITION that the gene's
        # molecules reach it, hence every position before it; so a short stretch that nothing vouches for, with reads carrying
        # the path all through it and one of those rules vouching again beyond it for SUSTAIN bases without a break, does not
        # end the certain part. (Pairs from before a join are left out: that rule leans on the place before the join being
        # certain.) Without this the certain part ended at places where another locus's reads lie on the path for a few
        # dozen bases and the flow rule cannot stand clear of them (B-A10, aph(6)-Id: vouched 0-2627 and 2655-9999).
        # The stretch is kept short and the evidence beyond it long because a rule that holds at one position by chance
        # somewhere along ten kilobases must not count.
        nbridge = 0
        if BRIDGE and V < ln - g:
            okd = okF | okG | okB
            bb_ = base.copy()
            if broken <= ln:
                okd &= x < broken
                bb_ &= x < broken
            while V < ln - g:
                a_ = g + V                                            # the first position nothing vouches for
                lim = min(ln, a_ + BRIDGE)
                nb_ = np.nonzero(~bb_[a_:lim + 1])[0]
                endb = a_ + int(nb_[0]) if len(nb_) else lim + 1      # reads carry the path up to here
                found = -1
                for b_ in range(a_ + 1, endb):
                    if okd[b_] and b_ + SUSTAIN <= ln + 1 and okd[b_:b_ + SUSTAIN].all() and bb_[b_:b_ + SUSTAIN].all():
                        found = b_
                        break
                if found < 0:
                    break
                # ... but not over a fork: in the stretches this is for, another locus's reads lie on the path and nothing
                # is charged to the budget (B-A10 aph(6)-Id, E-G1 tet(40): 0 charged); where most of the gene's reads are
                # seen leaving, the flow rule alone can be fooled beyond by the other locus's depth swinging up (made-up
                # case, 20 % swing: the certain part went from 1,489 to 1,844 bases with the fork at 1,500)
                w0_ = max(g, a_ - 26)
                if float(Lcum[min(ln, found)] - Lcum[w0_]) > BRQ * max(0.0, float(own[w0_])):
                    break
                bad2 = np.nonzero(~ok[found:])[0]
                V = min((found + (int(bad2[0]) if len(bad2) else ln + 1 - found)) - g, ln - g)
                nbridge += 1
        # why it ends there (information for the tables; no rule reads it)
        if V >= ln - g:
            why = "end_of_candidate"
        else:
            xv = g + V
            if hole[xv]:
                why = "reads_do_not_cross"
            elif s31[xv] < 2 or thin[xv]:
                why = "too_few_reads"
            elif broken <= xv:
                why = "pairs_do_not_cross"
            elif injoin[xv]:
                why = "another_locus_on_the_path"
            else:
                why = "own_reads_spent"
        g_dep = round(float(np.median(dex[max(0, g - 200):g])), 1) if g > 0 else 0.0
        d_bins = [round(float(np.median(dex[pos:min(ln + 1, pos + DEPTH_BIN)])), 1)
                  for pos in range(g, ln, DEPTH_BIN)] if ln > g else []
        out.append({"name": name, "flank_len": ln - g, "verified": V, "joins": len(joins),
                    "open_joins": sum(1 for p0, n, u, nraw, fa, isopen, anc in joins if u >= BIG and isopen),
                    "depth": float(np.median(dex[g:g + 200])) if ln > g else 0.0,
                    "gene_depth": g_dep, "depth_bins": d_bins, "stop": why,
                    "verified_end": V, "stop_end": why, "shared_end": (g - min(deep)) if deep else 0, "bridged": nbridge})
        if detail and detail in name:
            import sys
            sys.stderr.write("%s verified %d (%s)\n joins (flank bp, reads, on the path until, supply, open): %s\n" % (
                name, V, why, [(p0 - g, round(n, 1), (u - g) if u < BIG else "open", round(fa, 1), isopen) for p0, n, u, nraw, fa, isopen, anc in joins]))
            for xx in range(g, min(ln, g + 3000), 25):
                sys.stderr.write("  %5d depth %4.0f supply %5.1f low %4.0f join %5d tau %4d tg %3d budget %5.1f (start %.0f)  %s\n" % (
                    xx - g, ds[xx], Fa[xx], clow[xx], pstar[xx] - g if pstar[xx] >= 0 else -1, tau[xx], tg[xx], O[xx], c0,
                    " ".join(k for k, v in (("F", okF), ("G", okG), ("J", okJ), ("B", okB)) if v[xx])))
    return out
