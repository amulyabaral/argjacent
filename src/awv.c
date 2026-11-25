// awv.c (ARGjacent 0.5.0) - the record program: lays the reads kept near the genes on every candidate context of a gene
// and writes, per candidate position, what the reads say. Supports both standalone single-measure runs and single-pass dual
// strict + tolerant runs (-Z vt_jobs.tsv -q) sharing read loading, candidate 21-mer indexing, read seeding, and mismatch scans.
//
//   awv jobs.tsv reads.tsv                  the STRICT measure of 0.3.0: a read matches a candidate only without any difference in
//                                           the flank part; depth, templates and windows count matching reads only
//   awv -T 50 -B 1 -I 12 -g 30 -w 9 jobs.tsv reads.tsv   the TOLERANT measure: a read also matches with a few single differences (at most
//                                           1 + length/50 outside patches, 1 + length/12 differing bases in all), with small patches
//                                           (a group of differences of at most 40 bp and 12 differences from which the read comes
//                                           back for 15 bases), and with one insertion or deletion of up to 12 bases; the windows of
//                                           31, 61 ... bases that one read carries run on over a single difference only if no other
//                                           lies within 30 bases of it
//   awv -q -Z vt_jobs.tsv -T 50 -B 1 -I 12 -g 30 -w 9 v_jobs.tsv reads.tsv
//                                           single-pass DUAL mode used by argjacent: writes strict records to v_jobs.tsv prefixes
//                                           and tolerant records to vt_jobs.tsv prefixes in one pass; -q omits computing unused
//                                           prototype arrays (0, 4, 5, 8, 9..13) while keeping the 20-array .i32 layout intact.
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <pthread.h>
#define AWV_VERSION "0.5.1"
#define K 21
#define NW 4
#define NARR 20
#define QMIN 15
#define BIG 0x3fffffff
#define RET 15
#define MAXPL 256
static const int WS[NW] = {31, 61, 91, 121};
static const int PER = 0, MINON = 40, BOTH = 1, DW = 30, PERM = 100, INDEL = 60, HIST = 0;

typedef struct { char *name; uint8_t *s; int len, g; } Cand;
typedef struct { uint64_t km; int32_t pos; uint16_t c; } Post;
typedef struct { int c, d, st, i0, i1, A, B, ls, le, n, ng, cr, cl, a40, rgr, rgl, mm_off, ok; } BaseHit;
typedef struct { int c, d, st, A, B, ls, le, n, cr, cl, a40, ok, indel, exact, en, lv, enp, lvp, ee2, ls2, cr2, a402, bub, onp, sfp, rgr, rgl, match, A2, B2, nodense, i0, i1, lcls, ecls, npat, plo[4], phi[4], npd, nqq, ps1, pe1, ps2, pe2, gl; } Pl;

typedef struct {
    int tol, kd, maxi, bub, wiso, patchmax, patchdiff, totper, w100, quick;
    int32_t ***arr;                      /* arr[c][j] for c in 0..nc-1, j in 0..NARR-1 (NULL for skipped arrays when quick=1) */
    FILE *f_ent, *f_lev, *f_tpl;
    uint64_t n_tol, n_stitch, n_orph, n_lvo, n_eno, n_isl, n_bub, n_placed, n_exact, n_pairs;
    long isz[2001];
} VMode;

static int TOL = 0, KD = 3, MAXI = 0, BUB = 0, WISO = 12;
static int PATCHMAX = 40, PATCHDIFF = 12, TOTPER = 12, W100 = 0, QUICK = 0;
static int8_t code[256];

static __thread Cand *cd; static __thread int nc = 0;
static __thread Post *post; static __thread size_t npost = 0;
static __thread uint32_t *ht; static __thread uint64_t hmask;
static __thread uint8_t t_fs[1024], t_q[1024], t_r[1024], t_rq[1024];

static int cmp_post(const void *x, const void *y) {
    const Post *a = x, *b = y;
    if (a->km != b->km) return a->km < b->km ? -1 : 1;
    if (a->c != b->c) return a->c < b->c ? -1 : 1;
    return a->pos - b->pos;
}
static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static inline long find(uint64_t km) {
    uint64_t h = mix(km) & hmask;
    while (ht[h]) { if (post[ht[h] - 1].km == km) return (long)ht[h] - 1; h = (h + 1) & hmask; }
    return -1;
}
static inline void radd(int32_t *d, int a, int b) { d[a]++; d[b]--; }
static void low3(int32_t *b1, int32_t *b2, int32_t *b3, int A, int B, int s) {
    for (int x = A; x < B; x++) {
        if (s < b1[x]) { if (b3) b3[x] = b2[x]; b2[x] = b1[x]; b1[x] = s; }
        else if (s < b2[x]) { if (b3) b3[x] = b2[x]; b2[x] = s; }
        else if (b3 && s < b3[x]) b3[x] = s;
    }
}

static inline int arr_used(int quick, int j) {
    if (!quick) return 1;
    return (j >= 1 && j <= 3) || j == 6 || j == 7 || (j >= 14 && j <= 19);
}

/* Step 1 (shared across strict and tolerant): probe 21-mer hits directly on packed read `b`, and if any hit is found,
   unpack forward/reverse read buffers once and extract mismatches `mm[]` and exact stretches for every placement. */
static int seed_read(const uint8_t *b, int L, BaseHit *bh, int *mm_pool, int *any_out) {
    int np = 0, offs[128], no = 0;
    for (int o = 0; o + K <= L && no < 126; o += 11) offs[no++] = o;
    if (no && offs[no - 1] != L - K) offs[no++] = L - K;
    for (int st = 0; st < 2; st++) {
        for (int oi = 0; oi < no; oi++) {
            int off = offs[oi];
            uint64_t km = 0; int bad = 0;
            if (!st) {
                for (int i = 0; i < K; i++) { uint8_t x = b[off + i] & 7; if (x > 3) { bad = 1; break; } km = (km << 2) | x; }
            } else {
                for (int i = 0; i < K; i++) { uint8_t x = b[L - 1 - (off + i)] & 7; if (x > 3) { bad = 1; break; } km = (km << 2) | (3 - x); }
            }
            if (bad) continue;
            long p = find(km); if (p < 0) continue;
            for (size_t j = (size_t)p; j < npost && post[j].km == km; j++) {
                int c = post[j].c, d = post[j].pos - off, seen = 0;
                for (int t = 0; t < np; t++) if (bh[t].c == c && bh[t].d == d && bh[t].st == st) { seen = 1; break; }
                if (!seen && np < MAXPL) { bh[np].c = c; bh[np].d = d; bh[np].st = st; bh[np].ok = 0; np++; }
            }
        }
    }
    *any_out = 0;
    if (!np) return 0;
    for (int i = 0; i < L; i++) {
        uint8_t fb = b[i] & 7, fq = (b[i] & 8) ? 33 : 70;
        t_fs[i] = fb; t_q[i] = fq;
    }
    for (int i = 0; i < L; i++) {
        t_r[i] = t_fs[L - 1 - i] < 4 ? 3 - t_fs[L - 1 - i] : 4;
        t_rq[i] = t_q[L - 1 - i];
    }
    int mm_tot = 0, any = 0;
    for (int t = 0; t < np; t++) {
        Cand *C = &cd[bh[t].c]; int d = bh[t].d;
        const uint8_t *s = bh[t].st ? t_r : t_fs, *qq = bh[t].st ? t_rq : t_q;
        int i0 = d < 0 ? -d : 0, i1 = L < C->len - d ? L : C->len - d;
        bh[t].i0 = i0; bh[t].i1 = i1;
        if (i1 - i0 < 40) continue;
        int *mm = mm_pool + mm_tot, n = 0, ng = 0;
        for (int i = i0; i < i1; i++) {
            uint8_t a = s[i], cb = C->s[d + i];
            if (a != cb && a < 4 && cb < 4 && qq[i] >= 33 + QMIN) { mm[n++] = d + i; if (d + i < C->g) ng++; }
        }
        int A = d + i0, B = d + i1, ls = A, le = A, ps = A;
        for (int j = 0; j <= n; j++) {
            int pe = j < n ? mm[j] : B;
            if (pe - ps > le - ls) { ls = ps; le = pe; }
            ps = pe + 1;
        }
        if (le - ls < 31) continue;
        mm_tot += n;
        any = 1; bh[t].ok = 1; bh[t].A = A; bh[t].B = B; bh[t].ls = ls; bh[t].le = le;
        bh[t].n = n; bh[t].ng = ng; bh[t].mm_off = (int)(mm - mm_pool);
        bh[t].rgr = bh[t].rgl = -1;
        bh[t].cr = bh[t].cl = bh[t].a40 = 0;
        if (n > 0) {
            int jf = -1; for (int j = 0; j < n; j++) if (mm[j] >= B - 11) { jf = j; break; }
            if (jf >= 0) { int prev = jf > 0 ? mm[jf - 1] + 1 : A; if (mm[jf] - prev >= 40) bh[t].rgr = mm[jf]; }
            int jl = -1; for (int j = n - 1; j >= 0; j--) if (mm[j] < A + 11) { jl = j; break; }
            if (jl >= 0) { int next = jl + 1 < n ? mm[jl + 1] : B; if (next - (mm[jl] + 1) >= 40) bh[t].rgl = mm[jl]; }
            for (int j = 0; j < n; j++) {
                if (mm[j] >= le && mm[j] < le + 30) bh[t].cr++;
                if (mm[j] < ls && mm[j] >= ls - 30) bh[t].cl++;
                if (mm[j] >= le && mm[j] < le + 40) bh[t].a40++;
            }
        }
    }
    *any_out = any;
    return np;
}

static void eval_mode(VMode *M, int L, const BaseHit *bh, int np, const int *mm_pool, int any, Pl *pl) {
    const int tol = M->tol, kd = M->kd, maxi = M->maxi, bub = M->bub, wiso = M->wiso;
    const int patchmax = M->patchmax, patchdiff = M->patchdiff, totper = M->totper, w100 = M->w100, quick = M->quick;
    const int nw_lim = quick ? 2 : NW;
    for (int t = 0; t < np; t++) {
        memset(&pl[t], 0, sizeof(Pl));
        pl[t].c = bh[t].c; pl[t].d = bh[t].d; pl[t].st = bh[t].st;
        pl[t].i0 = bh[t].i0; pl[t].i1 = bh[t].i1;
        if (!bh[t].ok) continue;
        Cand *C = &cd[bh[t].c]; int32_t **ca = M->arr[bh[t].c];
        int A = bh[t].A, B = bh[t].B, n = bh[t].n, ng = bh[t].ng;
        pl[t].ok = 1; pl[t].A = A; pl[t].B = B; pl[t].ls = bh[t].ls; pl[t].le = bh[t].le; pl[t].n = n;
        pl[t].rgr = bh[t].rgr; pl[t].rgl = bh[t].rgl;
        pl[t].cr = bh[t].cr; pl[t].cl = bh[t].cl; pl[t].a40 = bh[t].a40;
        if (!quick && n * 10 <= B - A) radd(ca[0], A, B);
        if (n == 0) {
            /* Fast path for zero-mismatch placements (60-85 % of placed reads): exact match, no dense places, no events. */
            pl[t].lvp = pl[t].enp = pl[t].sfp = -1;
            pl[t].ee2 = B; pl[t].ls2 = A;
            pl[t].nqq = 1; pl[t].ps1 = A; pl[t].pe1 = B; pl[t].ps2 = A; pl[t].pe2 = B;
            pl[t].onp = 1; pl[t].nodense = 1;
            pl[t].A2 = A; pl[t].B2 = B; pl[t].exact = 1; pl[t].match = 1; M->n_exact++;
            radd(ca[1], A, B); radd(ca[19], A, B);
            if (!quick) { low3(ca[9], ca[10], NULL, A, B, A); low3(ca[11], ca[12], ca[13], A, B, A); }
            { int32_t r3[3] = { pl[t].c, A, B }; fwrite(r3, 4, 3, M->f_tpl); }
            for (int w = 0; w < nw_lim; w++) if (B - A >= WS[w]) { ca[2 + w][A + WS[w] - 1]++; ca[2 + w][B]--; }
            continue;
        }
        const int *mm = mm_pool + bh[t].mm_off;
        const uint8_t *rseq_ = bh[t].st ? t_r : t_fs;
        int d = bh[t].d;
        pl[t].lvp = pl[t].enp = -1; pl[t].cr2 = 0; pl[t].bub = 0; pl[t].sfp = -1;
        {
            static __thread uint8_t dn[1024], pf[1024], amb[1024];
            for (int j = 0; j < n; j++) dn[j] = pf[j] = amb[j] = 0;
            for (int j = 0; j < n; j++) { int k = 0; for (int i = j; i < n && mm[i] < mm[j] + DW; i++) k++; if (k >= kd) for (int i = j; i < j + k; i++) dn[i] = 1; }
            pl[t].npat = 0; pl[t].npd = 0; pl[t].lcls = pl[t].ecls = 0;
            if (bub) {
                for (int j = 0; j < n; ) {
                    int k = j, dense = dn[j];
                    while (k + 1 < n && mm[k + 1] - mm[k] < RET) { k++; dense |= dn[k]; }
                    if (dense) {
                        int lo = mm[j], hi = mm[k];
                        int small = (hi - lo + 1 <= patchmax) && (k - j + 1 <= patchdiff);
                        if (small && lo - A >= RET && B - 1 - hi >= RET && lo >= C->g) {
                            for (int i = j; i <= k; i++) pf[i] = 1;
                            pl[t].npd += k - j + 1;
                            if (pl[t].npat < 4) { pl[t].plo[pl[t].npat] = lo; pl[t].phi[pl[t].npat] = hi; pl[t].npat++; }
                        } else if (small) for (int i = j; i <= k; i++) amb[i] = 1;
                    }
                    j = k + 1;
                }
                if (pl[t].npat) {
                    for (int j = 0; j < n; j++) dn[j] = 0;
                    for (int j = 0; j < n; j++) {
                        if (pf[j]) continue;
                        int k = 0; for (int i = j; i < n && mm[i] < mm[j] + DW; i++) if (!pf[i]) k++;
                        if (k >= kd) for (int i = j; i < n && mm[i] < mm[j] + DW; i++) if (!pf[i]) dn[i] = 1;
                    }
                }
                if (w100 > 0) {
                    int changed = 0;
                    for (int j = 0; j < n; j++) {
                        if (mm[j] < C->g) continue;
                        int k = 0; for (int i = j; i < n && mm[i] < mm[j] + 100; i++) k++;
                        if (k >= w100) for (int i = j; i < j + k; i++) { if (pf[i]) changed = 1; dn[i] = 1; pf[i] = 0; }
                    }
                    if (changed) {
                        int w = 0; pl[t].npd = 0;
                        for (int q = 0; q < pl[t].npat; q++) {
                            int keep = 1, cnt = 0;
                            for (int i = 0; i < n; i++) if (mm[i] >= pl[t].plo[q] && mm[i] <= pl[t].phi[q]) { cnt++; if (!pf[i]) keep = 0; }
                            if (keep) { pl[t].plo[w] = pl[t].plo[q]; pl[t].phi[w] = pl[t].phi[q]; w++; pl[t].npd += cnt; }
                        }
                        pl[t].npat = w;
                    }
                }
            }
            int s = A, iso = 0, first_s = -1, first_e = -1, last_s = -1, last_e = -1, nq = 0, fiso = -1;
            int PERX = tol ? tol : PER;
            for (int j = 0; j <= n; j++) {
                if (j < n && !dn[j]) { if (!pf[j]) iso++; if (fiso < 0 && mm[j] >= C->g) fiso = mm[j]; continue; }
                int e = j < n ? mm[j] : B;
                if (e - s >= MINON && (PERX == 0 || iso <= 1 + (e - s) / PERX)) { if (!nq) { first_s = s; first_e = e; } last_s = s; last_e = e; nq++; if (pl[t].sfp < 0 && fiso >= 0) pl[t].sfp = fiso; }
                s = e + 1; iso = 0; fiso = -1;
            }
            pl[t].ee2 = B; pl[t].ls2 = A;
            if (nq == 1 || (nq > 1 && !BOTH)) {
                if (first_s > A && first_s - A >= 12) pl[t].enp = first_s;
                if (last_e < B && B - last_e >= 12) pl[t].lvp = last_e;
                if (pl[t].lvp >= 0) pl[t].ee2 = pl[t].lvp;
                if (pl[t].enp >= 0) pl[t].ls2 = pl[t].enp;
            } else if (nq > 1) {
                pl[t].lvp = first_e; if (BOTH == 1) pl[t].enp = last_s;
                if (first_s > A && first_s - A >= 12) pl[t].ls2 = first_s;
                if (BOTH == 1 && last_e < B && B - last_e >= 12) pl[t].ee2 = last_e;
                if (BOTH == 2 && first_s > A && first_s - A >= 12) { pl[t].enp = first_s; pl[t].ee2 = first_e; }
            }
            if (pl[t].lvp >= 0) for (int j = 0; j < n; j++) if (mm[j] >= pl[t].lvp && mm[j] < pl[t].lvp + 30) pl[t].cr2++;
            pl[t].nqq = nq; pl[t].ps1 = first_s; pl[t].pe1 = first_e; pl[t].ps2 = last_s; pl[t].pe2 = last_e;
            pl[t].bub = nq > 1;
            { int anyd = 0; for (int j = 0; j < n; j++) anyd |= dn[j]; pl[t].onp = !anyd && (PERM == 0 || n <= 1 + (B - A) / (tol ? tol : PERM)); pl[t].nodense = !anyd; }
            pl[t].a402 = 0; for (int j = 0; j < n; j++) if (mm[j] >= pl[t].ee2 && mm[j] < pl[t].ee2 + 40) pl[t].a402++;
            pl[t].gl = 0;
            if (pl[t].enp >= 0) {
                int lo = pl[t].enp - 60 > A ? pl[t].enp - 60 : A, len = pl[t].enp - lo, nb = 0;
                for (int j = 0; j < n; j++) if (mm[j] >= lo && mm[j] < pl[t].enp) nb++;
                if (len >= 20) {
                    pl[t].gl = nb * 100 <= 15 * len ? 1 : (nb * 100 > 30 * len ? 2 : 0);
                    if (pl[t].gl != 1) {
                        for (int sh = -60; sh <= 60 && pl[t].gl != 1; sh++) {
                            if (!sh) continue;
                            int cmpn = 0, same = 0;
                            for (int p = lo; p < pl[t].enp; p++) { int i = p - d, q = p + sh; if (i < 0 || i >= L || q < 0 || q >= C->len || rseq_[i] > 3 || C->s[q] > 3) continue; cmpn++; same += rseq_[i] == C->s[q]; }
                            if (cmpn >= 20 && same * 100 >= 85 * cmpn) pl[t].gl = 1;
                        }
                    }
                }
            }
            if (bub) {
                for (int j = 0; j < n; j++) {
                    if (pl[t].lvp >= 0 && mm[j] == pl[t].lvp && amb[j] && B - pl[t].lvp < patchmax + RET) pl[t].lcls = 2;
                    if (pl[t].enp >= 0 && mm[j] == pl[t].enp - 1 && amb[j] && pl[t].enp - A < patchmax + RET) pl[t].ecls = 2;
                }
            }
            if (HIST && (pl[t].lvp >= 0 || pl[t].enp >= 0)) fprintf(stderr, "H\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", pl[t].c, A, B, pl[t].enp, pl[t].lvp, nq, n);
        }
        int gov = (C->g < B ? C->g : B) - A; if (gov < 0) gov = 0;
        pl[t].A2 = A; pl[t].B2 = B; pl[t].match = 0;
        if (n == ng && ng <= gov / 25) { pl[t].exact = 1; pl[t].match = 1; M->n_exact++; }
        else if (tol && pl[t].nodense && ng <= gov / 25) {
            int A2 = A, B2 = B, nf = 0;
            for (int j = 0; j < n; j++) if (mm[j] >= C->g) { if (mm[j] < A + 12) A2 = mm[j] + 1; }
            for (int j = n - 1; j >= 0; j--) if (mm[j] >= C->g) { if (mm[j] >= B - 12) B2 = mm[j]; }
            for (int j = 0; j < n; j++) if (mm[j] >= C->g && mm[j] >= A2 && mm[j] < B2) nf++;
            if (B2 - A2 >= 40 && nf - pl[t].npd <= 1 + (B2 - A2) / tol && nf <= 1 + (B2 - A2) / totper) { pl[t].match = 1; pl[t].A2 = A2; pl[t].B2 = B2; M->n_tol++; }
        }
        if (pl[t].match) {
            int A2 = pl[t].A2, B2 = pl[t].B2;
            radd(ca[1], A2, B2); radd(ca[19], A2, B2);
            if (!quick) { low3(ca[9], ca[10], NULL, A2, B2, A2); low3(ca[11], ca[12], ca[13], A2, B2, A2); }
            { int32_t r3[3] = { pl[t].c, A2, B2 }; fwrite(r3, 4, 3, M->f_tpl); }
        }
        int ps = pl[t].match ? pl[t].A2 : A;
        {
            int lim = pl[t].match ? pl[t].B2 : B;
            for (int j = 0; j <= n; j++) {
                if (j < n && (mm[j] < ps || mm[j] >= lim)) continue;
                if (j < n && tol && pl[t].match && !pl[t].exact) {
                    int iso = (j == 0 || mm[j] - mm[j - 1] > wiso) && (j == n - 1 || mm[j + 1] - mm[j] > wiso);
                    if (iso) continue;
                }
                int pe = j < n ? mm[j] : lim;
                for (int w = 0; w < nw_lim; w++) if (pe - ps >= WS[w]) { ca[2 + w][ps + WS[w] - 1]++; ca[2 + w][pe]--; }
                ps = pe + 1;
            }
        }
    }
    if (tol && maxi > 0) for (int t = 0; t < np; t++) if (pl[t].ok && !pl[t].match) for (int u = 0; u < np; u++) {
        if (u == t || !pl[u].ok || pl[u].match || pl[t].match || pl[u].c != pl[t].c || pl[u].st != pl[t].st) continue;
        int dd = pl[u].d - pl[t].d; if (dd == 0 || dd > maxi || dd < -maxi) continue;
        Cand *C = &cd[pl[t].c]; int32_t **ca = M->arr[pl[t].c];
        const uint8_t *s = pl[t].st ? t_r : t_fs, *qq = pl[t].st ? t_rq : t_q;
        int dp = pl[t].d, ds = pl[u].d, ins = dd < 0 ? -dd : 0;
        int i0 = pl[t].i0, i1 = pl[u].i1; if (i1 - i0 < 60) continue;
        static __thread int pm[1026], sm[1026];
        pm[i0] = 0;
        for (int i = i0; i < i1; i++) { int x = dp + i; int df = (x >= 0 && x < C->len && s[i] < 4 && C->s[x] < 4 && s[i] != C->s[x] && qq[i] >= 33 + QMIN); pm[i + 1] = pm[i] + df; }
        sm[i1] = 0;
        for (int i = i1 - 1; i >= i0; i--) { int x = ds + i; int df = (x >= 0 && x < C->len && s[i] < 4 && C->s[x] < 4 && s[i] != C->s[x] && qq[i] >= 33 + QMIN); sm[i] = sm[i + 1] + df; }
        int best = 1 << 30, bb = -1;
        for (int b = i0 + 20; b + ins + 20 <= i1; b++) { if (dp + b > C->len || ds + b + ins < 0) continue; int cst = pm[b] + sm[b + ins]; if (cst < best) { best = cst; bb = b; } }
        if (bb < 0) continue;
        int tot = best + (dd < 0 ? -dd : dd);
        int ca_ = i0, cb_ = i1;
        for (int i = i0; i < i0 + 12 && i < bb; i++) if (pm[i + 1] - pm[i]) ca_ = i + 1;
        for (int i = i1 - 1; i >= i1 - 12 && i >= bb + ins; i--) if (sm[i] - sm[i + 1]) cb_ = i;
        tot -= pm[ca_] + sm[cb_];
        int A2 = dp + ca_, B2 = ds + cb_;
        if (A2 < C->g) continue;
        if (B2 - A2 < 60 || tot > 1 + (B2 - A2) / tol) continue;
        pl[t].match = 2; pl[t].A2 = A2; pl[t].B2 = B2; pl[t].indel = pl[u].indel = 1; pl[u].nqq = -1; M->n_stitch++;
        radd(ca[1], A2, B2); radd(ca[19], A2, B2);
        if (!quick) { low3(ca[9], ca[10], NULL, A2, B2, A2); low3(ca[11], ca[12], ca[13], A2, B2, A2); }
        { int32_t r3[3] = { pl[t].c, A2, B2 }; fwrite(r3, 4, 3, M->f_tpl); }
        break;
    }
    for (int t = 0; t < np; t++) if (pl[t].ok && !pl[t].exact && !pl[t].match) for (int u = 0; u < np; u++) {
        if (u == t || !pl[u].ok || pl[u].exact || pl[u].match || pl[u].c != pl[t].c || pl[u].st != pl[t].st) continue;
        int dd = pl[u].d - pl[t].d; if (dd > INDEL || dd < -INDEL || dd == 0) continue;
        int gap = pl[u].ls - pl[t].le - (dd > 0 ? dd : 0);
        if (gap >= -12 && gap <= 12) { pl[t].indel = pl[u].indel = 1; }
    }
    for (int t = 0; t < np; t++) {
        if (!pl[t].ok || pl[t].match) continue;
        int32_t **ca = M->arr[pl[t].c];
        if (pl[t].indel) { if (pl[t].nqq >= 0) radd(ca[19], pl[t].ls, pl[t].le); continue; }
        if (pl[t].nqq >= 1) radd(ca[19], pl[t].ps1, pl[t].pe1);
        if (pl[t].nqq >= 2 && pl[t].ps2 > pl[t].pe1) radd(ca[19], pl[t].ps2, pl[t].pe2);
    }
    for (int t = 0; t < np; t++) {
        if (!pl[t].ok || pl[t].exact || pl[t].indel) continue;
        int32_t **ca = M->arr[pl[t].c];
        pl[t].lv = pl[t].lvp >= 0;
        if (pl[t].bub) M->n_bub++;
        if (pl[t].lvp >= 0 && pl[t].enp >= 0) M->n_isl++; else if (pl[t].lvp >= 0) M->n_lvo++; else if (pl[t].enp >= 0) M->n_eno++;
        pl[t].en = pl[t].enp >= 0;
        if (pl[t].lv) ca[7][pl[t].lvp]++;
        if (pl[t].rgr >= 0 && pl[t].lvp < 0) ca[17][pl[t].rgr]++;
        if (pl[t].rgl >= 0 && pl[t].enp < 0) ca[18][pl[t].rgl]++;
        if (pl[t].en) ca[6][pl[t].enp]++;
    }
    M->n_placed += any;
}

static void events(VMode *M, Pl *a, int na, Pl *b, int nb, int have_mate) {
    const int bub = M->bub;
    for (int t = 0; t < na; t++) {
        int sf = a[t].sfp >= 0 && a[t].ok && !a[t].exact && !a[t].indel && a[t].lvp < 0;
        if (!a[t].en && !a[t].lv && !sf && !(bub && a[t].npat && a[t].ok)) continue;
        int up = -1, down = -1;
        int mA = have_mate ? -1 : -3, mB = mA;
        if (!have_mate && (a[t].en || a[t].lv)) M->n_orph++;
        if (have_mate) {
            if (a[t].st == 1) up = 0; else down = 0;
            for (int u = 0; u < nb; u++) {
                if (b[u].c != a[t].c || b[u].st == a[t].st || !b[u].ok) continue;
                int near_up = a[t].st == 1 && b[u].A <= a[t].A && a[t].B - b[u].A <= 1500;
                int near_dn = a[t].st == 0 && b[u].B >= a[t].B && b[u].B - a[t].A <= 1500;
                if (!near_up && !near_dn) continue;
                if (!b[u].match) { int cd_ = b[u].onp ? -2 : -4; if (mA == -1 || (mA == -4 && cd_ == -2)) mA = mB = cd_; continue; }
                if (near_up) up = 1;
                if (near_dn) down = 1;
                if (mA < 0 || (near_up ? b[u].A2 > mA : b[u].A2 < mA)) { mA = b[u].A2; mB = b[u].B2; }
            }
        }
        if (sf) { int32_t **ca = M->arr[a[t].c]; ca[14][a[t].sfp]++; if (have_mate && a[t].st == 1) { ca[15][a[t].sfp]++; if (mA >= 0) ca[16][a[t].sfp]++; } }
        if (a[t].en) { int32_t r8[8] = { a[t].c, a[t].enp, a[t].ee2, a[t].B, a[t].a402 + (bub ? 1000 * a[t].ecls : 0) + 10000 * a[t].gl, up, mA, mB }; fwrite(r8, 4, 8, M->f_ent); }
        if (bub && a[t].ok) for (int k = 0; k < a[t].npat; k++) {
            if (a[t].plo[k] - a[t].A >= MINON) { int32_t r8[8] = { a[t].c, a[t].A, a[t].A, a[t].plo[k], 1000, down, mA, mB }; fwrite(r8, 4, 8, M->f_lev); }
            if (a[t].B - a[t].phi[k] - 1 >= MINON) { int32_t r8[8] = { a[t].c, a[t].phi[k] + 1, a[t].B, a[t].B, 1000, up, mA, mB }; fwrite(r8, 4, 8, M->f_ent); }
        }
        if (a[t].lv) { int32_t r8[8] = { a[t].c, a[t].A, a[t].ls2, a[t].lvp, a[t].cr2 + (bub ? 1000 * a[t].lcls : 0), down, mA, mB }; fwrite(r8, 4, 8, M->f_lev); }
    }
}

static void finish_pair(VMode *M, Pl pl[2][MAXPL], int np[2], int have[2]) {
    if (have[0] && np[0]) events(M, pl[0], np[0], pl[1], np[1], have[1]);
    if (have[1] && np[1]) events(M, pl[1], np[1], pl[0], np[0], have[0]);
    if (have[0] && have[1] && np[0] && np[1]) {
        for (int c = 0; c < nc; c++) {
            int best = -1, bA = 0, bB = 0;
            for (int i = 0; i < np[0]; i++) if (pl[0][i].match && pl[0][i].c == c) for (int j = 0; j < np[1]; j++) if (pl[1][j].match && pl[1][j].c == c && pl[0][i].st != pl[1][j].st) {
                Pl *fw = pl[0][i].st == 0 ? &pl[0][i] : &pl[1][j], *rv = pl[0][i].st == 0 ? &pl[1][j] : &pl[0][i];
                int span = rv->B2 - fw->A2;
                if (fw->A2 <= rv->A2 && rv->B2 >= fw->B2 && span <= 1500 && (best < 0 || span < best)) { best = span; bA = fw->A2; bB = rv->B2; }
            }
            if (best >= 0) {
                int32_t r3[3] = { c, bA, bB }; fwrite(r3, 4, 3, M->f_tpl);
                int32_t **ca = M->arr[c];
                if (!M->quick) { radd(ca[8], bA, bB); low3(ca[11], ca[12], ca[13], bA, bB, bA); }
                if (c == 0) { M->n_pairs++; M->isz[best > 2000 ? 2000 : best]++; }
            }
        }
    }
}

static int mode_init(VMode *M, const char *outp) {
    M->arr = calloc((size_t)(nc ? nc : 1), sizeof(int32_t **));
    for (int c = 0; c < nc; c++) {
        M->arr[c] = calloc(NARR, sizeof(int32_t *));
        for (int j = 0; j < NARR; j++) {
            if (!arr_used(M->quick, j)) continue;
            M->arr[c][j] = calloc((size_t)cd[c].len + 2, 4);
            if (j >= 9 && j <= 13) for (int x = 0; x <= cd[c].len; x++) M->arr[c][j][x] = BIG;
        }
    }
    char fn[4096];
    snprintf(fn, sizeof fn, "%s.ent", outp); M->f_ent = fopen(fn, "wb");
    snprintf(fn, sizeof fn, "%s.lev", outp); M->f_lev = fopen(fn, "wb");
    snprintf(fn, sizeof fn, "%s.tpl", outp); M->f_tpl = fopen(fn, "wb");
    if (!M->f_ent || !M->f_lev || !M->f_tpl) { fprintf(stderr, "awv: cannot write %s.*: ", outp); perror(NULL); return 1; }
    setvbuf(M->f_ent, NULL, _IOFBF, 1 << 18);
    setvbuf(M->f_lev, NULL, _IOFBF, 1 << 18);
    setvbuf(M->f_tpl, NULL, _IOFBF, 1 << 18);
    return 0;
}

static int mode_finish(VMode *M, const char *outp, int jobno) {
    char fn[4096];
    snprintf(fn, sizeof fn, "%s.idx", outp); FILE *fi = fopen(fn, "w");
    snprintf(fn, sizeof fn, "%s.i32", outp); FILE *fo = fopen(fn, "wb");
    if (!fi || !fo) { fprintf(stderr, "awv: cannot write %s.*: ", outp); perror(NULL); return 1; }
    setvbuf(fo, NULL, _IOFBF, 1 << 20);
    int maxlen = 0;
    for (int c = 0; c < nc; c++) if (cd[c].len > maxlen) maxlen = cd[c].len;
    int32_t *zbuf = M->quick ? calloc((size_t)maxlen + 2, 4) : NULL;
    long off = 0;
    for (int c = 0; c < nc; c++) {
        Cand *C = &cd[c]; int32_t **ca = M->arr[c];
        for (int j = 0; j < NARR; j++) {
            if (ca[j] && ((j < 9 && j != 6 && j != 7) || j == 19)) {
                int32_t s = 0;
                for (int x = 0; x <= C->len; x++) { s += ca[j][x]; ca[j][x] = s; }
            }
        }
        fprintf(fi, "%s\t%d\t%ld\n", C->name, C->len, off);
        for (int j = 0; j < NARR; j++) {
            fwrite(ca[j] ? ca[j] : zbuf, 4, (size_t)C->len + 1, fo);
            off += 4L * (C->len + 1);
        }
    }
    free(zbuf);
    int bad = ferror(fi) | ferror(fo) | ferror(M->f_ent) | ferror(M->f_lev) | ferror(M->f_tpl);
    if (fclose(fi) | fclose(fo) | fclose(M->f_ent) | fclose(M->f_lev) | fclose(M->f_tpl) | bad) {
        fprintf(stderr, "awv: error writing %s.* (disk full?)\n", outp);
        return 1;
    }
    long tot = 0, acc = 0; int med = 0;
    for (int i = 0; i <= 2000; i++) tot += M->isz[i];
    for (int i = 0; i <= 2000; i++) { acc += M->isz[i]; if (acc * 2 >= tot) { med = i; break; } }
    fprintf(stderr, "awv: job %d: %d candidates, %lu reads laid, %lu matching placements, %lu matching pairs on candidate 0 (median outer distance %d); "
                    "placements that leave only %lu, come on only %lu, both %lu; events of reads without a mate %lu; tolerant matches %lu, stitched %lu\n",
            jobno, nc, (unsigned long)M->n_placed, (unsigned long)M->n_exact, (unsigned long)M->n_pairs, med,
            (unsigned long)M->n_lvo, (unsigned long)M->n_eno, (unsigned long)M->n_isl, (unsigned long)M->n_orph,
            (unsigned long)M->n_tol, (unsigned long)M->n_stitch);
    for (int c = 0; c < nc; c++) {
        for (int j = 0; j < NARR; j++) free(M->arr[c][j]);
        free(M->arr[c]);
    }
    free(M->arr);
    return 0;
}

/* ---------- the reads, held in memory once; the genes, each by one thread ---------- */
typedef struct { uint64_t off; long pair; uint16_t len; uint8_t mate; } RRec;
static uint8_t *RS; static RRec *RR; static size_t NRR = 0;
typedef struct { char *cfa, *outp, *outp2; } Job;
static Job *jobs; static int njobs = 0, next_job = 0, failed = 0, dual_mode = 0;
static pthread_mutex_t jm = PTHREAD_MUTEX_INITIALIZER;

static int one_gene(const Job *J, int jobno) {
    nc = 0; cd = NULL; npost = 0;
    char *line = NULL; size_t cap = 0; ssize_t n; int capc = 0;
    FILE *f = fopen(J->cfa, "r"); if (!f) { fprintf(stderr, "awv: cannot open %s: ", J->cfa); perror(NULL); return 1; }
    while ((n = getline(&line, &cap, f)) > 0) {
        while (n && (line[n - 1] == '\n' || line[n - 1] == '\r')) line[--n] = 0;
        if (line[0] == '>') {
            if (nc == capc) { capc = capc ? capc * 2 : 64; cd = realloc(cd, (size_t)capc * sizeof(Cand)); }
            cd[nc].name = strdup(line + 1); cd[nc].s = NULL; cd[nc].len = 0; cd[nc].g = 0;
            char *p = line; for (int k = 0; k < 3 && p; k++) { p = strchr(p, '|'); if (p) p++; }
            if (p) cd[nc].g = atoi(p);
            nc++;
        } else if (nc) {
            Cand *C = &cd[nc - 1]; C->s = realloc(C->s, (size_t)C->len + (size_t)n + 1);
            for (ssize_t i = 0; i < n; i++) C->s[C->len++] = (uint8_t)code[(unsigned char)line[i]];
        }
    }
    fclose(f); free(line);
    for (int c = 0; c < nc; c++) npost += cd[c].len >= K ? (size_t)(cd[c].len - K + 1) : 0;
    post = malloc((npost + 1) * sizeof(Post)); size_t m = 0;
    const uint64_t kmask = (1ULL << (2 * K)) - 1;
    for (int c = 0; c < nc; c++) {
        Cand *C = &cd[c];
        uint64_t km = 0; int valid = 0;
        for (int p = 0; p < C->len; p++) {
            uint8_t b = C->s[p];
            if (b > 3) { valid = 0; km = 0; } else { km = ((km << 2) | b) & kmask; valid++; }
            if (valid >= K) post[m++] = (Post){ km, p - K + 1, (uint16_t)c };
        }
    }
    npost = m; qsort(post, npost, sizeof(Post), cmp_post);
    uint64_t hs = 1024; while (hs < npost * 2 + 2) hs <<= 1; hmask = hs - 1; ht = calloc(hs, 4);
    for (size_t j = 0; j < npost; j++) {
        if (j && post[j].km == post[j - 1].km) continue;
        uint64_t h = mix(post[j].km) & hmask; while (ht[h]) h = (h + 1) & hmask; ht[h] = (uint32_t)j + 1;
    }

    VMode M1, M2;
    memset(&M1, 0, sizeof M1); memset(&M2, 0, sizeof M2);
    if (dual_mode) {
        /* M1 = strict (0.3.0 defaults), M2 = tolerant (options given on command line) */
        M1.kd = 3; M1.wiso = 12; M1.patchmax = 40; M1.patchdiff = 12; M1.totper = 12; M1.quick = QUICK;
        M2.tol = TOL; M2.kd = KD; M2.maxi = MAXI; M2.bub = BUB; M2.wiso = WISO;
        M2.patchmax = PATCHMAX; M2.patchdiff = PATCHDIFF; M2.totper = TOTPER; M2.w100 = W100; M2.quick = QUICK;
        if (mode_init(&M1, J->outp) || mode_init(&M2, J->outp2)) return 1;
    } else {
        M1.tol = TOL; M1.kd = KD; M1.maxi = MAXI; M1.bub = BUB; M1.wiso = WISO;
        M1.patchmax = PATCHMAX; M1.patchdiff = PATCHDIFF; M1.totper = TOTPER; M1.w100 = W100; M1.quick = QUICK;
        if (mode_init(&M1, J->outp)) return 1;
    }

    static __thread BaseHit bh[MAXPL];
    static __thread int mm_pool[MAXPL * 1024];
    static __thread Pl pl1[2][MAXPL], pl2[2][MAXPL];
    int np1[2] = {0, 0}, np2[2] = {0, 0}, have[2] = {0, 0};
    long cur = -1;
    for (size_t ri = 0; ; ri++) {
        int more = ri < NRR ? 1 : 0;
        long pair = -2; int mate = 0, L = 0;
        if (more) { pair = RR[ri].pair; mate = RR[ri].mate; L = RR[ri].len; }
        if (pair != cur) {
            finish_pair(&M1, pl1, np1, have);
            if (dual_mode) finish_pair(&M2, pl2, np2, have);
            have[0] = have[1] = 0; np1[0] = np1[1] = 0; np2[0] = np2[1] = 0; cur = pair;
        }
        if (!more) break;
        int any = 0;
        int np = seed_read(RS + RR[ri].off, L, bh, mm_pool, &any);
        np1[mate] = np;
        if (np) eval_mode(&M1, L, bh, np, mm_pool, any, pl1[mate]);
        if (dual_mode) {
            np2[mate] = np;
            if (np) eval_mode(&M2, L, bh, np, mm_pool, any, pl2[mate]);
        }
        have[mate] = 1;
    }

    int rc = mode_finish(&M1, J->outp, jobno);
    if (dual_mode && mode_finish(&M2, J->outp2, jobno)) rc = 1;
    for (int c = 0; c < nc; c++) { free(cd[c].name); free(cd[c].s); }
    free(cd); free(post); free(ht); cd = NULL; post = NULL; ht = NULL;
    return rc;
}

static void *worker(void *arg) {
    (void)arg;
    for (;;) {
        pthread_mutex_lock(&jm); int j = next_job < njobs ? next_job++ : -1; pthread_mutex_unlock(&jm);
        if (j < 0) break;
        int rc = one_gene(&jobs[j], j);
        pthread_mutex_lock(&jm);
        if (rc) failed = 1; else { printf("%d\n", j); fflush(stdout); }
        pthread_mutex_unlock(&jm);
    }
    return NULL;
}

int main(int argc, char **argv) {
    int nt = 1, c; const char *dual_fn = NULL;
    while ((c = getopt(argc, argv, "t:vT:K:I:g:B:m:n:o:w:Z:q")) != -1) {
        if (c == 't') nt = atoi(optarg);
        else if (c == 'T') TOL = atoi(optarg);
        else if (c == 'K') KD = atoi(optarg);
        else if (c == 'I') MAXI = atoi(optarg);
        else if (c == 'B') BUB = atoi(optarg);
        else if (c == 'm') PATCHMAX = atoi(optarg);
        else if (c == 'n') PATCHDIFF = atoi(optarg);
        else if (c == 'o') TOTPER = atoi(optarg);
        else if (c == 'w') W100 = atoi(optarg);
        else if (c == 'g') WISO = atoi(optarg);
        else if (c == 'Z') dual_fn = optarg;
        else if (c == 'q') QUICK = 1;
        else if (c == 'v') { printf("awv %s\n", AWV_VERSION); return 0; }
        else { argc = 0; break; }
    }
    if (argc - optind != 2) {
        fprintf(stderr, "awv %s\nusage: awv [-t 1] [-q] [-Z vt_jobs.tsv] [-T 50 -B 1 -I 12 -g 30 -w 9] jobs.tsv reads.tsv\n", AWV_VERSION);
        return 1;
    }
    if (nt < 1) nt = 1; if (nt > 64) nt = 64;
    memset(code, 4, sizeof code); code['A'] = code['a'] = 0; code['C'] = code['c'] = 1; code['G'] = code['g'] = 2; code['T'] = code['t'] = 3;
    char *line = NULL; size_t cap = 0; ssize_t n;
    FILE *fj = fopen(argv[optind], "r"); if (!fj) { fprintf(stderr, "awv: cannot open %s: ", argv[optind]); perror(NULL); return 1; }
    int capj = 0;
    while ((n = getline(&line, &cap, fj)) > 0) {
        while (n && (line[n - 1] == '\n' || line[n - 1] == '\r')) line[--n] = 0;
        char *t = strchr(line, '\t'); if (!t) continue; *t = 0;
        if (njobs == capj) { capj = capj ? capj * 2 : 64; jobs = realloc(jobs, (size_t)capj * sizeof(Job)); }
        jobs[njobs].cfa = strdup(line); jobs[njobs].outp = strdup(t + 1); jobs[njobs].outp2 = NULL; njobs++;
    }
    fclose(fj);
    if (dual_fn) {
        dual_mode = 1;
        FILE *fj2 = fopen(dual_fn, "r"); if (!fj2) { fprintf(stderr, "awv: cannot open %s: ", dual_fn); perror(NULL); return 1; }
        int k = 0;
        while ((n = getline(&line, &cap, fj2)) > 0 && k < njobs) {
            while (n && (line[n - 1] == '\n' || line[n - 1] == '\r')) line[--n] = 0;
            char *t = strchr(line, '\t'); if (!t) continue;
            jobs[k++].outp2 = strdup(t + 1);
        }
        fclose(fj2);
        if (k != njobs) { fprintf(stderr, "awv: %s and %s have different numbers of jobs\n", argv[optind], dual_fn); return 1; }
    }
    const char *reads = argv[optind + 1];
    FILE *fr = strcmp(reads, "-") ? fopen(reads, "r") : stdin; if (!fr) { fprintf(stderr, "awv: cannot open %s: ", reads); perror(NULL); return 1; }
    size_t caps = 1 << 24, ns = 0, capr = 1 << 16; RS = malloc(caps); RR = malloc(capr * sizeof(RRec));
    while ((n = getline(&line, &cap, fr)) > 0) {
        char *p = line; long pair = strtol(p, &p, 10); int mate = (int)strtol(p + 1, &p, 10) - 1; char *seq = p + 1;
        char *t1 = strchr(seq, '\t'); if (!t1) continue; int L = (int)(t1 - seq);
        char *t2 = strchr(t1 + 1, '\t'); char *qual = t2 ? t2 + 1 : NULL;
        if (mate < 0 || mate > 1) continue;
        if (L > 1000) L = 1000;
        if (ns + (size_t)L + 1 > caps) { while (ns + (size_t)L + 1 > caps) caps *= 2; RS = realloc(RS, caps); }
        if (NRR == capr) { capr *= 2; RR = realloc(RR, capr * sizeof(RRec)); }
        if (!RS || !RR) { fprintf(stderr, "awv: out of memory while loading the reads\n"); return 2; }
        int ql = 0; if (qual) while (qual[ql] && qual[ql] != '\n' && qual[ql] != '\t') ql++;
        for (int i = 0; i < L; i++) RS[ns + (size_t)i] = (uint8_t)code[(unsigned char)seq[i]] | ((qual && ql >= L && (uint8_t)qual[i] < 33 + QMIN) ? 8 : 0);
        RR[NRR++] = (RRec){ ns, pair, (uint16_t)L, (uint8_t)mate }; ns += (size_t)L;
    }
    if (fr != stdin) fclose(fr);
    free(line);
    fprintf(stderr, "awv: %zu reads, %d genes, %d threads%s\n", NRR, njobs, nt, dual_mode ? " (strict+tolerant)" : "");
    if (nt > njobs) nt = njobs;
    pthread_t *th = malloc(sizeof(pthread_t) * (size_t)(nt > 0 ? nt : 1));
    for (int i = 0; i < nt; i++) if (pthread_create(&th[i], NULL, worker, NULL)) { perror("awv: cannot start a thread"); return 1; }
    for (int i = 0; i < nt; i++) pthread_join(th[i], NULL);
    return failed ? 1 : 0;
}
