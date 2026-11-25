// awx.c (ARGjacent 0.5.0) - the proposer: from every gene end a path is extended base by base with the reads kept near
// the genes, and where the reads part, up to -c arms per gene side are followed, the best-supported first.
// Since 0.4 (19 Nov 2026) two things differ from 0.3.0, each an option that is on by default:
//   -P 1  an arm that carries another version of a few bases and agrees with the chosen arm again 31 to 60 bases on is a variant of
//         it, not a candidate of its own (0.3.0 knew single bases and small insertions or deletions as variants; a gene side's sixteen
//         candidates could be versions of ONE context that differ in such patches);
//   -Q 4  a read that disagrees with the path at a place classed as a variant STAYS on the path, up to 4 such places per read. 0.3.0
//         dropped a read at its first disagreement, so after a base in which two versions differ every read of the other version was
//         gone for the rest of its length, and a fork that only that version takes was seen by few of its reads or by none. The base
//         of the path is still chosen by the reads without a difference; arms are looked for among all reads on the path; the path
//         remembers its variant places, and a read found later that differs only there is laid on it too.
//   -P 0 -Q 0 gives the candidates of 0.3.0, byte for byte.
// On three gut samples with long reads a candidate then matches 83 / 74 / 54 % of the gene-carrying molecules with 1 / 2 / 5 kb of
// flank (0.3.0's proposer: 75 / 66 / 46 %); the rules decide as before what is certain.
// awx9.c - awx7, and the last REDO bases of the gene are taken from the sample's reads instead of the database sequence
// (argv[6], default 0): reads are laid on the whole gene part as before, then the path starts REDO bases before the gene end.
// A sample's own version of a gene often differs from the database allele; candidates that begin with the database
// sequence then have no read matching at the gene, and the verifier vouches for nothing.
// awx7.c - awx6 without the limit of two other arms per 50 bp (ALTMAX).
// awx6.c - awx5, and forks are looked for among all reads that match the path end (overlap of 31 bp or more) at every
// place, not only after another locus came on: an arm whose reads differ from the path further back was never seen.
// awx5.c - awx3 with another order of exploration (16 Nov 2026): the arm that carries the largest share of the gene's reads is
// extended next (share = product of the arms' shares of reads at the forks passed, times the flank still to go), the pool of
// waiting arms keeps the best when it is full instead of dropping late arms, and the number of candidates is an option.
// awx.c - prototype base-level proposer. From the end of the gene the path is extended one base at a time by the reads that
// match its end without difference. Where the reads disagree, those that share the most of the path decide (read-backed
// phasing); where they still disagree the path forks: the arm linked by read pairs to the part before another locus came on
// is followed first, the other arms become further candidates (structural arms only, not single-base or small indel variants).
//   usage: awx seeds.fa out.fa [max_cands=16] [flank=10000] < reads.tsv
//   seeds.fa : >gene|side|GLEN  + the gene part running outward (last GLEN bp of the gene, or reverse complement of the first)
//   reads.tsv: pair \t mate \t sequence \t anchors \t qualities ; mates on adjacent lines
//   out.fa   : >gene|side|n|GLEN|0|0  candidates (gene part + flank); a log of every candidate goes to stderr
// Build: cc -O3 -o awx awx.c
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <time.h>
#include <pthread.h>
#define AWX_VERSION "0.5.1"
#define K 31
#define SAMP 8
#define QMIN 15
#define REDO 0              /* the gene part is the seed as given; reads are laid on it with differences tolerated */
#define MAXQ 512
#define LA 60
static int PATCH = 1;
#ifndef ALTMAX
#define ALTMAX 64
#endif
typedef struct { uint64_t off; uint16_t len; } Read;
typedef struct __attribute__((packed)) { uint64_t km; uint32_t r; uint16_t o; } Ent;      /* 14 bytes: there are about 18 per read */
typedef struct { uint32_t r; int32_t d; uint8_t s; uint8_t nt; int32_t skip; } Act;   /* awxq: nt = variant places passed with another version; skip = no say before this path position */
typedef struct { uint32_t k; int32_t d; } PE;         /* a read placed on a path: k = ((read + 1) << 1) | strand, 0 = free; d = its diagonal */
typedef struct {
    uint8_t *P; uint8_t *V; int len, capP;        /* awxq: V[p] = 1: a place where versions differ in this one base; 2: inside a small insertion, deletion or patch */
    Act *act; int na, capa;
    PE *pm; uint32_t pmask, pn;        /* the reads that were laid on this path, with diagonal and strand of their matching placement */
    int jfirst, jpos, jcnt;            /* earliest place where another locus came on; current cluster */
    int pen, nvar, nfork, nthin, rep, born, altpos, altn;
    double mass, pri;                  /* share of the gene's reads that came this way; priority while waiting */
    char why[48];
} St;
static uint8_t *rs; static Read *rd; static uint32_t nr = 0;
static Ent *ix; static size_t nix = 0;
#define PB 24                                     /* the index is entered through a table of the K-mers' first 12 bases */
static uint32_t *pt;
static int8_t code[256];
/* what belongs to the gene end being extended: one per thread */
static __thread int EX0; static __thread int redo = 0;
static __thread uint64_t work = 0, n_verify = 0, n_hot = 0, n_clone = 0;
#define MAXHIT 20000
#define WORKMAX 12000000000ULL                       /* positions before EX0 are the seed gene: differences tolerated */
static inline uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static int cmp_ent(const void *a, const void *b) { const Ent *x = a, *y = b; if (x->km != y->km) return x->km < y->km ? -1 : 1; if (x->r != y->r) return x->r < y->r ? -1 : 1; return (int)x->o - (int)y->o; }
static long find(uint64_t km) {                     /* the first index entry with this K-mer, or -1 */
    size_t lo = pt[km >> (2 * K - PB)], hi = pt[(km >> (2 * K - PB)) + 1];
    while (lo < hi) { size_t mid = (lo + hi) / 2; if (ix[mid].km < km) lo = mid + 1; else hi = mid; }
    return lo < pt[(km >> (2 * K - PB)) + 1] && ix[lo].km == km ? (long)lo : -1;
}
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static inline uint8_t rb(uint32_t r, int s, int i) {      /* base i of read r in orientation s; bit 3 = low quality */
    const Read *R = &rd[r]; uint8_t c = s ? rs[R->off + R->len - 1 - i] : rs[R->off + i]; uint8_t b = c & 7;
    if (s && b < 4) b = 3 - b;
    return b | (c & 8);
}
static St *st_new(int cap) {
    St *s = calloc(1, sizeof(St)); s->capP = cap; s->P = malloc(cap); s->V = calloc((size_t)cap, 1); s->capa = 256; s->act = malloc(s->capa * sizeof(Act));
    s->pmask = 1023; s->pm = calloc(s->pmask + 1, sizeof(PE)); s->jfirst = -1; s->jpos = -1000; s->altpos = -1000; return s;
}
static St *st_clone(const St *a) {
    St *s = malloc(sizeof(St)); *s = *a; s->P = malloc(a->capP); memcpy(s->P, a->P, a->len); s->V = malloc(a->capP); memcpy(s->V, a->V, a->capP);
    s->act = malloc(a->capa * sizeof(Act)); memcpy(s->act, a->act, a->na * sizeof(Act));
    s->pm = malloc(((size_t)a->pmask + 1) * sizeof(PE)); memcpy(s->pm, a->pm, ((size_t)a->pmask + 1) * sizeof(PE)); return s;
}
static void st_free(St *s) { free(s->P); free(s->V); free(s->act); free(s->pm); free(s); }
static void st_slim(St *s) { free(s->act); free(s->pm); free(s->V); s->act = NULL; s->pm = NULL; s->V = NULL; }   /* a finished state needs only its path */
static inline const PE *pm_find(const St *s, uint32_t r) {       /* the placement of read r on this path, or NULL */
    for (uint32_t h = (r * 2654435761u) & s->pmask; ; h = (h + 1) & s->pmask) {
        const PE *e = &s->pm[h];
        if (!e->k) return NULL;
        if ((e->k >> 1) == r + 1) return e;
    }
}
static void pm_set(St *s, uint32_t r, int sd, int32_t d) {
    uint32_t h = (r * 2654435761u) & s->pmask;
    while (s->pm[h].k && (s->pm[h].k >> 1) != r + 1) h = (h + 1) & s->pmask;
    if (!s->pm[h].k) s->pn++;
    s->pm[h] = (PE){ ((r + 1) << 1) | (uint32_t)(sd & 1), d };
    if (s->pn * 2 > s->pmask) {
        PE *o = s->pm; uint32_t om = s->pmask;
        s->pmask = om * 2 + 1; s->pm = calloc((size_t)s->pmask + 1, sizeof(PE));
        for (uint32_t i = 0; i <= om; i++) if (o[i].k) { uint32_t t = (((o[i].k >> 1) - 1) * 2654435761u) & s->pmask; while (s->pm[t].k) t = (t + 1) & s->pmask; s->pm[t] = o[i]; }
        free(o);
    }
}
static int pool_cap = 64;
static void act_add(St *s, uint32_t r, int sd, int d, int nt) {
    if (s->na == s->capa) { s->capa *= 2; s->act = realloc(s->act, s->capa * sizeof(Act)); }
    s->act[s->na++] = (Act){ r, d, (uint8_t)sd, (uint8_t)nt, 0 }; pm_set(s, r, sd, d);
}
static int QTOL = 4;                                     /* awxq: variant places a read may pass with another version and stay on the path */
static void vflag(St *s, int ty, int sh) {               /* the path remembers a variant place that begins at its next position */
    int a = s->len, e = a + (ty == 2 ? 31 : 1) + (sh < 0 ? -sh : sh);
    if (ty == 1 && sh == 0) { if (a < s->capP && !s->V[a]) s->V[a] = 1; return; }
    for (int p = a; p <= e && p < s->capP; p++) s->V[p] = 2;
}
/* does read r (orientation sd, diagonal d) match the path up to its last base? 1 yes; 0 no; 2 no, and it comes on from elsewhere at *where */
static int verify(const St *s, uint32_t r, int sd, int d, int *where, int *ntol) {
    int L = rd[r].len, x = s->len - 1, i1 = x - d, i0 = d < 0 ? -d : 0;
    if (i1 >= L) i1 = L - 1;
    if (i1 - i0 + 1 < K) return 0;
    int mmg = 0, ng = 0, first = -1, cl = 0, nt = 0;
    for (int i = i1; i >= i0; i--) {
        uint8_t b = rb(r, sd, i); int p = d + i;
        if (p < EX0) { ng++; if (!(b & 8) && (b & 7) < 4 && (b & 7) != s->P[p]) mmg++; continue; }
        if ((b & 8) || (b & 7) > 3) continue;
        if ((b & 7) != s->P[p]) {
            if (QTOL && first < 0 && nt < QTOL && s->V[p] == 1) { nt++; continue; }     /* awxq: another version of one base */
            if (QTOL && first < 0 && nt < QTOL && s->V[p] == 2) { nt++; break; }        /* ... of a few: what lies before is not compared */
            if (first < 0) first = p;
            if (p >= first - 30) cl++; else break;
        }
    }
    if (first < 0) { *ntol = nt; return mmg * 10 <= ng; }
    int e = x - first;                                   /* stretch without difference up to the end of the path */
    if (e >= K && e < K + SAMP && cl >= 3 && first - (d + i0) >= 12) { *where = first + 1; return 2; }
    return 0;
}
static void discover(St *s) {                            /* reads that carry the last K bases of the path */
    if (s->len < K) return;
    uint64_t km = 0, rc = 0;
    for (int i = 0; i < K; i++) { uint8_t b = s->P[s->len - K + i]; if (b > 3) return; km = (km << 2) | b; rc = (rc >> 2) | ((uint64_t)(3 - b) << (2 * (K - 1))); }
    for (int o2 = 0; o2 < 2; o2++) {
        long p = find(o2 ? rc : km); if (p < 0) continue;
        { size_t e = (size_t)p + MAXHIT; if (e < nix && ix[e].km == (o2 ? rc : km)) { n_hot++; continue; } }   /* a K-mer in more than MAXHIT reads says nothing */
        for (size_t j = (size_t)p; j < nix && ix[j].km == (o2 ? rc : km); j++) {
            uint32_t r = ix[j].r; int L = rd[r].len, o = o2 ? L - K - ix[j].o : ix[j].o, d = (s->len - K) - o;
            { const PE *e = pm_find(s, r); if (e && (int)(e->k & 1) == o2 && e->d == d) continue; }
            int where = 0, ntol = 0, v = verify(s, r, o2, d, &where, &ntol); n_verify++; work += 60;
            if (v == 1) act_add(s, r, o2, d, ntol);
            else if (v == 2) {
                uint32_t m = r ^ 1;                       /* a variant that comes back (its mate matches the path before) is not another locus */
                if (m < nr) { const PE *e = pm_find(s, m); if (e && e->d < where - 60) continue; }
                if (where - s->jpos <= 3) s->jcnt++; else { s->jpos = where; s->jcnt = 1; }
                if (s->jcnt >= 2 && s->jfirst < 0 && where >= EX0) s->jfirst = s->jpos;
            }
        }
    }
}
static const int TIER[4] = { 110, 80, 55, K };
static int flank_target = 10000; static __thread int glen;
/* extend one state to its end; arms not taken are pushed onto q */
static void extend(St *s, St **q, int *nq, int maxq) {
    int target = glen + flank_target;
    static __thread uint64_t *seen = NULL; static __thread int seen_cap = 0;
    if (!seen) { seen_cap = 1 << 16; while (seen_cap < 4 * (target + 64)) seen_cap <<= 1; seen = malloc((size_t)seen_cap * 8); }   /* the K-mers of one path: never more than a quarter full */
    memset(seen, 0xff, (size_t)seen_cap * 8);
    for (int p = K - 1; p < s->len; p++) { uint64_t km = 0; int bad = 0; for (int i = 0; i < K; i++) { uint8_t b = s->P[p - K + 1 + i]; if (b > 3) { bad = 1; break; } km = (km << 2) | b; } if (bad) continue; uint64_t h = mix(km) & (seen_cap - 1); while (seen[h] != ~0ULL && seen[h] != km) h = (h + 1) & (seen_cap - 1); seen[h] = km; }
    strcpy(s->why, "reached target");
    while (s->len < target) {
        if (work > WORKMAX) { strcpy(s->why, "work cap"); break; }
        work += (uint64_t)s->na;
        int cnt[4][4]; memset(cnt, 0, sizeof cnt); int cntA[4] = {0};
        int w = 0;
        for (int a = 0; a < s->na; a++) {
            Act A = s->act[a]; int i = s->len - A.d;
            if (i >= rd[A.r].len) continue;              /* the read ends here */
            s->act[w++] = A;
            if (i < 0) continue;                         /* awx9: laid on the seed further on; not reached yet */
            if (A.skip > s->len) continue;               /* awxq: inside a variant it does not carry: no say */
            uint8_t b = rb(A.r, A.s, i); if ((b & 8) || (b & 7) > 3) continue;
            int ov = s->len - (A.d > 0 ? A.d : 0);
            cntA[b & 7]++;
            if (A.nt) continue;                          /* awxq: the base is chosen by the reads without a difference */
            for (int t = 0; t < 4; t++) if (ov >= TIER[t]) cnt[t][b & 7]++;
        }
        s->na = w;
        int tot3 = cnt[3][0] + cnt[3][1] + cnt[3][2] + cnt[3][3];
        int totA = cntA[0] + cntA[1] + cntA[2] + cntA[3];
        if (QTOL && tot3 < 2 && totA > tot3) { memset(cnt, 0, sizeof cnt); for (int b = 0; b < 4; b++) cnt[3][b] = cntA[b]; tot3 = totA; }   /* awxq: (nearly) none without a difference: all decide */
        if (tot3 < 1) { strcpy(s->why, "ran out of reads"); break; }
        int t = 3;
        for (int u = 0; u < 3; u++) if (cnt[u][0] + cnt[u][1] + cnt[u][2] + cnt[u][3] >= 4) { t = u; break; }
        int b1 = 0, b2 = -1;
        for (int b = 1; b < 4; b++) if (cnt[t][b] > cnt[t][b1]) b1 = b;
        for (int b = 0; b < 4; b++) if (b != b1 && (b2 < 0 || cnt[t][b] > cnt[t][b2])) b2 = b;
        int n1 = cnt[t][b1], n2 = cnt[t][b2], choose = b1;
        int split = n2 >= 2 && n2 * 100 >= 15 * (n1 + n2);
        int anyvar = 0; for (int b = 0; b < 4; b++) if (b != b1 && cnt[3][b] >= 2 && cnt[3][b] * 100 >= 15 * tot3) anyvar = 1;
        if (QTOL) for (int b = 0; b < 4; b++) if (b != b1 && cntA[b] >= 2 && cntA[b] * 100 >= 15 * totA) anyvar = 1;
        int isvar[4] = {0}, vsh[4] = {0}, vty[4] = {0};
        if (tot3 == 1) s->nthin++;
        if (split || anyvar) {            /* awx6: an arm seen only by reads with a short overlap is looked at everywhere, not only after another locus came on */
            /* arms: what the reads of each base say next, and how many are linked by their mates to the path before the join */
            int n[4] = {0}, anch[4] = {0}; uint8_t la[4][LA + 1]; int lc[4][LA + 1][4]; memset(lc, 0, sizeof lc);
            for (int a = 0; a < s->na; a++) {
                Act A = s->act[a]; int i = s->len - A.d; if (i < 0 || A.skip > s->len) continue; uint8_t b = rb(A.r, A.s, i); if ((b & 8) || (b & 7) > 3) continue; b &= 7;
                n[b]++;
                uint32_t m = A.r ^ 1;
                if (m < nr) { const PE *e = pm_find(s, m); if (e && (e->k & 1) != A.s && e->d < A.d && A.d - e->d <= 1200 && (s->jfirst < 0 || (e->d > 0 ? e->d : 0) + 40 <= s->jfirst)) anch[b]++; }
                for (int k = 1; k <= LA && i + k < rd[A.r].len; k++) { uint8_t c = rb(A.r, A.s, i + k); if (!(c & 8) && (c & 7) < 4) lc[b][k][c & 7]++; }
            }
            for (int b = 0; b < 4; b++) for (int k = 1; k <= LA; k++) { int m = 0; for (int c = 1; c < 4; c++) if (lc[b][k][c] > lc[b][k][m]) m = c; la[b][k] = lc[b][k][m] >= (k > 25 ? 2 : 1) ? (uint8_t)m : 9; }
            int tot = n[0] + n[1] + n[2] + n[3];
            if (s->jfirst >= 0) { int best = -1; for (int b = 0; b < 4; b++) if (anch[b] >= 2 && (best < 0 || anch[b] > anch[best])) best = b; if (best >= 0) choose = best; }
            int nstruct = 0; int isarm[4] = {0}; int armtot = n[choose];
            for (int b = 0; b < 4; b++) {
                if (b == choose || n[b] < 2) continue;
                if (!(n[b] * 100 >= 15 * tot || anch[b] >= 2)) continue;
                /* a single-base or small insertion/deletion variant of the chosen arm? then it is not another candidate */
                int variant = 0;
                for (int si = 0; si <= 24 && !variant; si++) {
                    int sh = (si & 1) ? (si + 1) / 2 : -(si / 2);
                    int same = 0, cmpn = 0;
                    for (int k = 1; k <= 25; k++) { int k2 = k + sh; if (k2 < 1 || k2 > 25 || la[b][k] > 3 || la[choose][k2] > 3) continue; cmpn++; same += la[b][k] == la[choose][k2]; }
                    if (cmpn >= 10 && same >= cmpn - 1) { variant = 1; vsh[b] = sh; vty[b] = 1; }
                }
                if (!variant && PATCH) for (int si = 0; si <= 24 && !variant; si++) {
                    /* the arms agree again: bases 31 to 60 after this place, but for one */
                    int sh = (si & 1) ? (si + 1) / 2 : -(si / 2);
                    int same = 0, cmpn = 0;
                    for (int k = 31; k <= LA; k++) { int k2 = k + sh; if (k2 < 26 || k2 > LA || la[b][k] > 3 || la[choose][k2] > 3) continue; cmpn++; same += la[b][k] == la[choose][k2]; }
                    if (cmpn >= 15 && same >= cmpn - 1) { variant = 1; vsh[b] = sh; vty[b] = 2; }
                }
                if (variant) { s->nvar++; if (QTOL) { isvar[b] = 1; vflag(s, vty[b], vsh[b]); } continue; }
                s->nfork++; nstruct++; isarm[b] = 1; armtot += n[b];
            }
            double mass0 = s->mass;
            for (int b = 0; b < 4; b++) {
                if (!isarm[b]) continue;
                if (s->len - s->altpos > 50) { s->altpos = s->len; s->altn = 0; }
                if (s->altn >= ALTMAX) continue;         /* awx7: was 2; where many contexts part at one place (the end of a shared element) every arm is wanted */
                double cm = mass0 * (double)n[b] / (double)armtot, cp = cm * (double)(target - s->len) / (double)flank_target;
                int slot = -1;
                if (*nq < pool_cap && *nq < maxq) slot = (*nq)++;
                else { int wi = 0; for (int i = 1; i < *nq; i++) if (q[i]->pri < q[wi]->pri) wi = i; if (q[wi]->pri < cp) { st_free(q[wi]); slot = wi; } }
                if (slot < 0) continue;
                s->altn++;
                St *c = st_clone(s); int w2 = 0; n_clone++; work += (uint64_t)nr / 8;
                for (int a = 0; a < c->na; a++) { Act A = c->act[a]; if (c->len - A.d < 0 || A.skip > c->len) { c->act[w2++] = A; continue; } uint8_t x = rb(A.r, A.s, c->len - A.d); if ((x & 8) || (x & 7) > 3 || (x & 7) == b) c->act[w2++] = A; }
                c->na = w2; c->P[c->len++] = (uint8_t)b; c->pen = s->pen + (anch[b] >= 2 ? 1 : 2) + (n[b] * 100 < 15 * tot); c->born = s->len - glen;
                c->mass = cm; c->pri = cp;
                discover(c); q[slot] = c;
            }
            if (nstruct) s->mass = mass0 * (double)n[choose] / (double)armtot;
        } else if (anyvar) s->nvar++;
        /* take the base: reads that say otherwise are off the path */
        w = 0;
        for (int a = 0; a < s->na; a++) {
            Act A = s->act[a]; if (s->len - A.d < 0 || A.skip > s->len) { s->act[w++] = A; continue; }
            uint8_t x = rb(A.r, A.s, s->len - A.d);
            if ((x & 8) || (x & 7) > 3 || (x & 7) == choose) s->act[w++] = A;
            else if (QTOL && isvar[x & 7] && A.nt < QTOL) {      /* awxq: another version of a few bases: the read stays on the path */
                int sh = vsh[x & 7]; A.nt++;
                if (vty[x & 7] == 2 || sh) { A.d += sh; A.skip = s->len + 1 + (vty[x & 7] == 2 ? 30 : 0) + (sh > 0 ? sh : 0); pm_set(s, A.r, A.s, A.d); }
                s->act[w++] = A;
            }
        }
        s->na = w;
        if (s->len == s->capP) { s->capP *= 2; s->P = realloc(s->P, s->capP); s->V = realloc(s->V, s->capP); memset(s->V + s->capP / 2, 0, (size_t)s->capP / 2); }
        s->P[s->len++] = (uint8_t)choose;
        /* a path that walks over itself is a loop */
        { uint64_t km = 0; int bad = 0; for (int i = 0; i < K; i++) { uint8_t b = s->P[s->len - K + i]; if (b > 3) { bad = 1; break; } km = (km << 2) | b; }
          if (!bad) { uint64_t h = mix(km) & (seen_cap - 1); while (seen[h] != ~0ULL && seen[h] != km) h = (h + 1) & (seen_cap - 1);
              if (seen[h] == km) { if (++s->rep >= 150) { strcpy(s->why, "loop"); s->len -= 150; break; } } else { seen[h] = km; s->rep = 0; } } }
        discover(s);
    }
}
/* ---------- the gene ends, each by one thread ---------- */
typedef struct { char *gene, *side, *seq; int glen, n; } Seed;
typedef struct { const Seed *sd; int redo, canon, nout; char **cseq, **cmeta; char stats[192]; } Task;
static Seed *seeds_; static int nseeds = 0;
static Task *tasks; static int ntasks = 0, next_task = 0;
static int maxc = 16; static double minmass = 0.0;
static pthread_mutex_t tm = PTHREAD_MUTEX_INITIALIZER;
static void one_seed(Task *T) {
    const Seed *sd = T->sd; const char *line = sd->seq;
    redo = T->redo; glen = sd->glen;
    EX0 = glen > REDO + K ? glen - REDO : glen / 2;
    St *s0 = st_new(glen + flank_target + 64);
    for (int i = 0; i < EX0; i++) s0->P[i] = (uint8_t)code[(unsigned char)line[i]];
    /* reads laid on the seed gene (differences tolerated there), by every K-mer of it */
    for (s0->len = K; s0->len <= EX0; s0->len++) discover(s0);
    if (redo > 0 && glen - redo >= K + 20) EX0 = glen - redo;      /* the path starts here; what the reads say from here on is taken */
    s0->len = EX0;
    { int w = 0; for (int a = 0; a < s0->na; a++) { Act A = s0->act[a]; if (s0->len - A.d < rd[A.r].len) s0->act[w++] = A; } s0->na = w; }
    St *q[MAXQ]; int nq = 0, done = 0; q[nq++] = s0; work = n_verify = n_hot = n_clone = 0; s0->mass = 1.0; s0->pri = 1.0;
    St *fin[MAXQ]; int nfin = 0;
    while (nq && done < maxc && work <= WORKMAX) {
        int bi = 0; for (int i = 1; i < nq; i++) if (q[i]->pri > q[bi]->pri) bi = i;
        if (done && q[bi]->mass < minmass) break;
        St *s = q[bi]; q[bi] = q[--nq];
        extend(s, q, &nq, MAXQ - 1);
        st_slim(s); fin[nfin++] = s; done++;
    }
    for (int i = 0; i < nq; i++) st_free(q[i]);
    T->cseq = calloc((size_t)(nfin ? nfin : 1), sizeof(char *));
    T->cmeta = calloc((size_t)(nfin ? nfin : 1), sizeof(char *));
    if (!T->cseq || !T->cmeta) { fprintf(stderr, "awx: out of memory\n"); exit(2); }
    int outn = 0;
    for (int i = 0; i < nfin; i++) {
        St *s = fin[i]; int dup = 0;
        for (int j = 0; j < nfin && !dup; j++) if (j != i && fin[j]->len >= s->len && (fin[j]->len > s->len || j < i) && !memcmp(fin[j]->P, s->P, s->len)) dup = 1;
        if (!dup && s->len - glen >= 30) {
            char *sq = malloc((size_t)s->len + 1), meta[192];
            if (!sq) { fprintf(stderr, "awx: out of memory\n"); exit(2); }
            for (int p = 0; p < s->len; p++) sq[p] = "ACGTN"[s->P[p] > 3 ? 4 : s->P[p]];
            sq[s->len] = 0;
            snprintf(meta, sizeof meta, "%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", s->len - glen, s->why, s->born, s->pen, s->nvar, s->nfork, s->nthin, s->jfirst >= 0 ? s->jfirst - glen : -99999, s->mass);
            T->cseq[outn] = sq; T->cmeta[outn] = strdup(meta);
            if (!T->cmeta[outn]) { fprintf(stderr, "awx: out of memory\n"); exit(2); }
            outn++;
        }
    }
    T->nout = outn;
    snprintf(T->stats, sizeof T->stats, "\t%d candidates\twork %llu\tverify %llu\thot_kmers %llu\tclones %llu%s\n", outn, (unsigned long long)work, (unsigned long long)n_verify, (unsigned long long)n_hot, (unsigned long long)n_clone, work > WORKMAX ? "\t(work cap)" : "");
    for (int i = 0; i < nfin; i++) st_free(fin[i]);
}
static void *worker(void *arg) {
    (void)arg;
    for (;;) {
        pthread_mutex_lock(&tm);
        while (next_task < ntasks && tasks[next_task].canon != next_task) next_task++;
        int t = next_task < ntasks ? next_task++ : -1;
        pthread_mutex_unlock(&tm);
        if (t < 0) break;
        one_seed(&tasks[t]);
    }
    return NULL;
}
static size_t sort_next = 0;
static void *sort_worker(void *arg) {
    (void)arg;
    const size_t CHUNK = 1 << 14, TOT = (size_t)1 << PB;
    for (;;) {
        size_t p0 = __sync_fetch_and_add(&sort_next, CHUNK);
        if (p0 >= TOT) break;
        size_t p1 = p0 + CHUNK < TOT ? p0 + CHUNK : TOT;
        for (size_t p = p0; p < p1; p++) {
            size_t a = pt[p], b = pt[p + 1];
            if (b - a < 2) continue;
            if (b - a > 32) { qsort(ix + a, b - a, sizeof(Ent), cmp_ent); continue; }
            for (size_t i = a + 1; i < b; i++) {            /* few entries: insertion keeps the order of read and offset among equal K-mers */
                Ent e = ix[i]; size_t j = i;
                while (j > a && ix[j - 1].km > e.km) { ix[j] = ix[j - 1]; j--; }
                ix[j] = e;
            }
        }
    }
    return NULL;
}
int main(int argc, char **argv) {
    int c, nt = 1; const char *rspec = "0";
    while ((c = getopt(argc, argv, "c:D:r:m:t:vP:Q:")) != -1) {
        if (c == 'c') maxc = atoi(optarg); else if (c == 'P') PATCH = atoi(optarg); else if (c == 'Q') QTOL = atoi(optarg); else if (c == 'D') flank_target = atoi(optarg); else if (c == 'r') rspec = optarg;
        else if (c == 't') nt = atoi(optarg);
        else if (c == 'm') minmass = atof(optarg);          /* arms that carry a smaller share of the gene's reads are not extended */
        else if (c == 'v') { printf("awx %s\n", AWX_VERSION); return 0; }
        else { argc = 0; break; }
    }
    if (argc - optind != 3 || maxc < 1 || flank_target < 30) {
        fprintf(stderr, "awx %s\nusage: awx [-c 16] [-D 10000] [-r 0,300] [-t 1] [-P 1] [-Q 4] seeds.fa out_prefix reads.tsv   (reads.tsv: - for standard input; -P 0 -Q 0: the proposer of 0.3.0)\n", AWX_VERSION);
        return 1;
    }
    const char *seeds = argv[optind], *outp = argv[optind + 1], *reads = argv[optind + 2];
    if (maxc > MAXQ - 2) maxc = MAXQ - 2;
    memset(code, 4, sizeof code); code['A'] = code['a'] = 0; code['C'] = code['c'] = 1; code['G'] = code['g'] = 2; code['T'] = code['t'] = 3;
    FILE *fr = strcmp(reads, "-") ? fopen(reads, "r") : stdin; if (!fr) { fprintf(stderr, "awx: cannot open %s: ", reads); perror(NULL); return 1; }
    double t_start = now();
    /* the reads */
    size_t caps = 1 << 24, ns = 0, capr = 1 << 16; rs = malloc(caps); rd = malloc(capr * sizeof(Read));
    char *line = NULL; size_t cap = 0; ssize_t n; long prevpair = -1; int prevmate = -1;
    while ((n = getline(&line, &cap, fr)) > 0) {
        char *p = line; long pair = strtol(p, &p, 10); int mate = (int)strtol(p + 1, &p, 10) - 1; char *seq = p + 1;
        char *t1 = strchr(seq, '\t'); if (!t1 || mate < 0 || mate > 1) continue; int L = (int)(t1 - seq);
        char *t2 = strchr(t1 + 1, '\t'); char *qual = t2 ? t2 + 1 : NULL; if (L > 1000) L = 1000;
        /* mates must sit at indices 2i, 2i+1: a read without its mate gets an empty partner */
        if (mate == 0 || pair != prevpair || prevmate != 0) { if (nr & 1) { if (nr == capr) { capr *= 2; rd = realloc(rd, capr * sizeof(Read)); } rd[nr++] = (Read){ ns, 0 }; } }
        if (mate == 1 && (pair != prevpair || prevmate != 0)) { if (nr == capr) { capr *= 2; rd = realloc(rd, capr * sizeof(Read)); } rd[nr++] = (Read){ ns, 0 }; }
        if (ns + L + 1 > caps) { while (ns + L + 1 > caps) caps *= 2; rs = realloc(rs, caps); }
        if (nr == capr) { capr *= 2; rd = realloc(rd, capr * sizeof(Read)); }
        int ql = 0; if (qual) while (qual[ql] && qual[ql] != '\n' && qual[ql] != '\t') ql++;
        for (int i = 0; i < L; i++) rs[ns + i] = (uint8_t)code[(unsigned char)seq[i]] | ((qual && ql >= L && qual[i] < 33 + QMIN) ? 8 : 0);
        rd[nr++] = (Read){ ns, (uint16_t)L }; ns += L; prevpair = pair; prevmate = mate;
    }
    if (nr & 1) { if (nr == capr) { capr *= 2; rd = realloc(rd, capr * sizeof(Read)); } rd[nr++] = (Read){ ns, 0 }; }
    /* index: K-mers at every SAMP-th offset of every read, and the last one */
    double t_load = now();
    /* in two goes over the reads: how many entries each first-12-bases bucket gets, then each entry into its bucket. Within a
       bucket the entries then lie in the order of the reads, and are sorted by K-mer, read and offset. */
    pt = calloc(((size_t)1 << PB) + 2, 4);
    if (!pt || !rs || !rd) { fprintf(stderr, "awx: out of memory while loading the reads\n"); return 2; }
    const uint64_t kmask = (1ULL << (2 * K)) - 1;
    for (int go = 0; go < 2; go++) {
        for (uint32_t r = 0; r < nr; r++) {
            int L = rd[r].len; if (L < K) continue;
            const uint8_t *rp = rs + rd[r].off;
            uint64_t km = 0; int valid = 0;
            for (int p = 0; p < L; p++) {
                uint8_t b = rp[p] & 7;
                if (b > 3) { valid = 0; km = 0; } else { km = ((km << 2) | b) & kmask; valid++; }
                int o = p - (K - 1);
                if (o >= 0 && ((o & (SAMP - 1)) == 0 || o == L - K) && valid >= K) {
                    if (go) ix[pt[(km >> (2 * K - PB)) + 1]++] = (Ent){ km, r, (uint16_t)o };
                    else { pt[(km >> (2 * K - PB)) + 2]++; nix++; }
                }
            }
        }
        if (go) break;
        if (nix >= 0xFFFFFFF0ULL) { fprintf(stderr, "awx: too many reads for one index\n"); return 2; }
        for (size_t p = 0; p < ((size_t)1 << PB); p++) pt[p + 2] += pt[p + 1];      /* pt[p + 1]: where bucket p begins; it runs on as the bucket fills */
        ix = malloc((nix + 1) * sizeof(Ent));
        if (!ix) { fprintf(stderr, "awx: out of memory (the reads near the genes and their index need about %.1f GB)\n", (ns + nix * sizeof(Ent) + (double)nr * sizeof(Read)) / 1073741824.0); return 2; }
    }
    /* now pt[p] is where bucket p begins and pt[p + 1] where it ends: independent buckets, sorted across threads */
    if (nt < 1) nt = 1; if (nt > 64) nt = 64;
    sort_next = 0;
    if (nt == 1) sort_worker(NULL);
    else {
        pthread_t sth[64];
        for (int i = 0; i < nt; i++) if (pthread_create(&sth[i], NULL, sort_worker, NULL)) { perror("awx: cannot start a thread"); return 1; }
        for (int i = 0; i < nt; i++) pthread_join(sth[i], NULL);
    }
    double t_index = now();
    { double per = 5.0 * (double)nr + 4096.0; pool_cap = (int)(600e6 / per); if (pool_cap < 16) pool_cap = 16; if (pool_cap > MAXQ - 1) pool_cap = MAXQ - 1; }
    if (fr != stdin) fclose(fr);
    fprintf(stderr, "awx: %u reads, %zu index entries, pool of %d waiting arms; reading %.1f s, index %.1f s\n", nr, nix, pool_cap, t_load - t_start, t_index - t_load);
    /* the gene ends */
    {
        FILE *f = fopen(seeds, "r"); if (!f) { fprintf(stderr, "awx: cannot open %s: ", seeds); perror(NULL); return 1; }
        char name[1024] = ""; int caps_ = 0;
        while ((n = getline(&line, &cap, f)) > 0) {
            while (n && (line[n - 1] == '\n' || line[n - 1] == '\r')) line[--n] = 0;
            if (line[0] == '>') { strncpy(name, line + 1, sizeof name - 1); name[sizeof name - 1] = 0; continue; }
            if (!name[0] || n < K + 10) continue;
            char *b1 = strchr(name, '|'), *b2 = b1 ? strchr(b1 + 1, '|') : NULL; if (!b2) continue;
            *b1 = 0; *b2 = 0;
            if (nseeds == caps_) { caps_ = caps_ ? caps_ * 2 : 64; seeds_ = realloc(seeds_, (size_t)caps_ * sizeof(Seed)); }
            Seed *sd = &seeds_[nseeds++]; sd->gene = strdup(name); sd->side = strdup(b1 + 1); sd->seq = strdup(line); sd->n = (int)n;
            sd->glen = atoi(b2 + 1); if (sd->glen > n) sd->glen = (int)n;
            name[0] = 0;
        }
        fclose(f);
    }
    /* every gene end with every setting (the gene's last `redo` bases taken from the sample's reads) is one piece of work */
    int redos[16], nredo = 0;
    { char *rs_ = strdup(rspec), *tok, *sv = NULL;
      for (tok = strtok_r(rs_, ",", &sv); tok && nredo < 16; tok = strtok_r(NULL, ",", &sv)) { redos[nredo] = atoi(tok); if (redos[nredo] < 0) redos[nredo] = 0; nredo++; }
      free(rs_); }
    ntasks = nredo * nseeds; tasks = calloc((size_t)(ntasks ? ntasks : 1), sizeof(Task));
    int nuniq = 0;
    for (int r = 0; r < nredo; r++) for (int i = 0; i < nseeds; i++) {
        int t = r * nseeds + i, c0 = t;
        tasks[t].sd = &seeds_[i]; tasks[t].redo = redos[r];
        for (int j = 0; j < i; j++) {
            if (seeds_[j].glen == seeds_[i].glen && seeds_[j].n == seeds_[i].n && !strcmp(seeds_[j].seq, seeds_[i].seq)) { c0 = r * nseeds + j; break; }
        }
        tasks[t].canon = c0;
        if (c0 == t) nuniq++;
    }
    if (nt > nuniq) nt = nuniq ? nuniq : 1;
    pthread_t th[64];
    for (int i = 0; i < nt; i++) if (pthread_create(&th[i], NULL, worker, NULL)) { perror("awx: cannot start a thread"); return 1; }
    for (int i = 0; i < nt; i++) pthread_join(th[i], NULL);
    for (int r = 0; r < nredo; r++) {
        char fn[4096];
        snprintf(fn, sizeof fn, "%s.r%d.fa", outp, redos[r]); FILE *fo = fopen(fn, "w"); if (!fo) { fprintf(stderr, "awx: cannot write %s: ", fn); perror(NULL); return 1; }
        snprintf(fn, sizeof fn, "%s.r%d.tsv", outp, redos[r]); FILE *ft = fopen(fn, "w"); if (!ft) { fprintf(stderr, "awx: cannot write %s: ", fn); perror(NULL); return 1; }
        fprintf(ft, "gene\tside\tcand\tflank\tstop\tborn_at\tpenalty\tvariants\tforks\tone_read\tjoin_first\tmass_end\n");
        for (int i = 0; i < nseeds; i++) {
            const Seed *sd = &seeds_[i];
            Task *T = &tasks[tasks[r * nseeds + i].canon];
            for (int k = 0; k < T->nout; k++) {
                fprintf(fo, ">%s|%s|%d|%d|0|0\n%s\n", sd->gene, sd->side, k + 1, sd->glen, T->cseq[k]);
                fprintf(ft, "%s\t%s\t%d\t%s", sd->gene, sd->side, k + 1, T->cmeta[k]);
            }
            fprintf(stderr, "awx: redo %d\t%s\t%s%s", redos[r], sd->gene, sd->side, T->stats);
        }
        for (int i = 0; i < nseeds; i++) {
            int t = r * nseeds + i;
            if (tasks[t].canon == t) {
                for (int k = 0; k < tasks[t].nout; k++) { free(tasks[t].cseq[k]); free(tasks[t].cmeta[k]); }
                free(tasks[t].cseq); free(tasks[t].cmeta);
            }
        }
        int bad = ferror(fo) | ferror(ft);
        if (fclose(fo) | fclose(ft) | bad) { fprintf(stderr, "awx: error writing %s.r%d.* (disk full?)\n", outp, redos[r]); return 1; }
    }
    fprintf(stderr, "awx: %d gene ends, %d settings, %d threads, %.1f s\n", nseeds, nredo, nt, now() - t_index);
    return 0;
}
