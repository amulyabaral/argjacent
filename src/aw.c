// aw.c - "anchor walk": the reads near genes, for any flank length and all genes at once, in two passes over the reads
//
// Pass 1  Every read -> its minimizers ("anchors": about one k-mer per (w+1)/2 bases).
//         Anchors seen at least twice become nodes of a map (a Bloom filter drops singletons).
//         Neighbouring anchors inside a read become a link with their exact distance.
//         Result: a thin map of the whole sample, far smaller than a k-mer table.
//         If the map does not fit in the memory it is given (-m), the pass is run again in parts: each part reads the
//         files once more and builds the anchors of one share of the hash space, which is written to a file; the walk
//         then reads the map from that file. The map is the same, link for link, however it was built.
// Walk    From each gene's anchors, go outward (never turning back) up to D bp per side. Every anchor reached gets its
//         shortest distance from the gene. Each branch carries the share of the gene's reads that goes that way (links
//         split a share in proportion to their read counts, shares add up again where branches meet); a branch carrying
//         less than -p of the share still on its way at that distance is not followed. This bounds the walk in
//         repeats: at any distance at most about 1/p branches are alive.
// Pass 2  Keep read pairs that contain a reached anchor.
// Then    awx proposes contexts from the kept reads; awv and the rules say how far each is certain (see argjacent/).
//
// Usage: aw [-k 31] [-w 20] [-D 10000] [-b 28] [-e 1] [-p 0.002] [-M 2000000] [-B] [-t 1] [-m MB] [-P parts] [-S map] [-G map] [-f]
//           -g genes.fa -o out reads_1.fq[.gz] reads_2.fq[.gz]
//        aw [-k 31] [-w 20] -A seqs.fa          (just print the anchors of some sequences)
//   -p  minimum share of the gene's reads a branch must carry to be followed (0 = follow everything)
//   -M  hard limit of anchors reached per gene side
//   -m  memory for building the map, in MB (0 = no limit, the default): if it is not enough, the map is built in parts
//   -P  build the map in this many parts whatever its size (for tests)
//   -S  save the pass-1 map of the sample to a file      -G  use a saved map instead of running pass 1
//   -f  also write the kept read pairs as FASTQ (out_1.fq, out_2.fq)
//   -t  threads for the read passes (default 1; the result is the same for any number)
//   -b  log2 of Bloom filter bits (28 = 32 MB)     -B  no Bloom filter (keep singleton anchors)
//   -e  min reads supporting a link used in the walk
// Output: out.reads.tsv (pair, mate, read, an empty field, base qualities; mates on adjacent lines), out.genes.tsv (per
//         gene: its anchors, how many were found, how far the walk got), out.stats.txt
// Build: cc -O3 -o aw aw.c -lz -lpthread -lm      (macOS: needs Xcode Command Line Tools)
// Paired 4-line FASTQ, plain or gzip.
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <zlib.h>
#include <pthread.h>
#include <errno.h>
#include <signal.h>
#include <fcntl.h>
#include <math.h>
#include <sys/mman.h>
#include <sys/stat.h>

#define AW_VERSION "0.6.0"
#define NONE 0xFFFFFFFFu
#define ENONE 0x7FFFFFFFu                          /* end of a list of links (31 bits; the top bit of `next` is a side) */
#define SLACK 800
#define MAXA 8192

static int K = 31, W = 20, D = 10000, MINEDGE = 1, use_bloom = 1, MATE_INS = 380;
static uint32_t MAXBALL = 2000000;
static int NTHREADS = 1;
static int LOWC = 1;                               /* 0 (-L): keep low-complexity k-mers as anchors, as before 0.2.0 */
static double EPS = 0.002;          /* minimum share of the gene's reads for a branch to be followed */
static uint64_t KMASK;

/* The map while it is built. A link names the anchor at its other end by that anchor's hash, not by a place in memory,
   so that the map can be built one share of the hash space at a time. `next`: the next link of the same anchor, its
   top bit is the side of the other end; `dist`: distance in bp, its top bit is the side of this end. */
typedef struct __attribute__((packed)) { uint64_t h; uint32_t head; } Node;            /* 12 bytes */
typedef struct { uint64_t to; uint32_t next; uint16_t dist; uint16_t count; } Edge;     /* 16 bytes */
typedef struct { uint64_t h; uint32_t pos; uint8_t st; } Anc;

/* anchors and links are held in pieces of four million: growing never copies, and nothing is taken that is not used */
#define CHB 22
#define CHN (1u << CHB)
static Node *nodec[512]; static uint32_t n_nodes = 0;
static Edge *edgec[512]; static uint32_t n_edges = 0;
#define NODE(i) (nodec[(i) >> CHB][(i) & (CHN - 1)])
#define EDGE(i) (edgec[(i) >> CHB][(i) & (CHN - 1)])
#define HUB_BIT 0x80000000u                        /* top bit of Node.head: anchor has > HUB_DEG links, indexed in hub_ht */
#define HEAD(u) (NODE(u).head & ~HUB_BIT)
#define HUB_DEG 16
typedef struct { uint32_t u1, e; } HubEnt;         /* 8 bytes: node + 1 (0 = empty), link index */
static HubEnt *hub_ht = NULL; static uint64_t hub_cap = 0, hub_n = 0;
static uint32_t *idx; static uint64_t idx_size = 0;
static uint64_t *bloom, *bloom1; static uint64_t bloom_mask = 0;
static int8_t code[256];
static uint32_t NPART = 1, PART = 0;               /* the map is built in NPART parts; this pass builds part PART */
#define MINE(h) (NPART == 1 || (uint32_t)(((h) >> 40) % NPART) == PART)
static uint64_t MEM_BUDGET = 0;                    /* bytes the map may take while it is built (0: no limit) */
static volatile int p1_over = 0;                   /* 1: the part being built does not fit; 2: it has too many anchors or links for one table */

/* ---------- small helpers ---------- */
static void die_oom(size_t n) {
    fprintf(stderr, "aw: out of memory (asked for %.1f MB). Try a larger -w (fewer anchors) or a smaller sample.\n", n / 1048576.0);
    exit(2);
}
static void *xmalloc(size_t n) { void *p = malloc(n ? n : 1); if (!p) die_oom(n); return p; }
static void *xcalloc(size_t n, size_t s) { void *p = calloc(n ? n : 1, s ? s : 1); if (!p) die_oom(n * s); return p; }
static void *xrealloc(void *q, size_t n) { void *p = realloc(q, n ? n : 1); if (!p) die_oom(n); return p; }
static FILE *xfopen(const char *fn, const char *mode) {
    FILE *f = fopen(fn, mode);
    if (!f) { fprintf(stderr, "aw: cannot open %s: ", fn); perror(NULL); exit(1); }
    if (mode[0] == 'w') setvbuf(f, NULL, _IOFBF, 1 << 20);
    return f;
}
static void xfclose(FILE *f, const char *what) {
    int bad = ferror(f);
    if (fclose(f) || bad) { fprintf(stderr, "aw: error writing %s (disk full?)\n", what); exit(1); }
}
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
static long rss_peak_mb(void) {
    struct rusage u;
    if (getrusage(RUSAGE_SELF, &u)) return 0;
#ifdef __APPLE__
    return (long)(u.ru_maxrss / (1024 * 1024));      /* macOS reports bytes */
#else
    return (long)(u.ru_maxrss / 1024);               /* Linux reports kilobytes */
#endif
}
static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x;
}

/* ---------- node table ---------- */
/* The large tables come straight from the system and go straight back to it (the memory allocator may keep what is
   freed for a while, and a map built in parts must not pay for the part before). They come filled with zeros. */
static void *big_alloc(size_t n) {
    void *p = mmap(NULL, n ? n : 1, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
    if (p == MAP_FAILED) die_oom(n);
    return p;
}
static void big_free(void *p, size_t n) { if (p) munmap(p, n ? n : 1); }
static void idx_rebuild(uint64_t size) {
    big_free(idx, idx_size * 4); idx = big_alloc(size * 4); idx_size = size;
    for (uint32_t i = 0; i < n_nodes; i++) {
        uint64_t s = NODE(i).h & (size - 1);
        while (idx[s]) s = (s + 1) & (size - 1);
        idx[s] = i + 1;
    }
}
static inline uint32_t node_find(uint64_t h) {
    uint64_t s = h & (idx_size - 1);
    while (idx[s]) { if (NODE(idx[s] - 1).h == h) return idx[s] - 1; s = (s + 1) & (idx_size - 1); }
    return NONE;
}
static inline uint32_t node_get(uint64_t h) {
    uint64_t s = h & (idx_size - 1);
    while (idx[s]) { if (NODE(idx[s] - 1).h == h) return idx[s] - 1; s = (s + 1) & (idx_size - 1); }
    if ((n_nodes & (CHN - 1)) == 0) {
        if (n_nodes >= 0x40000000u) { p1_over = 2; return n_nodes - 1; }        /* more anchors than one table holds: the map is built in parts */
        nodec[n_nodes >> CHB] = big_alloc((size_t)CHN * sizeof(Node));
    }
    NODE(n_nodes).h = h; NODE(n_nodes).head = ENONE;
    idx[s] = ++n_nodes;
    if ((uint64_t)n_nodes * 10 > idx_size * 7) idx_rebuild(idx_size * 2);       /* the table of places is at most seven tenths full */
    return n_nodes - 1;
}
static inline int bloom_tas_on(uint64_t *bl, uint64_t h) {
    uint64_t *w = &bl[(h >> 20) & bloom_mask];
    uint64_t bits = (1ULL << (h & 63)) | (1ULL << ((h >> 6) & 63)) | (1ULL << ((h >> 12) & 63));
    int seen = (*w & bits) == bits; *w |= bits; return seen;
}
static inline void bloom_set_on(uint64_t *bl, uint64_t h) {
    uint64_t *w = &bl[(h >> 20) & bloom_mask];
    *w |= (1ULL << (h & 63)) | (1ULL << ((h >> 6) & 63)) | (1ULL << ((h >> 12) & 63));
}
static inline int bloom_has(uint64_t h) {
    const uint64_t *w = &bloom[(h >> 20) & bloom_mask];
    uint64_t bits = (1ULL << (h & 63)) | (1ULL << ((h >> 6) & 63)) | (1ULL << ((h >> 12) & 63));
    return (*w & bits) == bits;
}
static inline uint64_t hub_hash(uint32_t u, uint64_t to, uint16_t ds, uint8_t sv) {
    return mix64(to ^ ((uint64_t)u * 0x9e3779b97f4a7c15ULL) ^ ((uint64_t)ds << 1) ^ (uint64_t)sv);
}
static void hub_rebuild(uint64_t cap) {
    HubEnt *old = hub_ht; uint64_t ocap = hub_cap;
    hub_ht = big_alloc(cap * sizeof(HubEnt)); hub_cap = cap;
    const uint64_t m = cap - 1;
    for (uint64_t i = 0; i < ocap; i++) if (old[i].u1) {
        const Edge *E = &EDGE(old[i].e);
        uint64_t s = hub_hash(old[i].u1 - 1, E->to, E->dist, (uint8_t)(E->next >> 31)) & m;
        while (hub_ht[s].u1) s = (s + 1) & m;
        hub_ht[s] = old[i];
    }
    big_free(old, ocap * sizeof(HubEnt));
}
static inline void hub_insert(uint32_t u, uint32_t e, uint64_t h) {
    if (__builtin_expect(!hub_cap, 0)) hub_rebuild(1 << 16);
    else if (__builtin_expect(hub_n * 10 > hub_cap * 7, 0)) hub_rebuild(hub_cap * 2);
    const uint64_t m = hub_cap - 1;
    uint64_t s = h & m;
    while (hub_ht[s].u1) s = (s + 1) & m;
    hub_ht[s] = (HubEnt){ u + 1, e };
    hub_n++;
}
static inline void edge_half(uint32_t u, uint64_t to, uint8_t su, uint8_t sv, uint16_t d, uint8_t mate) {
    /* a link of anchor u to the anchor with hash `to`, d bp away; the same link seen again is counted.
       Bit 14 of `dist` marks a paired-end mate bridge (used by the walk only at dead ends without within-read links).
       Anchors with more than HUB_DEG links (e.g. chimeric poly-G/adapter junction k-mers) are indexed in hub_ht so
       high-degree hubs stay O(1) while keeping the exact same linked list on NODE(u). */
    const uint16_t ds = (uint16_t)((d & 0x3FFFu) | ((uint16_t)mate << 14) | ((uint16_t)su << 15));
    uint64_t hh = 0;
    if (__builtin_expect((NODE(u).head & HUB_BIT) != 0, 0)) {
        hh = hub_hash(u, to, ds, sv);
        const uint64_t m = hub_cap - 1;
        for (uint64_t s = hh & m; hub_ht[s].u1; s = (s + 1) & m) {
            if (hub_ht[s].u1 == u + 1) {
                Edge *E = &EDGE(hub_ht[s].e);
                if (E->to == to && E->dist == ds && (E->next >> 31) == sv) {
                    if (E->count < 65535) E->count++;
                    return;
                }
            }
        }
    } else {
        uint32_t steps = 0;
        for (uint32_t e = HEAD(u); e != ENONE; e = EDGE(e).next & ENONE) {
            Edge *E = &EDGE(e);
            if (E->to == to && E->dist == ds && (E->next >> 31) == sv) {
                if (E->count < 65535) E->count++;
                return;
            }
            steps++;
        }
        if ((n_edges & (CHN - 1)) == 0 && n_edges >= 0x40000000u) { p1_over = 2; return; }
        if (__builtin_expect(steps >= HUB_DEG, 0)) {
            for (uint32_t e = HEAD(u); e != ENONE; e = EDGE(e).next & ENONE) {
                const Edge *E = &EDGE(e);
                hub_insert(u, e, hub_hash(u, E->to, E->dist, (uint8_t)(E->next >> 31)));
            }
            NODE(u).head |= HUB_BIT;
            hh = hub_hash(u, to, ds, sv);
        }
    }
    if ((n_edges & (CHN - 1)) == 0) {
        if (n_edges >= 0x40000000u) { p1_over = 2; return; }
        edgec[n_edges >> CHB] = big_alloc((size_t)CHN * sizeof(Edge));
    }
    uint32_t ne = n_edges++;
    EDGE(ne) = (Edge){ to, HEAD(u) | ((uint32_t)sv << 31), ds, 1 };
    NODE(u).head = ne | (NODE(u).head & HUB_BIT);
    if (__builtin_expect((NODE(u).head & HUB_BIT) != 0, 0)) hub_insert(u, ne, hh);
}
static uint64_t map_bytes(void) {                  /* what the part of the map being built takes, with the Bloom filter */
    return (uint64_t)n_nodes * sizeof(Node) + (uint64_t)n_edges * sizeof(Edge) + idx_size * 4 + hub_cap * sizeof(HubEnt) + (use_bloom ? (bloom_mask + 1) * 8 : 0);
}
static void map_free(void) {
    for (uint32_t i = 0; i < 512; i++) { big_free(nodec[i], (size_t)CHN * sizeof(Node)); big_free(edgec[i], (size_t)CHN * sizeof(Edge)); nodec[i] = NULL; edgec[i] = NULL; }
    big_free(idx, idx_size * 4); idx = NULL; n_nodes = n_edges = 0; idx_size = 0;
    big_free(hub_ht, hub_cap * sizeof(HubEnt)); hub_ht = NULL; hub_cap = hub_n = 0;
}

/* ---------- minimizers ---------- */
typedef struct { uint64_t *HB; uint8_t *SB; int cap; } Scr;      /* work space of one thread */
static volatile int warned_long = 0;
/* A k-mer that is a run of one base or of a 2- or 3-base unit (AAAA.., ACAC.., ACGACG..) is no anchor. Such k-mers
   come from junk reads (poly-G and poly-A tails, dephased reads) by the thousand, link to themselves at a distance
   of 1 to 3, and tie unrelated places together. A real repeat of that kind is at least 31 bp of microsatellite:
   it is left without anchors, like a run of N. */
static inline int low_complexity(const char *s) {
    for (int p = 1; p <= 3; p++) {
        int i = p;
        while (i < K && (s[i] | 32) == (s[i - p] | 32)) i++;
        if (i == K) return 1;
    }
    return 0;
}
static Scr scr0;
static int minimizers(Scr *sc, const char *s, int len, Anc *out, int maxout) {
    int nk = len - K + 1; if (nk < 1) return 0;
    if (nk > sc->cap) { sc->cap = nk * 2; sc->HB = xrealloc(sc->HB, (size_t)sc->cap * 8); sc->SB = xrealloc(sc->SB, (size_t)sc->cap); }
    uint64_t *HB = sc->HB; uint8_t *SB = sc->SB;
    uint64_t f = 0, r = 0; int valid = 0; const int sh = 2 * (K - 1);
    for (int i = 0; i < len; i++) {
        int c = code[(unsigned char)s[i]];
        if (__builtin_expect(c < 0, 0)) { valid = 0; f = r = 0; c = 0; } else valid++;
        f = ((f << 2) | (uint64_t)c) & KMASK; r = (r >> 2) | ((uint64_t)(3 - c) << sh);
        if (i >= K - 1) {
            uint64_t m = f < r ? f : r;                       /* canonical k-mer, no branch */
            HB[i - K + 1] = valid >= K ? mix64(m) : UINT64_MAX;
            SB[i - K + 1] = (uint8_t)(r < f);
        }
    }
    int w = W < nk ? W : nk, n = 0, mp = -1, last = -1;
    for (int i = 0; i + w <= nk; i++) {
        if (mp < i) { mp = i; for (int j = i + 1; j < i + w; j++) if (HB[j] < HB[mp]) mp = j; }
        else if (HB[i + w - 1] < HB[mp]) mp = i + w - 1;
        if (mp != last && HB[mp] != UINT64_MAX) {
            if (LOWC && low_complexity(s + mp)) { last = mp; continue; }
            if (n < maxout) { out[n].h = HB[mp]; out[n].pos = (uint32_t)mp; out[n].st = SB[mp]; n++; }
            else if (!__sync_lock_test_and_set(&warned_long, 1))
                fprintf(stderr, "aw: warning: a sequence has more than %d anchors; the rest of it is not used (aw is made for short reads)\n", maxout);
            last = mp;
        }
    }
    return n;
}

/* ---------- FASTQ block reader (plain or gzip, via zlib) ---------- */
typedef struct { gzFile fd; const char *fn; char *buf; size_t cap, beg, end; int eof; int piped, rfd, wfd; pthread_t th; int cut;
                 volatile uint64_t coff; uint64_t fsize; } LR;     /* coff: how far into the file (as stored) it has been read */
static int quiet_eof = 0;                          /* pass 2 re-reads the same files: warn only once */
static int have_pass1 = 1;                         /* 0 when the map came from a file (-G) */
static void warn_truncated(const char *fn) {
    if (!quiet_eof) fprintf(stderr, "aw: warning: %s ends in the middle of a gzip stream (truncated download?); using the reads before that point\n", fn);
}
/* A gzip file that was cut off: zlib up to 1.3.1 keeps saying "unexpected end of file" from the read that met the end;
   zlib 1.3.2 says it for that one read only (which still returns data) and nothing at the read that returns 0. So the
   state is looked at after every read and remembered. */
static inline void lr_note_cut(LR *r) { int e; gzerror(r->fd, &e); if (e == Z_BUF_ERROR) r->cut = 1; }
/* with threads: one thread per file reads (and unpacks) it and hands the text on through a pipe */
static void *lr_pump(void *arg) {
    LR *r = arg; char *b = xmalloc(1 << 20);
    for (;;) {
        int n = gzread(r->fd, b, 1 << 20);
        r->coff = (uint64_t)gzoffset(r->fd);
        if (n < 0) { int e; fprintf(stderr, "aw: error reading %s: %s\n", r->fn, gzerror(r->fd, &e)); exit(1); }
        lr_note_cut(r);
        if (n == 0) { if (r->cut) warn_truncated(r->fn); break; }
        int off = 0, stop = 0;
        while (off < n) {
            ssize_t w = write(r->wfd, b + off, (size_t)(n - off));
            if (w < 0) { if (errno == EINTR) continue; stop = 1; break; }    /* reader has closed: the other file ran out */
            off += (int)w;
        }
        if (stop) break;
    }
    close(r->wfd); free(b);
    return NULL;
}
static void lr_open(LR *r, const char *fn) {
    r->fd = gzopen(fn, "rb");
    if (!r->fd) { fprintf(stderr, "aw: cannot open %s: ", fn); perror(NULL); exit(1); }
    gzbuffer(r->fd, 1 << 20);
    r->fn = fn; r->cap = 1 << 22; r->buf = xmalloc(r->cap); r->beg = r->end = 0; r->eof = 0; r->piped = 0; r->cut = 0;
    { struct stat sb; r->coff = 0; r->fsize = stat(fn, &sb) ? 0 : (uint64_t)sb.st_size; }
    if (NTHREADS > 1) {
        int fds[2];
        if (pipe(fds)) { perror("aw: cannot make a pipe"); exit(1); }
#ifdef F_SETPIPE_SZ
        (void)fcntl(fds[1], F_SETPIPE_SZ, 1 << 20);
#endif
        r->rfd = fds[0]; r->wfd = fds[1]; r->piped = 1;
        if (pthread_create(&r->th, NULL, lr_pump, r)) { perror("aw: cannot start a reader thread"); exit(1); }
    }
}
static int lr_record(LR *r, char *ln[4], int len[4]) {
    for (;;) {
        size_t p = r->beg; int k = 0;
        while (k < 4) {
            char *nl = memchr(r->buf + p, '\n', r->end - p);
            if (!nl) { if (r->eof && k == 3 && p < r->end) nl = r->buf + r->end; else break; }
            ln[k] = r->buf + p; len[k] = (int)(nl - (r->buf + p)); p = (size_t)(nl - r->buf) + 1; k++;
        }
        if (k == 4) {
            r->beg = p < r->end ? p : r->end;
            for (int j = 0; j < 4; j++) if (len[j] && ln[j][len[j] - 1] == '\r') len[j]--;     /* Windows line ends */
            if (__builtin_expect(ln[0][0] != '@', 0)) {
                fprintf(stderr, "aw: %s does not look like 4-line FASTQ (a record starts with '%c', not '@')\n", r->fn, ln[0][0]);
                exit(1);
            }
            return 1;
        }
        if (r->eof) return 0;
        memmove(r->buf, r->buf + r->beg, r->end - r->beg); r->end -= r->beg; r->beg = 0;
        if (r->end == r->cap) { r->cap *= 2; r->buf = xrealloc(r->buf, r->cap); }
        size_t want = r->cap - r->end; if (want > (1u << 30)) want = 1u << 30;
        if (r->piped) {
            ssize_t n = read(r->rfd, r->buf + r->end, want);
            if (n < 0) { if (errno == EINTR) continue; perror("aw: error reading from the reader thread"); exit(1); }
            if (n == 0) r->eof = 1; else r->end += (size_t)n;
            continue;
        }
        int n = gzread(r->fd, r->buf + r->end, (unsigned)want);
        r->coff = (uint64_t)gzoffset(r->fd);
        if (n < 0) { int e; fprintf(stderr, "aw: error reading %s: %s\n", r->fn, gzerror(r->fd, &e)); exit(1); }
        lr_note_cut(r);
        if (n == 0) {
            if (r->cut) warn_truncated(r->fn);
            r->eof = 1;
        } else r->end += (size_t)n;
    }
}
static void lr_close(LR *r) {
    if (r->piped) { close(r->rfd); pthread_join(r->th, NULL); }
    gzclose(r->fd); free(r->buf);
}
/* read pairs in step; warn once if one file runs out first */
static int pair_record(LR *r1, LR *r2, char *l1[4], int n1[4], char *l2[4], int n2[4], int *warned) {
    int g1 = lr_record(r1, l1, n1), g2 = lr_record(r2, l2, n2);
    if (g1 && g2) return 1;
    if ((g1 || g2) && !*warned) {
        fprintf(stderr, "aw: warning: %s has more reads than %s; extra reads ignored (are these really a pair of files?)\n",
                g1 ? r1->fn : r2->fn, g1 ? r2->fn : r1->fn);
        *warned = 1;
    }
    return 0;
}

/* ---------- pass 1 ---------- */
static uint64_t st_reads = 0, st_anchors = 0, st_solid = 0;
static int write_fq = 0;
static void links_from(Anc *a, int n, uint32_t *id, uint8_t *solid) {
    /* anchors of one read -> nodes (if seen before, and of the part of the map being built) and exact links between
       neighbours. Whether an anchor has been seen before is asked of the Bloom filter for every anchor in every part,
       so that each part sees exactly the links a map built in one piece has. */
    if (PART == 0) { st_reads++; st_anchors += (uint64_t)n; }
    for (int j = 0; j < n; j++) {
        if (use_bloom) __builtin_prefetch(&bloom[(a[j].h >> 20) & bloom_mask], 0);
        if (MINE(a[j].h)) __builtin_prefetch(&idx[a[j].h & (idx_size - 1)]);
    }
    for (int j = 0; j < n; j++) {
        id[j] = NONE;
        if (!use_bloom || bloom_has(a[j].h)) {
            solid[j] = 1;
            if (PART == 0) st_solid++;
            if (MINE(a[j].h)) id[j] = node_get(a[j].h);
        } else solid[j] = 0;
    }
    int p = -1;
    for (int j = 0; j < n; j++) {
        if (!solid[j]) continue;
        if (p >= 0) {
            uint8_t sa = a[p].st ? 0 : 1, sb = a[j].st ? 1 : 0;   /* side: 1 = after the canonical k-mer, 0 = before it */
            uint32_t gap = a[j].pos - a[p].pos;
            if (gap <= 16000) {                                    /* distances are stored in 14 bits */
                if (id[p] != NONE) edge_half(id[p], a[j].h, sa, sb, (uint16_t)gap, 0);
                if (id[j] != NONE && (a[p].h != a[j].h || sa != sb)) edge_half(id[j], a[p].h, sb, sa, (uint16_t)gap, 0);
            }
        }
        p = j;
    }
}
static void pair_links(Anc *a1, int m1, Anc *a2, int m2) {
    /* one read pair into the map: exact within-read neighbour links, plus a paired-end FR bridge between disjoint mates */
    static uint32_t id1[MAXA], id2[MAXA];
    static uint8_t solid1[MAXA], solid2[MAXA];
    links_from(a1, m1, id1, solid1);
    links_from(a2, m2, id2, solid2);
    if (MATE_INS > K && m1 > 0 && m2 > 0) {
        int p1 = -1, p2 = -1;
        for (int i = m1 - 1; i >= 0; i--) if (solid1[i]) { p1 = i; break; }
        for (int j = m2 - 1; j >= 0; j--) if (solid2[j]) { p2 = j; break; }
        if (p1 >= 0 && p2 >= 0 && a1[p1].pos < (uint32_t)MATE_INS && a2[p2].pos < (uint32_t)MATE_INS) {
            int shared = 0;
            for (int i = 0; i < m1 && !shared; i++)
                for (int j = 0; j < m2; j++)
                    if (a1[i].h == a2[j].h) { shared = 1; break; }
            if (!shared) {
                uint8_t s1 = a1[p1].st ? 0 : 1, s2 = a2[p2].st ? 0 : 1;
                int g = (MATE_INS - K) - (int)a1[p1].pos - (int)a2[p2].pos;
                if (g < 80) g = 80;
                if (g <= 16000) {
                    if (id1[p1] != NONE) edge_half(id1[p1], a2[p2].h, s1, s2, (uint16_t)g, 1);
                    if (id2[p2] != NONE && (a1[p1].h != a2[p2].h || s1 != s2)) edge_half(id2[p2], a1[p1].h, s2, s1, (uint16_t)g, 1);
                }
            }
        }
    }
}

/* ---------- batches of read pairs ----------
   The reads are cut into batches and each batch goes through three steps: (1) take the records from the files,
   (2) work that needs no shared state (the anchors of every read; in pass 2 also the output text), (3) apply, IN
   ORDER (update the map, or write the output). With -t 1 one thread does the steps one after the other. With more
   threads the files are unpacked by a thread each, one thread fills the batches, N threads do step 2 and the main
   thread does step 3. The result is the same, byte for byte, for any number of threads. */
#define BATCH 2048
typedef struct {
    int n, done;                                     /* pairs in the batch; step 2 has finished */
    char *txt; size_t tcap, tlen;                    /* the records (all four lines of both mates) */
    uint32_t off[2 * BATCH][4]; int len[2 * BATCH][4];
    Anc *anc; size_t acap; uint32_t aoff[2 * BATCH + 1];   /* anchors of read i: anc[aoff[i] .. aoff[i+1]) */
    char *ob[3]; size_t ol[3];                       /* pass 2: text for reads.tsv, _1.fq, _2.fq */
    uint64_t first, rec;                             /* number of the first pair; pairs kept */
    uint32_t *lh;                                    /* pass 2: lengths of the kept reads */
} Batch;
typedef void (*WorkFn)(Batch *, Scr *, Anc *, Anc *);
static Batch *ring; static uint64_t NB;
static uint64_t q_fill, q_take, q_apply; static int q_stop;      /* batches filled, taken by a worker, applied */
static pthread_mutex_t qm = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t qc_work = PTHREAD_COND_INITIALIZER, qc_done = PTHREAD_COND_INITIALIZER, qc_room = PTHREAD_COND_INITIALIZER;
static WorkFn batch_work;
static LR *fill_r1, *fill_r2; static int *fill_warned;

static void batch_add(Batch *b, int i, char *ln[4], int len[4]) {
    size_t need = b->tlen; for (int k = 0; k < 4; k++) need += (size_t)len[k];
    if (need > b->tcap) { b->tcap = need * 2 + 65536; b->txt = xrealloc(b->txt, b->tcap); }
    for (int k = 0; k < 4; k++) {
        b->off[i][k] = (uint32_t)b->tlen; b->len[i][k] = len[k];
        memcpy(b->txt + b->tlen, ln[k], (size_t)len[k]); b->tlen += (size_t)len[k];
    }
}
static int batch_fill(Batch *b, LR *r1, LR *r2, int *warned, uint64_t *pairs, int *more) {
    /* -> pairs put into the batch; *more = 0 once the files have ended */
    char *l1[4], *l2[4]; int n1[4], n2[4];
    b->n = 0; b->tlen = 0; b->done = 0; b->first = *pairs; b->rec = 0;
    if (p1_over) { *more = 0; return 0; }            /* the part of the map being built does not fit: stop reading */
    while (b->n < BATCH && b->tlen < ((size_t)1 << 27)) {       /* offsets into the text are 32 bits */
        if (!pair_record(r1, r2, l1, n1, l2, n2, warned)) { *more = 0; break; }
        batch_add(b, 2 * b->n, l1, n1); batch_add(b, 2 * b->n + 1, l2, n2); b->n++;
    }
    *pairs += (uint64_t)b->n;
    return b->n;
}
static void *batch_filler(void *arg) {
    (void)arg; uint64_t pairs = 0; int more = 1;
    while (more) {
        pthread_mutex_lock(&qm);
        while (q_fill - q_apply == NB) pthread_cond_wait(&qc_room, &qm);
        pthread_mutex_unlock(&qm);
        int n = batch_fill(&ring[q_fill % NB], fill_r1, fill_r2, fill_warned, &pairs, &more);
        pthread_mutex_lock(&qm);
        if (n) { q_fill++; pthread_cond_signal(&qc_work); }
        if (!more) { q_stop = 1; pthread_cond_broadcast(&qc_work); pthread_cond_broadcast(&qc_done); }
        pthread_mutex_unlock(&qm);
    }
    return NULL;
}
static void *batch_worker(void *arg) {
    (void)arg;
    Scr sc = { NULL, NULL, 0 };
    Anc *t1 = xmalloc(sizeof(Anc) * MAXA), *t2 = xmalloc(sizeof(Anc) * MAXA);
    for (;;) {
        pthread_mutex_lock(&qm);
        while (q_take == q_fill && !q_stop) pthread_cond_wait(&qc_work, &qm);
        if (q_take == q_fill) { pthread_mutex_unlock(&qm); break; }
        Batch *b = &ring[q_take++ % NB];
        pthread_mutex_unlock(&qm);
        batch_work(b, &sc, t1, t2);
        pthread_mutex_lock(&qm);
        b->done = 1; pthread_cond_broadcast(&qc_done);
        pthread_mutex_unlock(&qm);
    }
    free(sc.HB); free(sc.SB); free(t1); free(t2);
    return NULL;
}
static void run_batches(LR *r1, LR *r2, int *warned, WorkFn work, void (*apply)(Batch *)) {
    if (NTHREADS < 2) {
        Batch *b = xcalloc(1, sizeof(Batch)); uint64_t pairs = 0; int more = 1;
        Anc *t1 = xmalloc(sizeof(Anc) * MAXA), *t2 = xmalloc(sizeof(Anc) * MAXA);
        while (more)
            if (batch_fill(b, r1, r2, warned, &pairs, &more)) { work(b, &scr0, t1, t2); apply(b); }
        free(b->txt); free(b->anc); free(b->lh); free(b); free(t1); free(t2);
        return;
    }
    NB = 2 * (uint64_t)NTHREADS + 6;
    ring = xcalloc((size_t)NB, sizeof(Batch));
    q_fill = q_take = q_apply = 0; q_stop = 0; batch_work = work;
    fill_r1 = r1; fill_r2 = r2; fill_warned = warned;
    pthread_t filler, *th = xmalloc(sizeof(pthread_t) * (size_t)NTHREADS);
    if (pthread_create(&filler, NULL, batch_filler, NULL)) { perror("aw: cannot start a thread"); exit(1); }
    for (int i = 0; i < NTHREADS; i++)
        if (pthread_create(&th[i], NULL, batch_worker, NULL)) { perror("aw: cannot start a thread"); exit(1); }
    for (;;) {
        pthread_mutex_lock(&qm);
        while (q_apply == q_fill && !q_stop) pthread_cond_wait(&qc_done, &qm);
        if (q_apply == q_fill) { pthread_mutex_unlock(&qm); break; }       /* all read, all applied */
        Batch *a = &ring[q_apply % NB];
        while (!a->done) pthread_cond_wait(&qc_done, &qm);
        pthread_mutex_unlock(&qm);
        apply(a);
        pthread_mutex_lock(&qm);
        q_apply++; pthread_cond_signal(&qc_room);
        pthread_mutex_unlock(&qm);
    }
    pthread_join(filler, NULL);
    for (int i = 0; i < NTHREADS; i++) pthread_join(th[i], NULL);
    for (uint64_t i = 0; i < NB; i++) { free(ring[i].txt); free(ring[i].anc); free(ring[i].lh); }
    free(ring); free(th);
}
static LR *p1_r1;
static double p1_frac = 1.0, p1_half = 0;           /* how much of the first file had been read when the part being built stopped fitting, and when it was half full */
static double p1_now(void) { double f = p1_r1 && p1_r1->fsize ? (double)p1_r1->coff / (double)p1_r1->fsize : 1.0; return f > 1 ? 1 : f < 0.001 ? 0.001 : f; }
static void batch_anchors(Batch *b, Scr *sc, Anc *t1, Anc *t2) {
    (void)t2; size_t n = 0;
    for (int i = 0; i < 2 * b->n; i++) {
        int m = minimizers(sc, b->txt + b->off[i][1], b->len[i][1], t1, MAXA);
        if (n + (size_t)m > b->acap) { b->acap = (n + (size_t)m) * 2 + 4096; b->anc = xrealloc(b->anc, b->acap * sizeof(Anc)); }
        memcpy(b->anc + n, t1, (size_t)m * sizeof(Anc));
        b->aoff[i] = (uint32_t)n; n += (size_t)m;
    }
    b->aoff[2 * b->n] = (uint32_t)n;
}
static void apply_links(Batch *b) {
    /* The anchors of the whole batch are known, so the memory they will touch is fetched ahead of use, in steps:
       the Bloom word and the index slot; then the anchor's record; then its first link; then the second. Hints only. */
    if (p1_over) return;                             /* this part does not fit: the rest of the pass is only read off */
    const Anc *A = b->anc; const uint32_t tot = b->aoff[2 * b->n];
    uint32_t c1 = 0, c2 = 0, c3 = 0, c4 = 0;
    for (int p = 0; p < b->n; p++) {
        uint32_t lo = b->aoff[2 * p], mid = b->aoff[2 * p + 1], hi = b->aoff[2 * p + 2];
        uint32_t e1 = hi + 128 < tot ? hi + 128 : tot, e2 = hi + 80 < tot ? hi + 80 : tot, e3 = hi + 40 < tot ? hi + 40 : tot;
        uint32_t e4 = hi + 12 < tot ? hi + 12 : tot;
        for (; c1 < e1; c1++) {
            if (use_bloom) __builtin_prefetch(&bloom[(A[c1].h >> 20) & bloom_mask], 1);
            if (MINE(A[c1].h)) __builtin_prefetch(&idx[A[c1].h & (idx_size - 1)]);
        }
        for (; c2 < e2; c2++) { if (!MINE(A[c2].h)) continue; uint32_t i = idx[A[c2].h & (idx_size - 1)]; if (i) __builtin_prefetch(&NODE(i - 1), 1); }
        for (; c3 < e3; c3++) {
            if (!MINE(A[c3].h)) continue;
            uint32_t i = idx[A[c3].h & (idx_size - 1)];
            if (i && HEAD(i - 1) != ENONE) __builtin_prefetch(&EDGE(HEAD(i - 1)), 1);
        }
        for (; c4 < e4; c4++) {
            if (!MINE(A[c4].h)) continue;
            uint32_t i = idx[A[c4].h & (idx_size - 1)], e;
            if (i && (e = HEAD(i - 1)) != ENONE && (EDGE(e).next & ENONE) != ENONE) __builtin_prefetch(&EDGE(EDGE(e).next & ENONE), 1);
        }
        pair_links(b->anc + lo, (int)(mid - lo), b->anc + mid, (int)(hi - mid));
    }
    if (MEM_BUDGET && !p1_over) {
        const uint64_t mb = map_bytes(), fixed = use_bloom ? (bloom_mask + 1) * 8 : 0;
        if (p1_half == 0 && (mb - fixed) * 2 > MEM_BUDGET - fixed) p1_half = p1_now();      /* half full: to tell how the map grows with the reads */
        if (mb > MEM_BUDGET) p1_over = 1;
    }
}
static void apply_bloom(Batch *b) {
    const Anc *A = b->anc; const uint32_t tot = b->aoff[2 * b->n];
    for (uint32_t c = 0; c < tot; c++) {
        if (c + 64 < tot) {
            __builtin_prefetch(&bloom1[(A[c + 64].h >> 20) & bloom_mask], 1);
            __builtin_prefetch(&bloom[(A[c + 64].h >> 20) & bloom_mask], 1);
        }
        if (bloom_tas_on(bloom1, A[c].h)) bloom_set_on(bloom, A[c].h);
    }
}
static void pass0(const char *fn1, const char *fn2) {
    LR r1, r2; lr_open(&r1, fn1); lr_open(&r2, fn2);
    int warned = 0; quiet_eof = 0; p1_over = 0;
    run_batches(&r1, &r2, &warned, batch_anchors, apply_bloom);
    lr_close(&r1); lr_close(&r2);
}
static void pass1(const char *fn1, const char *fn2) {
    LR r1, r2; lr_open(&r1, fn1); lr_open(&r2, fn2);
    int warned = use_bloom || PART > 0;              /* pass 0 (or the first part) has given the warnings about the files */
    quiet_eof = use_bloom || PART > 0;
    p1_r1 = &r1; p1_half = 0;
    run_batches(&r1, &r2, &warned, batch_anchors, apply_links);
    if (p1_over) p1_frac = p1_now();
    p1_r1 = NULL;
    lr_close(&r1); lr_close(&r2);
}

/* ---------- the map in a file ----------
   Written when the map is built in parts or is to be saved (-S); read through the page cache by the walk. Per part:
   the links of all its anchors (those of one anchor next to each other, in the order the walk goes through them),
   the anchors sorted by hash, and a table with the first anchor of every LMSTEP-th and the place of its first link. */
typedef struct { uint64_t h; uint32_t deg; uint32_t pad; } NodeD;
typedef struct { uint64_t to; uint16_t dist; uint16_t count; uint8_t su, sv, mate, pad; } EdgeD;
typedef struct { uint64_t h; uint64_t eoff; } Lm;
#define LMSTEP 512
typedef struct { uint64_t n_nodes, n_edges, n_lm, off_edges, off_nodes, off_lm; } PartDir;
typedef struct { char magic[8]; int32_t k, w; uint32_t nparts, pad; uint64_t reads, anchors, solid, off_dir; } MapHead;
#define MAP_MAGIC "AWMAP03"
typedef char check_sizes[(sizeof(Node) == 12 && sizeof(NodeD) == 16 && sizeof(Edge) == 16 && sizeof(EdgeD) == 16 && sizeof(Lm) == 16 && sizeof(PartDir) == 48 && sizeof(MapHead) == 56) ? 1 : -1];
static FILE *mapf; static PartDir *mdir; static const char *mapf_name;
static void map_wfail(void) { fprintf(stderr, "aw: error writing %s (disk full?)\n", mapf_name); exit(1); }
static void fpad16(FILE *f) { long long o = (long long)ftello(f); while (o & 15) { fputc(0, f); o++; } }
static void map_begin(const char *fn, uint32_t nparts) {
    mapf = xfopen(fn, "wb"); mapf_name = fn;
    MapHead h; memset(&h, 0, sizeof h);
    if (fwrite(&h, sizeof h, 1, mapf) != 1) map_wfail();
    mdir = xcalloc(nparts, sizeof(PartDir));
}
static int cmp_u64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }
static void map_part_save(uint32_t p) {
    /* the part of the map in memory -> the file, anchors in the order of their hashes. The table of places is not needed
       any more and makes room first; the list heads of the anchors are replaced by the numbers of their links on the way. */
    big_free(idx, idx_size * 4); idx = NULL; idx_size = 0;
    big_free(hub_ht, hub_cap * sizeof(HubEnt)); hub_ht = NULL; hub_cap = hub_n = 0;
    PartDir *d = &mdir[p]; d->n_nodes = n_nodes; d->n_lm = ((uint64_t)n_nodes + LMSTEP - 1) / LMSTEP;
    /* the anchors' numbers, sorted by hash: first into buckets by the top bits of the hash (the hashes are spread evenly:
       about 256 per bucket), then each bucket by itself */
    int bb = 0; while (bb < 30 && ((uint64_t)n_nodes >> (bb + 8))) bb++;
    const uint64_t nbk = 1ULL << bb;
    uint32_t *bst = xcalloc(nbk + 2, 4), *ord = xmalloc(((size_t)n_nodes + 1) * 4);
#define BKT(h) (bb ? (h) >> (64 - bb) : 0)
    for (uint32_t i = 0; i < n_nodes; i++) bst[BKT(NODE(i).h) + 2]++;
    for (uint64_t k = 0; k < nbk; k++) bst[k + 2] += bst[k + 1];
    for (uint32_t i = 0; i < n_nodes; i++) ord[bst[BKT(NODE(i).h) + 1]++] = i;     /* now bst[k] is where bucket k begins */
    Lm *lm = xmalloc((size_t)(d->n_lm + 1) * sizeof(Lm));
    fpad16(mapf); d->off_edges = (uint64_t)ftello(mapf);
    const uint32_t BUF = 65536; EdgeD *buf = xmalloc(sizeof(EdgeD) * BUF); uint32_t nb = 0; uint64_t ne = 0, out = 0;
    uint64_t *key = NULL; size_t key_cap = 0;
    for (uint64_t k = 0; k < nbk; k++) {
        const uint32_t a = bst[k], z = bst[k + 1], m = z - a;
        if (!m) continue;
        if (2 * (size_t)m > key_cap) { key_cap = 2 * (size_t)m + 1024; key = xrealloc(key, key_cap * 8); }
        for (uint32_t j = 0; j < m; j++) { key[2 * j] = NODE(ord[a + j]).h; key[2 * j + 1] = ord[a + j]; }
        qsort(key, m, 16, cmp_u64);
        for (uint32_t j = 0; j < m; j++, out++) {
            const uint32_t i = (uint32_t)key[2 * j + 1];
            ord[a + j] = i;
            if (out % LMSTEP == 0) lm[out / LMSTEP] = (Lm){ NODE(i).h, ne };
            uint32_t deg = 0;
            for (uint32_t e = HEAD(i); e != ENONE; e = EDGE(e).next & ENONE) {
                if (nb == BUF) { if (fwrite(buf, sizeof(EdgeD), nb, mapf) != nb) map_wfail(); nb = 0; }
                const Edge *E = &EDGE(e); EdgeD *o = &buf[nb++];
                o->to = E->to; o->dist = E->dist & 0x3FFF; o->count = E->count;
                o->mate = (uint8_t)((E->dist >> 14) & 1); o->su = (uint8_t)(E->dist >> 15); o->sv = (uint8_t)(E->next >> 31); o->pad = 0;
                deg++;
            }
            NODE(i).head = deg; ne += deg;
        }
    }
#undef BKT
    if (nb && fwrite(buf, sizeof(EdgeD), nb, mapf) != nb) map_wfail();
    d->n_edges = ne;
    d->off_nodes = (uint64_t)ftello(mapf);
    NodeD *nbuf = (NodeD *)buf;                       /* the same room, now for the anchors */
    nb = 0;
    for (uint32_t j = 0; j < n_nodes; j++) {
        if (nb == BUF) { if (fwrite(nbuf, sizeof(NodeD), nb, mapf) != nb) map_wfail(); nb = 0; }
        nbuf[nb++] = (NodeD){ NODE(ord[j]).h, NODE(ord[j]).head, 0 };
    }
    if (nb && fwrite(nbuf, sizeof(NodeD), nb, mapf) != nb) map_wfail();
    d->off_lm = (uint64_t)ftello(mapf);
    if (d->n_lm && fwrite(lm, sizeof(Lm), (size_t)d->n_lm, mapf) != d->n_lm) map_wfail();
    free(buf); free(lm); free(bst); free(ord); free(key);
}
static void map_end(uint32_t nparts) {
    fpad16(mapf);
    MapHead h; memset(&h, 0, sizeof h); memcpy(h.magic, MAP_MAGIC, 8);
    h.k = K; h.w = W; h.nparts = nparts; h.reads = st_reads; h.anchors = st_anchors; h.solid = st_solid; h.off_dir = (uint64_t)ftello(mapf);
    if (fwrite(mdir, sizeof(PartDir), nparts, mapf) != nparts) map_wfail();
    if (fseeko(mapf, 0, SEEK_SET) || fwrite(&h, sizeof h, 1, mapf) != 1) map_wfail();
    xfclose(mapf, mapf_name); mapf = NULL; free(mdir); mdir = NULL;
}
/* the map as the walk sees it: from a file (mm_base set), or the table in memory */
static uint8_t *mm_base = NULL; static size_t mm_len = 0; static uint32_t mm_parts = 0; static const PartDir *mm_dir;
static int mm_hint = 0; static uintptr_t page_mask;
static uint64_t tot_nodes = 0, tot_edges = 0;
static void map_bad(const char *fn, const char *why) { fprintf(stderr, "aw: %s %s\n", fn, why); exit(1); }
static void map_open(const char *fn, int is_temp) {
    int fd = open(fn, O_RDONLY);
    if (fd < 0) { fprintf(stderr, "aw: cannot open %s: ", fn); perror(NULL); exit(1); }
    struct stat sb;
    if (fstat(fd, &sb)) { perror("aw: fstat"); exit(1); }
    mm_len = (size_t)sb.st_size;
    MapHead h;
    if (mm_len < sizeof h || read(fd, &h, sizeof h) != (ssize_t)sizeof h || memcmp(h.magic, "AWMAP", 5)) map_bad(fn, "is not a map saved by aw -S");
    if (memcmp(h.magic, MAP_MAGIC, 8)) map_bad(fn, "is a map of another version of aw: make it again with this one");
    if (h.k != K || h.w != W) { fprintf(stderr, "aw: %s was made with -k %d -w %d, not -k %d -w %d\n", fn, h.k, h.w, K, W); exit(1); }
    if (h.nparts < 1 || h.nparts > 65536 || h.off_dir > mm_len || (h.off_dir & 15) || mm_len - h.off_dir < (uint64_t)h.nparts * sizeof(PartDir))
        map_bad(fn, "is shorter than its header says (incomplete copy?)");
    if (mm_len - h.off_dir != (uint64_t)h.nparts * sizeof(PartDir)) map_bad(fn, "is longer than its header says (not a map of this version?)");
    mm_base = mmap(NULL, mm_len, PROT_READ, MAP_SHARED, fd, 0);
    if (mm_base == MAP_FAILED) { perror("aw: cannot map the map file into memory"); exit(1); }
    close(fd);
    if (is_temp) unlink(fn);                         /* its room on the disk is given back when aw ends, however it ends */
    mm_parts = h.nparts; mm_dir = (const PartDir *)(mm_base + h.off_dir);
    st_reads = h.reads; st_anchors = h.anchors; st_solid = h.solid; tot_nodes = tot_edges = 0;
    for (uint32_t p = 0; p < mm_parts; p++) {
        const PartDir *d = &mm_dir[p];
        if ((d->off_edges & 15) || (d->off_nodes & 15) || (d->off_lm & 15) || d->off_edges > h.off_dir || d->off_nodes > h.off_dir || d->off_lm > h.off_dir ||
            d->n_edges > (h.off_dir - d->off_edges) / sizeof(EdgeD) || d->n_nodes > (h.off_dir - d->off_nodes) / sizeof(NodeD) ||
            d->n_lm > (h.off_dir - d->off_lm) / sizeof(Lm) || d->n_lm != (d->n_nodes + LMSTEP - 1) / LMSTEP) map_bad(fn, "is damaged (its table of contents does not fit the file)");
        const Lm *lm = (const Lm *)(mm_base + d->off_lm);
        for (uint64_t i = 0; i < d->n_lm; i++)
            if (lm[i].eoff > d->n_edges || (i && (lm[i].eoff < lm[i - 1].eoff || lm[i].h <= lm[i - 1].h))) map_bad(fn, "is damaged (its index is out of order)");
        tot_nodes += d->n_nodes; tot_edges += d->n_edges;
    }
    page_mask = ~((uintptr_t)sysconf(_SC_PAGESIZE) - 1);
    (void)madvise(mm_base, mm_len, MADV_RANDOM);
    /* a map much larger than what the page cache will hold is read from the disk as the walk asks for it: then the walk
       says ahead which anchors it will ask for */
    { long pages = sysconf(_SC_PHYS_PAGES), psz = sysconf(_SC_PAGESIZE); mm_hint = pages > 0 && (double)mm_len > 0.25 * (double)pages * (double)psz; }
}
static void map_close(void) {
    if (mm_base) { munmap(mm_base, mm_len); mm_base = NULL; }
    map_free();
}
typedef struct { uint32_t deg; const EdgeD *e; } NodeRef;
static EdgeD *ebuf; static uint32_t ebuf_cap = 0;
static inline const Lm *lm_block(uint64_t h, uint32_t *part, uint64_t *block) {
    uint32_t p = mm_parts == 1 ? 0 : (uint32_t)((h >> 40) % mm_parts);
    const PartDir *d = &mm_dir[p]; const Lm *lm = (const Lm *)(mm_base + d->off_lm);
    if (!d->n_lm || h < lm[0].h) return NULL;
    uint64_t lo = 0, hi = d->n_lm;                   /* lm[lo].h <= h < lm[hi].h */
    while (hi - lo > 1) { uint64_t mid = (lo + hi) / 2; if (lm[mid].h <= h) lo = mid; else hi = mid; }
    *part = p; *block = lo; return lm;
}
static int map_node(uint64_t h, NodeRef *r) {
    /* the anchor with hash h: how often it was seen and its links (valid until the next call); 0 if the map has no such anchor */
    if (mm_base) {
        uint32_t p; uint64_t b; const Lm *lm = lm_block(h, &p, &b);
        if (!lm) return 0;
        const PartDir *d = &mm_dir[p]; const NodeD *nd = (const NodeD *)(mm_base + d->off_nodes);
        uint64_t a = b * LMSTEP, z = a + LMSTEP < d->n_nodes ? a + LMSTEP : d->n_nodes, l = a, u = z;
        while (l < u) { uint64_t mid = (l + u) / 2; if (nd[mid].h < h) l = mid + 1; else u = mid; }
        if (l >= z || nd[l].h != h) return 0;
        uint64_t eo = lm[b].eoff; for (uint64_t i = a; i < l; i++) eo += nd[i].deg;
        if (eo > d->n_edges || nd[l].deg > d->n_edges - eo) return 0;      /* a damaged file must not lead outside it */
        r->deg = nd[l].deg; r->e = (const EdgeD *)(mm_base + d->off_edges) + eo;
        return 1;
    }
    uint32_t u = node_find(h); if (u == NONE) return 0;
    uint32_t n = 0;
    for (uint32_t e = HEAD(u); e != ENONE; e = EDGE(e).next & ENONE) {
        if (n == ebuf_cap) { ebuf_cap = ebuf_cap ? ebuf_cap * 2 : 1024; ebuf = xrealloc(ebuf, (size_t)ebuf_cap * sizeof(EdgeD)); }
        const Edge *E = &EDGE(e); EdgeD *o = &ebuf[n++];
        o->to = E->to; o->dist = E->dist & 0x3FFF; o->count = E->count;
        o->mate = (uint8_t)((E->dist >> 14) & 1); o->su = (uint8_t)(E->dist >> 15); o->sv = (uint8_t)(E->next >> 31); o->pad = 0;
    }
    r->deg = n; r->e = ebuf;
    return 1;
}
static int is_chimeric_hub(uint64_t h, uint8_t sv) {
    /* 1 if side `sv` of anchor h fans out into >= 24 within-read links yet none has count > 3
       (signature of a sequencing adapter / poly-tail junction entered from random genomic cut points,
       never a high-depth plasmid or multi-copy IS element which has a strong backbone edge on `sv`). */
    uint32_t deg_sv = 0, max_c = 0;
    if (mm_base) {
        uint32_t p; uint64_t b; const Lm *lm = lm_block(h, &p, &b);
        if (!lm) return 0;
        const PartDir *d = &mm_dir[p]; const NodeD *nd = (const NodeD *)(mm_base + d->off_nodes);
        uint64_t a = b * LMSTEP, z = a + LMSTEP < d->n_nodes ? a + LMSTEP : d->n_nodes, l = a, u = z;
        while (l < u) { uint64_t mid = (l + u) / 2; if (nd[mid].h < h) l = mid + 1; else u = mid; }
        if (l >= z || nd[l].h != h || nd[l].deg < 24) return 0;
        uint64_t eo = lm[b].eoff; for (uint64_t i = a; i < l; i++) eo += nd[i].deg;
        if (eo > d->n_edges || nd[l].deg > d->n_edges - eo) return 0;
        const EdgeD *es = (const EdgeD *)(mm_base + d->off_edges) + eo;
        for (uint32_t i = 0; i < nd[l].deg; i++) if (!es[i].mate && es[i].su == sv) {
            deg_sv++; if (es[i].count > max_c) max_c = es[i].count;
            if (max_c > 3) return 0;
        }
        return deg_sv >= 24;
    }
    uint32_t u = node_find(h);
    if (u == NONE || !(NODE(u).head & HUB_BIT)) return 0;
    for (uint32_t e = HEAD(u); e != ENONE; e = EDGE(e).next & ENONE) {
        const Edge *E = &EDGE(e);
        if (!((E->dist >> 14) & 1) && (E->dist >> 15) == sv) {
            deg_sv++; if (E->count > max_c) max_c = E->count;
            if (max_c > 3) return 0;
        }
    }
    return deg_sv >= 24;
}
static inline void map_hint(uint64_t h) {
    if (!mm_hint) return;
    uint32_t p; uint64_t b; const Lm *lm = lm_block(h, &p, &b);
    if (!lm) return;
    const PartDir *d = &mm_dir[p];
    uintptr_t a = (uintptr_t)(mm_base + d->off_nodes + b * LMSTEP * sizeof(NodeD)) & page_mask;
    (void)madvise((void *)a, LMSTEP * sizeof(NodeD) + (size_t)~page_mask, MADV_WILLNEED);
    uint64_t e0 = lm[b].eoff, e1 = b + 1 < d->n_lm ? lm[b + 1].eoff : d->n_edges;
    a = (uintptr_t)(mm_base + d->off_edges + e0 * sizeof(EdgeD)) & page_mask;
    (void)madvise((void *)a, (size_t)(e1 - e0) * sizeof(EdgeD) + (size_t)~page_mask, MADV_WILLNEED);
}

/* Bloom filter load: expected share of one-off anchors that slip through as if seen twice */
static double bloom_false_rate(const uint64_t *bl) {
    double s = 0; uint64_t nw = bloom_mask + 1;
    for (uint64_t i = 0; i < nw; i++) { double f = __builtin_popcountll(bl[i]) / 64.0; s += f * f * f; }
    return s / (double)nw;
}

/* ---------- a heap of anchors waiting to be placed, nearest first ---------- */
typedef struct { int32_t key; uint32_t node; int32_t coord; uint32_t pred; float delta; uint8_t rev, kind; } HE;
static HE *heap; static uint32_t heap_n = 0, heap_cap = 0;
static void hpush(HE e) {
    if (heap_n == heap_cap) { heap_cap = heap_cap ? heap_cap * 2 : 1024; heap = xrealloc(heap, (size_t)heap_cap * sizeof(HE)); }
    uint32_t i = heap_n++;
    while (i && heap[(i - 1) / 2].key > e.key) { heap[i] = heap[(i - 1) / 2]; i = (i - 1) / 2; }
    heap[i] = e;
}
static HE hpop(void) {
    HE top = heap[0], e = heap[--heap_n]; uint32_t i = 0;
    for (;;) {
        uint32_t c = 2 * i + 1; if (c >= heap_n) break;
        if (c + 1 < heap_n && heap[c + 1].key < heap[c].key) c++;
        if (heap[c].key >= e.key) break;
        heap[i] = heap[c]; i = c;
    }
    heap[i] = e; return top;
}

/* ---------- outward walk ----------
   What the walk knows about an anchor is kept for the anchors it meets only, by hash. state: 0 waiting, 1 followed (in
   the ball), 2 not followed (share too small), 3 not followed before, now waiting again. */
typedef struct { uint64_t h; int32_t coord; float mass; uint8_t placed, hasm, rev, seed, state, redo; } WS;
static WS *ws; static uint32_t ws_n = 0, ws_cap = 0; static uint32_t *wsi; static uint64_t wsi_size = 0;
static uint32_t *ball, *drop; static uint32_t nball, ball_cap = 0, ndrop, drop_cap = 0;
static double Tlive;         /* share that is still on its way: the sum over anchors waiting to be placed */
#define SAMEPOS 200          /* an anchor met again within this many bp of its place is the same place (bubble, indel) */
static void ws_reset(void) {
    if (!wsi_size) { wsi_size = 1 << 16; wsi = xcalloc(wsi_size, 4); } else memset(wsi, 0, wsi_size * 4);
    ws_n = 0;
}
static inline uint32_t ws_find(uint64_t h) {
    uint64_t s = (h >> 9) & (wsi_size - 1);
    while (wsi[s]) { if (ws[wsi[s] - 1].h == h) return wsi[s] - 1; s = (s + 1) & (wsi_size - 1); }
    return NONE;
}
static uint32_t ws_get(uint64_t h) {                /* the anchor's record, made if the walk has not met it yet */
    uint64_t s = (h >> 9) & (wsi_size - 1);
    while (wsi[s]) { if (ws[wsi[s] - 1].h == h) return wsi[s] - 1; s = (s + 1) & (wsi_size - 1); }
    if (ws_n == ws_cap) {
        if (ws_cap >= 0x7FFFFFF0u) { fprintf(stderr, "aw: the walk met too many anchors\n"); exit(2); }
        ws_cap = ws_cap ? (ws_cap < 0x40000000u ? ws_cap * 2 : 0x7FFFFFF0u) : 1 << 14; ws = xrealloc(ws, (size_t)ws_cap * sizeof(WS));
    }
    memset(&ws[ws_n], 0, sizeof(WS)); ws[ws_n].h = h;
    wsi[s] = ++ws_n;
    if ((uint64_t)ws_n * 2 > wsi_size) {
        free(wsi); wsi_size *= 2; wsi = xcalloc(wsi_size, 4);
        for (uint32_t i = 0; i < ws_n; i++) { uint64_t t = (ws[i].h >> 9) & (wsi_size - 1); while (wsi[t]) t = (t + 1) & (wsi_size - 1); wsi[t] = i + 1; }
    }
    map_hint(h);
    return ws_n - 1;
}
/* the anchors reached by any walk: pass 2 keeps the read pairs that hold one of them */
static uint64_t *bset; static uint64_t bset_mask = 0, n_reached = 0;
static inline int bset_has(uint64_t h) {
    uint64_t s = (h >> 7) & bset_mask;
    while (bset[s]) { if (bset[s] == h) return 1; s = (s + 1) & bset_mask; }
    return 0;
}
static void bset_add(uint64_t h) {
    if (!bset) { bset_mask = 1023; bset = xcalloc(bset_mask + 1, 8); }
    uint64_t s = (h >> 7) & bset_mask;
    while (bset[s]) { if (bset[s] == h) return; s = (s + 1) & bset_mask; }
    bset[s] = h; n_reached++;
    if (n_reached * 4 > bset_mask + 1) {
        uint64_t *o = bset, om = bset_mask;
        bset_mask = (bset_mask + 1) * 2 - 1; bset = xcalloc(bset_mask + 1, 8);
        for (uint64_t i = 0; i <= om; i++) if (o[i]) { uint64_t t = (o[i] >> 7) & bset_mask; while (bset[t]) t = (t + 1) & bset_mask; bset[t] = o[i]; }
        free(o);
    }
}

/* pass the share m of anchor ui on to its outward neighbours, in proportion to the read counts of the links */
static void pass_on(uint32_t ui, float m, int dir, int32_t lo, int32_t hi) {
    NodeRef nr; if (!map_node(ws[ui].h, &nr)) return;
    const int32_t cu = ws[ui].coord; const uint8_t out = ((dir > 0) == (ws[ui].rev == 0)) ? 1 : 0;
    double tot = 0; uint8_t want_mate = 0;
    for (uint32_t x = 0; x < nr.deg; x++) {
        const EdgeD *E = &nr.e[x];
        if (E->mate || E->su != out || E->count < MINEDGE) continue;
        int32_t cv = cu + dir * E->dist;
        if (cv > hi || cv < lo) continue;
        tot += E->count;
    }
    if (tot <= 0) {
        want_mate = 1;
        for (uint32_t x = 0; x < nr.deg; x++) {
            const EdgeD *E = &nr.e[x];
            if (!E->mate || E->su != out || E->count < MINEDGE) continue;
            int32_t cv = cu + dir * E->dist;
            if (cv > hi || cv < lo) continue;
            tot += E->count;
        }
    }
    if (tot <= 0) return;
    for (uint32_t x = 0; x < nr.deg; x++) {
        const EdgeD *E = &nr.e[x];
        if (E->mate != want_mate || E->su != out || E->count < MINEDGE) continue;
        if (EPS > 0 && E->count == 1 && tot >= 8 && is_chimeric_hub(E->to, E->sv)) continue;
        int32_t cv = cu + dir * E->dist;
        if (cv > hi || cv < lo) continue;
        float share = (float)(m * E->count / tot);
        uint32_t vi = ws_get(E->to); WS *v = &ws[vi];
        if (v->placed) {                                  /* placed before */
            if (v->seed) continue;                        /* gene anchors are sources, they take no share */
            int32_t dd = cv - v->coord; if (dd < 0) dd = -dd;
            if (dd > SAMEPOS) continue;                   /* the same anchor at another place of the flank: no share */
            v->mass += share;
            if (EPS <= 0 || v->redo >= 2) continue;       /* at most twice per anchor: tandem repeats (true cycles) must not keep shares circulating */
            if (v->state == 1 && share >= EPS * Tlive) {              /* v was followed before this share arrived (the longer
                                                                         arm of an indel bubble closes late): pass the extra on */
                v->redo++; Tlive += share; hpush((HE){ dir * v->coord, vi, v->coord, ui, share, v->rev, 1 });
            } else if (v->state == 2 && v->mass >= EPS * Tlive) {     /* not followed before, now it has enough */
                v->redo++; v->state = 3; Tlive += v->mass; hpush((HE){ dir * v->coord, vi, v->coord, ui, 0, v->rev, 2 });
            }
            continue;
        }
        if (!v->hasm) { v->hasm = 1; v->mass = 0; v->state = 0; v->redo = 0; }
        v->mass += share; Tlive += share;
        uint8_t rv = (dir > 0) ? (E->sv == 1) : (E->sv == 0);
        hpush((HE){ dir * cv, vi, cv, ui, 0, rv, 0 });
    }
}

typedef struct { int seeds, capped; uint32_t reached; int32_t far; } WalkRes;
static WalkRes walk(const char *gname, int dir, Anc *ga, int gn, int glen) {
    ws_reset(); heap_n = 0; nball = 0; ndrop = 0; Tlive = 0;
    int seeds = 0, capped = 0; int32_t far = 0;
    for (int j = 0; j < gn; j++) {
        NodeRef nr;
        if (!map_node(ga[j].h, &nr) || ws_find(ga[j].h) != NONE) continue;
        uint32_t u = ws_get(ga[j].h);
        ws[u].placed = ws[u].hasm = 1; ws[u].coord = (int32_t)ga[j].pos; ws[u].rev = ga[j].st; ws[u].seed = 1; ws[u].mass = 1.0f;
        hpush((HE){ dir * (int32_t)ga[j].pos, u, (int32_t)ga[j].pos, NONE, 0, ga[j].st, 0 });
        seeds++;
    }
    int32_t hi = glen + D + SLACK, lo = -(D + SLACK);
    while (heap_n) {
        HE e = hpop(); uint32_t u = e.node; float m;
        if (e.kind == 1) {                                   /* share that arrived after u was followed */
            Tlive -= e.delta; if (Tlive < 0) Tlive = 0;
            pass_on(u, e.delta, dir, lo, hi); continue;
        }
        if (e.kind == 2) {
            if (ws[u].state != 3) continue;
            m = ws[u].mass; Tlive -= m; if (Tlive < 0) Tlive = 0;
        } else {
            if (ws[u].placed && !(ws[u].seed && e.pred == NONE)) continue;   /* already placed */
            if (!ws[u].placed) { ws[u].placed = 1; ws[u].coord = e.coord; ws[u].rev = e.rev; ws[u].seed = 0; }
            if (ws[u].state != 0) continue;
            m = ws[u].mass;
            if (!ws[u].seed) {
                /* A branch is followed if it carries at least EPS of the share that is still on its way at this distance.
                   Relative, because share is lost at every dead end (read errors): an absolute limit would starve the
                   true path of a noisy sample. At most about 1/EPS branches can be alive at one distance. */
                double before = Tlive;
                Tlive -= m; if (Tlive < 0) Tlive = 0;
                if (EPS > 0 && m < EPS * before) {
                    ws[u].state = 2;
                    if (ndrop < MAXBALL) { if (ndrop == drop_cap) { drop_cap = drop_cap ? drop_cap * 2 : 4096; drop = xrealloc(drop, (size_t)drop_cap * 4); } drop[ndrop++] = u; }
                    continue;
                }
            }
        }
        if (nball >= MAXBALL) { capped = 1; break; }
        ws[u].state = 1;
        if (nball == ball_cap) { ball_cap = ball_cap ? ball_cap * 2 : 4096; ball = xrealloc(ball, (size_t)ball_cap * 4); }
        ball[nball++] = u;
        { int32_t o = dir > 0 ? ws[u].coord - glen : -ws[u].coord; if (o > far) far = o; }
        pass_on(u, m, dir, lo, hi);
    }
    for (uint32_t i = 0; i < nball; i++) bset_add(ws[ball[i]].h);
    uint32_t nd = 0;
    for (uint32_t i = 0; i < ndrop; i++) nd += ws[drop[i]].state == 2;       /* branches not followed */
    fprintf(stderr, "walk %s %c: gene anchors %d, found %d, reached %u anchors, %d bp out, %u branches below share%s\n", gname,
            dir > 0 ? 'R' : 'L', gn, seeds, nball, far, nd, capped ? " (CAPPED)" : "");
    return (WalkRes){ seeds, capped, nball, far };
}

/* ---------- pass 2 ---------- */
#define LHMAX 1000
static uint64_t st_pairs = 0, st_rec = 0, lenhist[LHMAX + 1];
static void put_read(FILE *ot, uint64_t pair, int mate, const char *s, int l, const char *q, int ql) {
    /* pair, mate, read, an empty field (the place of the anchor list of earlier versions), base qualities */
    fprintf(ot, "%lu\t%d\t", (unsigned long)pair, mate);
    fwrite(s, 1, (size_t)l, ot); fputc('\t', ot); fputc('\t', ot);
    if (ql >= l) fwrite(q, 1, (size_t)l, ot);                     /* else: no usable qualities, left empty */
    fputc('\n', ot);
}
static void put_fq(FILE *o, char *ln[4], int len[4]) {
    for (int i = 0; i < 4; i++) { fwrite(ln[i], 1, (size_t)len[i], o); fputc('\n', o); }
}
/* one pair in pass 2: keep it if a reached anchor is in it; -> 1 if kept */
static int keep_pair(Scr *sc, Anc *a1, Anc *a2, uint64_t number, char *l1[4], int n1[4], char *l2[4], int n2[4],
                     FILE *ot, FILE *o1, FILE *o2, uint32_t *lh) {
    int m1 = minimizers(sc, l1[1], n1[1], a1, MAXA), m2 = minimizers(sc, l2[1], n2[1], a2, MAXA), hit = 0;
    for (int j = 0; j < m1 && !hit; j++) hit = bset_has(a1[j].h);
    for (int j = 0; j < m2 && !hit; j++) hit = bset_has(a2[j].h);
    if (!hit) return 0;
    if (o1) { put_fq(o1, l1, n1); put_fq(o2, l2, n2); }
    put_read(ot, number, 1, l1[1], n1[1], l1[3], n1[3]);
    put_read(ot, number, 2, l2[1], n2[1], l2[3], n2[3]);
    lh[n1[1] < LHMAX ? n1[1] : LHMAX]++; lh[n2[1] < LHMAX ? n2[1] : LHMAX]++;
    return 1;
}
static FILE *p2_ot, *p2_o1, *p2_o2;
static void batch_keep(Batch *b, Scr *sc, Anc *t1, Anc *t2) {
    /* the output text of a batch: with threads into memory (written later, in order), else straight to the files */
    FILE *f[3] = { p2_ot, p2_o1, p2_o2 };
    if (!b->lh) b->lh = xcalloc(LHMAX + 1, 4); else memset(b->lh, 0, (LHMAX + 1) * 4);
    if (NTHREADS > 1) {
        f[0] = open_memstream(&b->ob[0], &b->ol[0]);
        if (write_fq) { f[1] = open_memstream(&b->ob[1], &b->ol[1]); f[2] = open_memstream(&b->ob[2], &b->ol[2]); }
        if (!f[0] || (write_fq && (!f[1] || !f[2]))) die_oom(0);
    }
    for (int p = 0; p < b->n; p++) {
        char *l1[4], *l2[4];
        for (int k = 0; k < 4; k++) { l1[k] = b->txt + b->off[2 * p][k]; l2[k] = b->txt + b->off[2 * p + 1][k]; }
        b->rec += (uint64_t)keep_pair(sc, t1, t2, b->first + (uint64_t)p + 1, l1, b->len[2 * p], l2, b->len[2 * p + 1], f[0], f[1], f[2], b->lh);
    }
    if (NTHREADS > 1) for (int k = 0; k < 3; k++) if (f[k] && fclose(f[k])) die_oom(0);
}
static void apply_keep(Batch *b) {
    st_pairs += (uint64_t)b->n; st_rec += b->rec;
    for (int i = 0; i <= LHMAX; i++) lenhist[i] += b->lh[i];
    if (NTHREADS < 2) return;
    fwrite(b->ob[0], 1, b->ol[0], p2_ot); free(b->ob[0]); b->ob[0] = NULL;
    if (write_fq) {
        fwrite(b->ob[1], 1, b->ol[1], p2_o1); fwrite(b->ob[2], 1, b->ol[2], p2_o2);
        free(b->ob[1]); free(b->ob[2]); b->ob[1] = b->ob[2] = NULL;
    }
}
static void pass2(const char *fn1, const char *fn2, const char *out) {
    char p[4096];
    LR r1, r2; lr_open(&r1, fn1); lr_open(&r2, fn2);
    FILE *o1 = NULL, *o2 = NULL;
    if (write_fq) {
        snprintf(p, sizeof p, "%s_1.fq", out); o1 = xfopen(p, "w");
        snprintf(p, sizeof p, "%s_2.fq", out); o2 = xfopen(p, "w");
    }
    snprintf(p, sizeof p, "%s.reads.tsv", out); FILE *ot = xfopen(p, "w");
    int warned = have_pass1;                                 /* pass 1 already warned */
    p2_ot = ot; p2_o1 = o1; p2_o2 = o2;
    run_batches(&r1, &r2, &warned, batch_keep, apply_keep);
    lr_close(&r1); lr_close(&r2);
    if (o1) { xfclose(o1, "the reads near the genes"); xfclose(o2, "the reads near the genes"); }
    xfclose(ot, "the reads near the genes");
}

/* ---------- FASTA (genes) ---------- */
typedef void (*seqfn)(const char *name, const char *seq, int len, void *ctx);
static int each_fasta(const char *fn, seqfn f, void *ctx) {
    FILE *fp = xfopen(fn, "r");
    size_t cap = 1 << 20, n = 0, lc = 0; char *sq = xmalloc(cap), nm[1024] = "", *ln = NULL; int cnt = 0; ssize_t l;
    for (;;) {
        l = getline(&ln, &lc, fp);
        if (l < 0 || ln[0] == '>') {
            if (n) { f(nm, sq, (int)n, ctx); cnt++; }
            if (l < 0) break;
            nm[0] = 0; sscanf(ln + 1, "%1023s", nm); n = 0; continue;
        }
        if (n + (size_t)l + 1 > cap) { while (n + (size_t)l + 1 > cap) cap *= 2; sq = xrealloc(sq, cap); }
        for (ssize_t i = 0; i < l; i++)                     /* blanks, tabs and line ends are not sequence */
            if (ln[i] != ' ' && ln[i] != '\t' && ln[i] != '\n' && ln[i] != '\r' && ln[i] != '\v' && ln[i] != '\f') sq[n++] = ln[i];
    }
    free(sq); free(ln); fclose(fp);
    return cnt;
}
static Anc *GA; static int ga_cap = 0;
static int gene_anchors(const char *seq, int len) {
    if (len > ga_cap) { ga_cap = len * 2; GA = xrealloc(GA, sizeof(Anc) * (size_t)ga_cap); }
    return minimizers(&scr0, seq, len, GA, ga_cap);
}
static void dump_seq(const char *name, const char *seq, int len, void *ctx) {
    (void)ctx; int m = gene_anchors(seq, len);
    for (int j = 0; j < m; j++) printf("%s\t%u\t%016lx\t%d\n", name, GA[j].pos, (unsigned long)GA[j].h, GA[j].st);
}
static void seed_bloom_seq(const char *name, const char *seq, int len, void *ctx) {
    (void)name; (void)ctx; int m = gene_anchors(seq, len);
    for (int j = 0; j < m; j++) bloom_set_on(bloom, GA[j].h);
}
static FILE *f_genes;
static void walk_gene(const char *name, const char *seq, int len, void *ctx) {
    (void)ctx; int m = gene_anchors(seq, len);
    WalkRes r = walk(name, +1, GA, m, len), l = walk(name, -1, GA, m, len);
    fprintf(f_genes, "%s\t%d\t%d\t%d\t%u\t%u\t%d\t%d\t%d\n", name, len, m, r.seeds, r.reached, l.reached, r.far, l.far, r.capped + 2 * l.capped);
}

int main(int argc, char **argv) {
    const char *genes = NULL, *out = "aw_out", *dump = NULL, *map_out = NULL, *map_in = NULL; int bloom_log2 = 28, c;
    uint32_t force_parts = 0; double mem_mb = 0;
    while ((c = getopt(argc, argv, "k:w:D:b:e:Bg:o:A:vp:M:S:G:ft:Lm:P:i:")) != -1) {
        if (c == 'v') { printf("aw %s (anchor walk)\n", AW_VERSION); return 0; }
        if (c == 'f') write_fq = 1;
        if (c == 't') NTHREADS = atoi(optarg);
        if (c == 'L') LOWC = 0;
        if (c == 'i') MATE_INS = atoi(optarg);
        if (c == 'p') EPS = atof(optarg);
        if (c == 'M') MAXBALL = (uint32_t)atol(optarg);
        if (c == 'S') map_out = optarg;
        if (c == 'G') map_in = optarg;
        if (c == 'm') mem_mb = atof(optarg);
        if (c == 'P') force_parts = (uint32_t)atol(optarg);
        if (c == 'k') K = atoi(optarg); else if (c == 'w') W = atoi(optarg); else if (c == 'D') D = atoi(optarg);
        else if (c == 'b') bloom_log2 = atoi(optarg); else if (c == 'e') MINEDGE = atoi(optarg); else if (c == 'B') use_bloom = 0;
        else if (c == 'g') genes = optarg; else if (c == 'o') out = optarg; else if (c == 'A') dump = optarg;
        else if (c == '?') { argc = 0; break; }
    }
    if (NTHREADS < 1) NTHREADS = 1;
    if (NTHREADS > 64) NTHREADS = 64;
    signal(SIGPIPE, SIG_IGN);
    if (argc == 0 || K > 32 || K < 11 || W < 1 || D < 1 || D > 100000000 || bloom_log2 < 12 || bloom_log2 > 36 || EPS < 0 || EPS >= 1 || MAXBALL < 1000 ||
        MAXBALL > 500000000u || mem_mb < 0 || force_parts > 4096 || (!dump && (!genes || argc - optind != 2))) {
        fprintf(stderr, "aw %s\nusage: aw [-k 31 -w 20 -D 10000 -b 28 -e 1 -p 0.002 -M 2000000 -t 1 -B -m MB -P parts -S map -G map -f] -g genes.fa -o out r1.fq[.gz] r2.fq[.gz]\n"
                        "       aw [-k 31 -w 20] -A seqs.fa\n", AW_VERSION);
        return 1;
    }
    KMASK = K < 32 ? (1ULL << (2 * K)) - 1 : ~0ULL;
    memset(code, -1, sizeof code);
    code['A'] = code['a'] = 0; code['C'] = code['c'] = 1; code['G'] = code['g'] = 2; code['T'] = code['t'] = 3;
    if (dump) { each_fasta(dump, dump_seq, NULL); return 0; }

    char p[4096], tmpmap[4096];
    snprintf(p, sizeof p, "%s.genes.tsv", out); f_genes = xfopen(p, "w");  /* fail early if out/ is not writable */
    { FILE *g = xfopen(genes, "r"); fclose(g); }
    if (map_out) { FILE *g = xfopen(map_out, "wb"); fclose(g); }       /* fail before the read pass, not after it */
    snprintf(tmpmap, sizeof tmpmap, "%s.map", out);

    bloom_mask = (1ULL << (bloom_log2 - 6)) - 1;
    if (map_in) use_bloom = 0;
    if (use_bloom) bloom = big_alloc((bloom_mask + 1) * 8);
    MEM_BUDGET = (uint64_t)(mem_mb * 1048576.0);
    const uint64_t fixed = use_bloom ? (bloom_mask + 1) * 8 : 0;
    if (!map_in && MEM_BUDGET && MEM_BUDGET < fixed + (1ULL << 20)) {
        fprintf(stderr, "aw: -m %.0f MB is too little: the Bloom filter alone takes %.0f MB (-b %d), and the map needs room beside it\n",
                mem_mb, fixed / 1048576.0, bloom_log2);
        return 2;
    }

    double t0 = now(), bfp = 0;
    uint32_t nparts = 1; int passes1 = 0;
    if (map_in) map_open(map_in, 0);
    else {
        if (use_bloom) {
            bloom1 = big_alloc((bloom_mask + 1) * 8);
            each_fasta(genes, seed_bloom_seq, NULL);
            pass0(argv[optind], argv[optind + 1]);
            bfp = bloom_false_rate(bloom1);
            big_free(bloom1, (bloom_mask + 1) * 8); bloom1 = NULL;
        }
        nparts = force_parts ? force_parts : 1;
        for (;;) {
            const int to_file = nparts > 1 || map_out != NULL;
            const char *mfn = map_out ? map_out : tmpmap;
            if (to_file) map_begin(mfn, nparts);
            NPART = nparts; p1_over = 0; tot_nodes = tot_edges = 0; st_reads = st_anchors = st_solid = 0;
            for (PART = 0; PART < nparts; PART++) {
                map_free(); idx_rebuild(1 << 16);
                pass1(argv[optind], argv[optind + 1]); passes1++;
                if (p1_over) break;
                tot_nodes += n_nodes; tot_edges += n_edges;
                if (nparts > 1) fprintf(stderr, "aw: part %u of %u of the map: %u anchors, %u links, %.0f MB\n", PART + 1, nparts, n_nodes, n_edges / 2, map_bytes() / 1048576.0);
                if (to_file) { map_part_save(PART); map_free(); }
            }
            if (!p1_over) {
                if (to_file) { map_end(nparts); map_open(mfn, map_out == NULL); }
                break;
            }
            /* the part did not fit: more parts */
            if (to_file) { fclose(mapf); mapf = NULL; free(mdir); mdir = NULL; if (!map_out) unlink(mfn); }
            uint32_t np;
            if (p1_over == 2) np = nparts * 2;
            else {
                /* What this part would have taken at the end of the files. A map grows with the reads by a power between 1
                   (every read brings its share) and 2 (anchors seen twice are few at first); which, is read off the two
                   moments at which the part was half full and full. Later the growth can only slow down. */
                double used = (double)(map_bytes() - fixed), g = 1.5;
                if (p1_half > 0 && p1_frac > 1.02 * p1_half) { g = log(2.0) / log(p1_frac / p1_half); if (g < 1) g = 1; if (g > 2) g = 2; }
                double est = 1.1 * used / pow(p1_frac, g);
                np = (uint32_t)ceil((double)nparts * est / (0.8 * (double)(MEM_BUDGET - fixed)));
                if (np <= nparts) np = nparts + 1;
            }
            if (force_parts || np > 4096) {
                fprintf(stderr, "aw: the map of this sample does not fit in %.0f MB of memory%s. Give more memory, or use a larger -w (fewer anchors).\n",
                        mem_mb, force_parts ? " with this number of parts (-P)" : ", not even in 4096 parts");
                return 2;
            }
            if (p1_over == 2) fprintf(stderr, "aw: the map of this sample has more anchors or links than one table holds; building it in %u parts, each with one more pass over the reads\n", np);
            else fprintf(stderr, "aw: the map of this sample does not fit in %.0f MB of memory (%s after %.0f %% of the reads); building it in %u parts, each with one more pass over the reads\n",
                         mem_mb, nparts > 1 ? "a part was full" : "full", 100 * p1_frac, np);
            map_free(); nparts = np;
        }
    }
    double t1 = now();
    long rss1 = rss_peak_mb();
    int b_suggest = bloom_log2;
    if (use_bloom && bfp > 0.05) {
        /* distinct anchors <= first sightings that missed the filter + anchors stored as nodes;
           about 9.8 anchors per 64-bit filter word keeps the false rate near 5% */
        double distinct = (double)(st_anchors - st_solid) + (double)tot_nodes;
        b_suggest = 12;
        while (b_suggest < 36 && 9.8 * (double)(1ULL << (b_suggest - 6)) < distinct) b_suggest++;
        fprintf(stderr, "aw: warning: Bloom filter too small for this sample: about %.0f%% of one-off anchors were kept as if seen twice.\n"
                        "    This costs memory and time. Next time use -b %d (%.0f MB).\n", 100 * bfp, b_suggest,
                (double)(1ULL << b_suggest) / 8 / 1048576);
    }
    big_free(bloom, (bloom_mask + 1) * 8); bloom = NULL;

    fprintf(f_genes, "gene\tlength\tanchors\tfound\treached_right\treached_left\tbp_right\tbp_left\tcapped\n");
    int ngenes = each_fasta(genes, walk_gene, NULL);
    xfclose(f_genes, "the walk output");
    if (!bset) { bset_mask = 1023; bset = xcalloc(bset_mask + 1, 8); }
    const uint64_t graph_bytes = tot_nodes * sizeof(Node) + tot_edges * sizeof(Edge) + (mm_base ? 0 : idx_size * 4);
    const int from_file = mm_base != NULL;
    map_close();                                         /* pass 2 needs the reached anchors only */
    free(ws); free(wsi); free(ball); free(drop); free(heap); free(ebuf);
    double t2 = now();

    have_pass1 = map_in ? 0 : 1;
    quiet_eof = have_pass1;
    pass2(argv[optind], argv[optind + 1], out);
    double t3 = now();
    if (st_pairs * 2 != st_reads) {
        if (map_in) fprintf(stderr, "aw: the map %s was made from %lu read pairs, but these files hold %lu: not the same reads\n",
                            map_in, (unsigned long)(st_reads / 2), (unsigned long)st_pairs);
        else fprintf(stderr, "aw: the read files gave %lu pairs in the first pass and %lu in the second. aw reads them twice: "
                             "give files, not pipes\n", (unsigned long)(st_reads / 2), (unsigned long)st_pairs);
        exit(1);
    }
    uint64_t kept = 0, acc = 0, sum = 0; int med = 0;
    for (int i = 0; i <= LHMAX; i++) { kept += lenhist[i]; sum += lenhist[i] * (uint64_t)i; }
    for (int i = 0; i <= LHMAX; i++) { acc += lenhist[i]; if (acc * 2 >= kept && kept) { med = i; break; } }

    snprintf(p, sizeof p, "%s.stats.txt", out); FILE *fs = xfopen(p, "w");
    for (int i = 0; i < 2; i++) {
        FILE *o = i ? fs : stderr;
        fprintf(o, "aw %s k=%d w=%d D=%d bloom=%s b=%d min_edge=%d min_share=%g genes=%d\n", AW_VERSION, K, W, D, use_bloom ? "on" : "off",
                bloom_log2, MINEDGE, EPS, ngenes);
        fprintf(o, "reads=%lu anchors_seen=%lu (%.2f/read) solid_occurrences=%lu bloom_false_rate=%.4f\n", (unsigned long)st_reads,
                (unsigned long)st_anchors, (double)st_anchors / (double)(st_reads ? st_reads : 1), (unsigned long)st_solid, bfp);
        fprintf(o, "graph_nodes=%lu graph_links=%lu map_parts=%u map_from_file=%d read_passes=%d\n", (unsigned long)tot_nodes, (unsigned long)(tot_edges / 2),
                map_in ? mm_parts : nparts, from_file, passes1 + 1);
        fprintf(o, "reached_anchors=%lu recruited_pairs=%lu of %lu\n", (unsigned long)n_reached, (unsigned long)st_rec, (unsigned long)st_pairs);
        fprintf(o, "kept_reads=%lu read_len_median=%d read_len_mean=%.1f\n", (unsigned long)kept, med, kept ? (double)sum / (double)kept : 0.0);
        fprintf(o, "pass1_s=%.2f walk_s=%.3f pass2_s=%.2f total_s=%.2f\n", t1 - t0, t2 - t1, t3 - t2, t3 - t0);
        fprintf(o, "peak_rss_mb_after_pass1=%ld peak_rss_mb=%ld graph_bytes=%lu bloom_bytes=%lu\n", rss1, rss_peak_mb(),
                (unsigned long)graph_bytes, (unsigned long)fixed);
    }
    xfclose(fs, "the statistics");
    return 0;
}
