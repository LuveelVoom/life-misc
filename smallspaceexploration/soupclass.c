/* soupclass.c — soup-based explosiveness statistics for range-1 Moore rules.
 *
 * Input (stdin), one rule per line:  <label> <128 hex digits>
 *   512-entry transition table, MAP index order (row-major from NW; NW=256 N=128 NE=64
 *   W=32 C=16 E=8 SW=4 S=2 SE=1), entry 0 first, 4 entries per hex digit, lowest index
 *   in the digit's MSB.  B0 rules (entry 0 set) are not supported (rejected).
 * Output (stdout), one line per rule: <label> key=value ...
 *
 * Each soup: w x w random cells at d% density, centred on an N x N bounded plane
 * (cells outside are permanently dead).  Soup i is identical for every rule.
 * Outcomes per soup:
 *   explode : occupied 16x16 tiles reached T before maxgen, or (checked every 32 gens)
 *             >= T/4 tiles with occ(t) >= 2.8 occ(t/2) and occ(t/2) >= 2.8 occ(t/4)
 *   settled : state became periodic (period <= 1024, incl. death)
 *   grow    : neither, and occupied tiles at maxgen >= 1.5 x (at maxgen/2) and +8
 *   active  : anything else (long-lived, or settled core with escaping ships)
 * Build: gcc -O3 -march=native -fopenmp soupclass.c -o soupclass -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define TS 16
#define PMAX 1024

enum { R_EXPLODE, R_SETTLED, R_GROW, R_ACTIVE };

typedef struct {
    int N, NT, S;               /* grid side, tiles per side, row stride (N+2) */
    uint8_t *buf[2];
    uint16_t *tpop;
    int32_t *stamp;
    int *list[2], nlist[2];
    int *edge;                  /* per active tile: edge-change bits */
    uint64_t *hist;
    int32_t *chg;               /* changed cells per generation */
    int32_t *occh;              /* occupied tiles per generation */
} World;

typedef struct {
    int outcome, gen, died;
    double speed;               /* explode only */
    int lifespan;               /* non-explode */
    int occ_half, occ_end;      /* occupied tiles at maxgen/2 and at end */
    int trend;                  /* explosion detected by growth trend, not threshold */
} SoupRes;

static int N_ = 1024, MAXGEN = 6000, NSOUP = 64, DENS = 50, SW = 16, TEXP = 512;
static uint64_t SEED = 12345;
static int DUMP = 0;   /* -D: run soup 0 for maxgen gens, print live cells (testing) */
static int VERB = 0;   /* -v: per-soup lines on stderr */
static int DUMPIDX = 0; /* -k: soup index for -D */
static int ERASE = 1;   /* -E0: disable absorbing border (dead outer cells only) */

static inline uint64_t mix64(uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static World *world_new(int N) {
    World *w = calloc(1, sizeof *w);
    w->N = N; w->NT = N / TS; w->S = N + 2;
    size_t cells = (size_t)(N + 2) * (N + 2), nt = (size_t)w->NT * w->NT;
    w->buf[0] = calloc(cells, 1); w->buf[1] = calloc(cells, 1);
    w->tpop = calloc(nt, sizeof *w->tpop);
    w->stamp = calloc(nt, sizeof *w->stamp);
    w->list[0] = malloc(nt * sizeof(int)); w->list[1] = malloc(nt * sizeof(int));
    w->edge = malloc(nt * sizeof(int));
    w->hist = malloc((MAXGEN + 2) * sizeof *w->hist);
    w->chg = malloc((MAXGEN + 2) * sizeof *w->chg);
    w->occh = malloc((MAXGEN + 2) * sizeof *w->occh);
    return w;
}

static int cmpint(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }

/* lifespan: last generation whose changed-cell count exceeds max(20, 3 x median of the
 * final 64 generations' counts).  Approximate; ignores steady ship/oscillator activity. */
static int lifespan_of(const int32_t *chg, int g) {
    if (g <= 0) return 0;
    int k = g < 64 ? g : 64, tmp[64];
    for (int i = 0; i < k; i++) tmp[i] = chg[g - i];
    qsort(tmp, k, sizeof(int), cmpint);
    int thr = 3 * tmp[k / 2]; if (thr < 20) thr = 20;
    for (int i = g; i >= 1; i--) if (chg[i] > thr) return i;
    return 0;
}

static SoupRes run_soup(World *w, const uint8_t *tab, int soup) {
    const int N = w->N, NT = w->NT, S = w->S;
    SoupRes r = {0};
    r.outcome = R_ACTIVE; r.gen = MAXGEN;
    memset(w->buf[0], 0, (size_t)S * S); memset(w->buf[1], 0, (size_t)S * S);
    memset(w->tpop, 0, (size_t)NT * NT * sizeof *w->tpop);
    memset(w->stamp, 0xff, (size_t)NT * NT * sizeof *w->stamp);   /* -1 */
    /* soup (depends only on SEED, soup index) */
    uint64_t st = mix64(SEED * 0x100000001b3ULL + (uint64_t)soup), hash = 0;
    long pop = 0; int occ = 0;
    int x0 = N / 2 - SW / 2, y0 = N / 2 - SW / 2;
    for (int y = 0; y < SW; y++)
        for (int x = 0; x < SW; x++) {
            st = mix64(st);
            if ((int)(st % 100) < DENS) {
                int X = x0 + x, Y = y0 + y;
                w->buf[0][(size_t)(Y + 1) * S + X + 1] = w->buf[1][(size_t)(Y + 1) * S + X + 1] = 1;
                hash ^= mix64(((uint64_t)Y << 20) | (uint64_t)X);
                pop++;
                int t = (Y / TS) * NT + X / TS;
                if (w->tpop[t]++ == 0) occ++;
            }
        }
    /* schedule tiles overlapping the soup, plus neighbours */
    int cur = 0, nl = 0;
    for (int ty = (y0 / TS) - 1; ty <= (y0 + SW - 1) / TS + 1; ty++)
        for (int tx = (x0 / TS) - 1; tx <= (x0 + SW - 1) / TS + 1; tx++)
            if (tx >= 0 && ty >= 0 && tx < NT && ty < NT) {
                int t = ty * NT + tx; w->stamp[t] = 0; w->list[0][nl++] = t;
            }
    w->nlist[0] = nl;
    w->hist[0] = hash; w->chg[0] = 0; w->occh[0] = occ;
    int Gend = 0;
    for (int g = 0; g < MAXGEN; g++) {
        const uint8_t *src = w->buf[cur]; uint8_t *dst = w->buf[cur ^ 1];
        int *L = w->list[g & 1], nL = w->nlist[g & 1];
        int *NL = w->list[(g & 1) ^ 1], nN = 0;
        long chg = 0;
        for (int li = 0; li < nL; li++) {
            int t = L[li], tx = t % NT, ty = t / NT, ebits = 0, dpop = 0;
            int bx = tx * TS, by = ty * TS;
            /* absorbing border: outer ring of tiles is forced dead (erases escaping ships) */
            const uint8_t emask = (ERASE && (tx == 0 || ty == 0 || tx == NT - 1 || ty == NT - 1)) ? 0 : 1;
            for (int yy = 0; yy < TS; yy++) {
                int Y = by + yy;
                const uint8_t *a = src + (size_t)Y * S + bx, *c = a + S, *b = c + S;
                uint8_t *o = dst + (size_t)(Y + 1) * S + bx + 1;
                unsigned idx = (a[0] << 7) | (a[1] << 6) | (c[0] << 4) | (c[1] << 3) | (b[0] << 1) | b[1];
                for (int xx = 0; xx < TS; xx++) {
                    idx = ((idx << 1) & 0x1B6) | (a[xx + 2] << 6) | (c[xx + 2] << 3) | b[xx + 2];
                    uint8_t nv = tab[idx] & emask;
                    o[xx] = nv;
                    if (nv != c[xx + 1]) {
                        chg++; dpop += nv ? 1 : -1;
                        hash ^= mix64(((uint64_t)Y << 20) | (uint64_t)(bx + xx));
                        if (yy == 0) ebits |= 1;
                        if (yy == TS - 1) ebits |= 2;
                        if (xx == 0) ebits |= 4;
                        if (xx == TS - 1) ebits |= 8;
                        ebits |= 16;
                    }
                }
            }
            if (ebits) {
                int before = w->tpop[t];
                w->tpop[t] = (uint16_t)(before + dpop);
                pop += dpop;
                if (before == 0 && w->tpop[t]) occ++;
                else if (before && w->tpop[t] == 0) occ--;
                /* schedule self and neighbours touched by edge changes */
                int up = ebits & 1, dn = ebits & 2, lf = ebits & 4, rt = ebits & 8;
                for (int dy = -1; dy <= 1; dy++) {
                    if ((dy == -1 && !up) || (dy == 1 && !dn)) continue;
                    int ny = ty + dy; if (ny < 0 || ny >= NT) continue;
                    for (int dx = -1; dx <= 1; dx++) {
                        if ((dx == -1 && !lf) || (dx == 1 && !rt)) continue;
                        int nx = tx + dx; if (nx < 0 || nx >= NT) continue;
                        int nt = ny * NT + nx;
                        if (w->stamp[nt] != g + 1) { w->stamp[nt] = g + 1; NL[nN++] = nt; }
                    }
                }
            }
        }
        w->nlist[(g & 1) ^ 1] = nN;
        cur ^= 1;
        int G = g + 1;
        w->hist[G] = hash; w->chg[G] = (int32_t)chg; w->occh[G] = occ;
        if (DUMP) continue;
        Gend = G;
        if (pop == 0) { r.outcome = R_SETTLED; r.gen = G; r.died = 1; break; }
        /* explosion: occupied tiles >= TEXP, or (checked every 32 gens) >= TEXP/4 with
         * near-quadratic growth over the last two doublings of time (ratio >= 2.8 each;
         * linear growth gives ~2).  speed: area-equivalent half-width growth over [G/2, G]. */
        int trend = 0;
        if (G >= 128 && (G & 31) == 0 && occ >= TEXP / 4) {
            int oh = w->occh[G / 2], oq = w->occh[G / 4];
            trend = occ >= 2.8 * oh && oh >= 2.8 * oq;
        }
        if (occ >= TEXP || trend) {
            r.outcome = R_EXPLODE; r.gen = G; r.trend = trend && occ < TEXP;
            r.speed = 8.0 * (sqrt((double)occ) - sqrt((double)w->occh[G / 2])) / (G - G / 2);
            return r;
        }
        if (nN == 0) { r.outcome = R_SETTLED; r.gen = G - 1; break; }  /* still life: state G == G-1 */
        /* periodicity: hash(G) == hash(G-p) and hash(G-1) == hash(G-1-p) */
        int pm = G - 1 < PMAX ? G - 1 : PMAX, found = 0;
        for (int p = 2; p <= pm; p++)
            if (w->hist[G - p] == hash && w->hist[G - 1 - p] == w->hist[G - 1]) { found = p; break; }
        if (found) {
            /* first gen of the cycle */
            int s = G - found;
            while (s - 1 >= 0 && w->hist[s - 1] == w->hist[s - 1 + found]) s--;
            r.outcome = R_SETTLED; r.gen = s; break;
        }
    }
    if (DUMP) {
        const uint8_t *b = w->buf[cur];
        for (int y = 0; y < N; y++) for (int x = 0; x < N; x++)
            if (b[(size_t)(y + 1) * S + x + 1]) printf("%d %d\n", x, y);
        printf("end pop=%ld\n", pop);
        return r;
    }
    if (r.outcome == R_SETTLED) {
        /* r.gen = first gen of the final cycle; includes time for ships to reach the border,
         * so lifespan uses the activity heuristic, capped by r.gen */
        int ls = lifespan_of(w->chg, Gend);
        r.lifespan = ls < r.gen ? ls : r.gen;
        return r;
    }
    int oe = w->occh[MAXGEN], oh = w->occh[MAXGEN / 2];
    r.occ_half = oh; r.occ_end = oe;
    if (oe >= 1.5 * oh && oe - oh >= 8) r.outcome = R_GROW;
    r.lifespan = lifespan_of(w->chg, MAXGEN);
    return r;
}

static int cmpd(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y);
}
static void pmed(const char *k, double *v, int n, const char *f) {
    if (!n) { printf(" %s=na", k); return; }
    qsort(v, n, sizeof *v, cmpd);
    double m = n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    printf(" %s=", k); printf(f, m);
}

int main(int argc, char **argv) {
    int opt, threads = 0;
    while ((opt = getopt(argc, argv, "n:g:s:d:w:N:T:t:Dvk:E:")) != -1) {
        switch (opt) {
        case 'n': NSOUP = atoi(optarg); break;
        case 'g': MAXGEN = atoi(optarg); break;
        case 's': SEED = strtoull(optarg, 0, 10); break;
        case 'd': DENS = atoi(optarg); break;
        case 'w': SW = atoi(optarg); break;
        case 'N': N_ = atoi(optarg); break;
        case 'T': TEXP = atoi(optarg); break;
        case 't': threads = atoi(optarg); break;
        case 'D': DUMP = 1; NSOUP = 1; break;
        case 'v': VERB = 1; break;
        case 'k': DUMPIDX = atoi(optarg); break;
        case 'E': ERASE = atoi(optarg); break;
        default: fprintf(stderr, "bad option\n"); return 2;
        }
    }
    if (N_ % TS || SW > N_ / 2 || TEXP < 4) { fprintf(stderr, "bad sizes\n"); return 2; }
#ifdef _OPENMP
    if (threads > 0) omp_set_num_threads(threads);
    int nth = omp_get_max_threads();
#else
    int nth = 1;
#endif
    World **ws = malloc(nth * sizeof *ws);
    for (int i = 0; i < nth; i++) ws[i] = world_new(N_);
    SoupRes *res = malloc(NSOUP * sizeof *res);
    char label[256], hex[1024];
    while (scanf("%255s %1023s", label, hex) == 2) {
        if (strlen(hex) != 128) { fprintf(stderr, "%s: need 128 hex digits\n", label); return 2; }
        uint8_t tab[512];
        for (int i = 0; i < 128; i++) {
            char ch = hex[i]; int v = ch <= '9' ? ch - '0' : (ch | 32) - 'a' + 10;
            for (int k = 0; k < 4; k++) tab[4 * i + k] = (v >> (3 - k)) & 1;
        }
        if (tab[0]) { fprintf(stderr, "%s: B0 not supported\n", label); return 2; }
        if (DUMP) { run_soup(ws[0], tab, DUMPIDX); continue; }
        #pragma omp parallel for schedule(dynamic, 1)
        for (int s = 0; s < NSOUP; s++) {
#ifdef _OPENMP
            World *w = ws[omp_get_thread_num()];
#else
            World *w = ws[0];
#endif
            res[s] = run_soup(w, tab, s);
        }
        int c[4] = {0}, died = 0, ne = 0, nl = 0, nsx = 0, ntr = 0;
        double *tx = malloc(NSOUP * sizeof(double)), *vv = malloc(NSOUP * sizeof(double)),
               *lf = malloc(NSOUP * sizeof(double)), *tsv = malloc(NSOUP * sizeof(double));
        int nv = 0;
        for (int s = 0; s < NSOUP; s++) {
            c[res[s].outcome]++; died += res[s].died; ntr += res[s].trend;
            if (res[s].outcome == R_EXPLODE) {
                tx[ne++] = res[s].gen;
                if (res[s].speed > 0) vv[nv++] = res[s].speed;
            } else lf[nl++] = res[s].lifespan;
            if (res[s].outcome == R_SETTLED) tsv[nsx++] = res[s].gen;
            if (VERB) fprintf(stderr, "%s soup=%d out=%d gen=%d life=%d occ_half=%d occ_end=%d v=%.3f\n", label, s,
                              res[s].outcome, res[s].gen, res[s].lifespan, res[s].occ_half, res[s].occ_end, res[s].speed);
        }
        printf("%s n=%d explode=%d xtrend=%d settled=%d died=%d grow=%d active=%d",
               label, NSOUP, c[R_EXPLODE], ntr, c[R_SETTLED], died, c[R_GROW], c[R_ACTIVE]);
        pmed("tx", tx, ne, "%.0f"); pmed("v", vv, nv, "%.4f");
        pmed("life", lf, nl, "%.0f"); pmed("ts", tsv, nsx, "%.0f");
        printf("\n"); fflush(stdout);
        free(tx); free(vv); free(lf); free(tsv);
    }
    return 0;
}
