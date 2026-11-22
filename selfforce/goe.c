// Garden of Eden search: beam search over output rows of width W.
// State = set of predecessor row pairs (r_{i}, r_{i+1}), each W+2 bits, consistent
// with the output rows chosen so far. Output row i is determined by predecessor rows
// i, i+1, i+2 (rows offset by one: predecessor covers the output plus a 1-cell border).
// Empty set => the chosen W x H pattern has no predecessor.
// Rule: 512 chars of 0/1, MAP index order (NW 256 ... SE 1). Bit x of a row = column x.
// The beam keeps the BEAM smallest distinct sets; SEED only breaks ties.
// usage: goe TABLE W BEAM MAXH SEED
//   TABLE from: python3 -c "import sfcore as S; print(''.join(map(str, S.rule_table('B3/S23-a5'))))"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int W, NB, NP;              // NB = W+2 bits per pred row, NP = 2^(2NB) pairs
static unsigned char T[512];
static int WORDS;                  // 64-bit words per state

typedef struct { uint64_t *bits; int *rows; int h; long pop; int par; } St;

static uint64_t rng = 88172645463325252ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static inline int outrow(int a, int b, int c) {
    int o = 0;
    for (int x = 0; x < W; x++) {
        int idx = (((a >> x) & 1) << 8) | (((a >> (x + 1)) & 1) << 7) | (((a >> (x + 2)) & 1) << 6)
                | (((b >> x) & 1) << 5) | (((b >> (x + 1)) & 1) << 4) | (((b >> (x + 2)) & 1) << 3)
                | (((c >> x) & 1) << 2) | (((c >> (x + 1)) & 1) << 1) | ((c >> (x + 2)) & 1);
        o |= T[idx] << x;
    }
    return o;
}

static unsigned char *OUT; // OUT[(a<<2NB)|(b<<NB)|c] if small enough, else NULL

static long popc(uint64_t *b) { long s = 0; for (int i = 0; i < WORDS; i++) s += __builtin_popcountll(b[i]); return s; }

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: goe TABLE W BEAM MAXH SEED\n"); return 1; }
    if (strlen(argv[1]) != 512) { fprintf(stderr, "table must be 512 chars\n"); return 1; }
    for (int i = 0; i < 512; i++) T[i] = argv[1][i] == '1';
    W = atoi(argv[2]); int BEAM = atoi(argv[3]); int MAXH = atoi(argv[4]);
    rng ^= (uint64_t)atoll(argv[5]) * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) xr();
    NB = W + 2; NP = 1 << (2 * NB); WORDS = (NP + 63) / 64;
    int NR = 1 << NB, NO = 1 << W;
    if ((long)NR * NR * NR <= (1L << 27)) {
        OUT = malloc((size_t)NR * NR * NR);
        for (int a = 0; a < NR; a++) for (int b = 0; b < NR; b++) for (int c = 0; c < NR; c++)
            OUT[((size_t)a << (2 * NB)) | (b << NB) | c] = (unsigned char)outrow(a, b, c);
    }
    St *beam = calloc(BEAM, sizeof(St)); int nb = 1;
    beam[0].bits = malloc(WORDS * 8); memset(beam[0].bits, 0xff, WORDS * 8);
    beam[0].rows = calloc(MAXH, sizeof(int)); beam[0].h = 0; beam[0].pop = NP;
    uint64_t **succ = malloc(NO * sizeof(uint64_t *));
    for (int o = 0; o < NO; o++) succ[o] = malloc(WORDS * 8);
    int *cpar = malloc(sizeof(int) * BEAM * NO), *crow = malloc(sizeof(int) * BEAM * NO);
    long *cpop = malloc(sizeof(long) * BEAM * NO); uint64_t *ckey = malloc(8 * (size_t)BEAM * NO);
    for (int h = 1; h <= MAXH; h++) {
        long nc = 0;
        for (int s = 0; s < nb; s++) {
            for (int o = 0; o < NO; o++) memset(succ[o], 0, WORDS * 8);
            uint64_t *B = beam[s].bits;
            for (int wi = 0; wi < WORDS; wi++) {
                uint64_t w = B[wi];
                while (w) {
                    int p = wi * 64 + __builtin_ctzll(w); w &= w - 1;
                    int a = p >> NB, b = p & (NR - 1);
                    for (int c = 0; c < NR; c++) {
                        int o = OUT ? OUT[((size_t)a << (2 * NB)) | (b << NB) | c] : outrow(a, b, c);
                        int q = (b << NB) | c;
                        succ[o][q >> 6] |= 1ULL << (q & 63);
                    }
                }
            }
            for (int o = 0; o < NO; o++) {
                long pc = popc(succ[o]);
                if (pc == 0) {
                    printf("GOE W=%d H=%d rows:", W, h);
                    for (int i = 0; i < beam[s].h; i++) printf(" %d", beam[s].rows[i]);
                    printf(" %d\n", o);
                    for (int i = 0; i < beam[s].h; i++) { for (int x = 0; x < W; x++) putchar((beam[s].rows[i] >> x) & 1 ? 'o' : '.'); putchar('\n'); }
                    for (int x = 0; x < W; x++) putchar((o >> x) & 1 ? 'o' : '.'); putchar('\n');
                    fflush(stdout); return 0;
                }
                uint64_t k = 1469598103934665603ULL;
                for (int i = 0; i < WORDS; i++) { k ^= succ[o][i]; k *= 1099511628211ULL; }
                cpar[nc] = s; crow[nc] = o; cpop[nc] = pc; ckey[nc] = k; nc++;
            }
        }
        // select BEAM best distinct (by pop, random tiebreak)
        long *ord = malloc(sizeof(long) * nc); uint64_t *rk = malloc(8 * nc);
        for (long i = 0; i < nc; i++) { ord[i] = i; rk[i] = xr(); }
        // shell sort by (pop, random)
        for (long gap = nc / 2; gap > 0; gap /= 2)
            for (long i = gap; i < nc; i++) {
                long t = ord[i], j = i;
                while (j >= gap && (cpop[ord[j - gap]] > cpop[t] || (cpop[ord[j - gap]] == cpop[t] && rk[ord[j - gap]] > rk[t]))) { ord[j] = ord[j - gap]; j -= gap; }
                ord[j] = t;
            }
        St *nbm = calloc(BEAM, sizeof(St)); int nn = 0;
        uint64_t *keys = malloc(8 * BEAM);
        for (long ii = 0; ii < nc && nn < BEAM; ii++) {
            long i = ord[ii]; int dup = 0;
            for (int j = 0; j < nn; j++) if (keys[j] == ckey[i]) { dup = 1; break; }
            if (dup) continue;
            keys[nn] = ckey[i];
            St *p = &beam[cpar[i]], *q = &nbm[nn++];
            q->rows = malloc(MAXH * sizeof(int)); memcpy(q->rows, p->rows, p->h * sizeof(int));
            q->rows[p->h] = crow[i]; q->h = p->h + 1; q->pop = cpop[i]; q->bits = NULL; q->par = cpar[i];
        }
        // recompute bits of chosen states from their parent + row
        for (int j = 0; j < nn; j++) {
            int par = nbm[j].par; int o = nbm[j].rows[nbm[j].h - 1];
            uint64_t *nbits = calloc(WORDS, 8), *B = beam[par].bits;
            for (int wi = 0; wi < WORDS; wi++) {
                uint64_t w = B[wi];
                while (w) {
                    int p = wi * 64 + __builtin_ctzll(w); w &= w - 1;
                    int a = p >> NB, b = p & (NR - 1);
                    for (int c = 0; c < NR; c++) {
                        int oo = OUT ? OUT[((size_t)a << (2 * NB)) | (b << NB) | c] : outrow(a, b, c);
                        if (oo == o) { int q = (b << NB) | c; nbits[q >> 6] |= 1ULL << (q & 63); }
                    }
                }
            }
            nbm[j].bits = nbits;
        }
        for (int s = 0; s < nb; s++) { free(beam[s].bits); free(beam[s].rows); }
        free(beam); beam = nbm; nb = nn; free(ord); free(rk); free(keys);
        fprintf(stderr, "h=%d beam=%d best=%ld worst=%ld\n", h, nb, beam[0].pop, beam[nb - 1].pop);
    }
    printf("none W=%d up to H=%d\n", W, MAXH);
    return 0;
}
