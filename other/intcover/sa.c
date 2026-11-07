// Annealer for intcover: W x H box, dead outside; cost = population + lambda * missing transitions.
// usage: ./sa W H seconds seed   (prints a line for each new best complete pattern, RLE-less: rows of .o)
// Heuristic: finds upper bounds only. Index: NW 256 ... SE 1, x right, y down.
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define MAXS 64
static int cls[512], cnt[102], W, H, bgc;
static uint8_t g[MAXS + 4][MAXS + 4];                      // cell (x,y) stored at [y+2][x+2]
static uint64_t st;
static inline uint64_t rnd(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static inline double urnd(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static void orbits(void) {
    memset(cls, -1, sizeof cls);
    int n = 0;
    for (int i = 0; i < 512; i++) {
        if (cls[i] >= 0) continue;
        int a[3][3], b[3][3];
        for (int k = 0; k < 9; k++) a[k / 3][k % 3] = (i >> (8 - k)) & 1;
        for (int t = 0; t < 8; t++) {
            int m = 0;
            for (int k = 0; k < 9; k++) m = 2 * m + a[k / 3][k % 3];
            cls[m] = n;
            if (t == 3) for (int r = 0; r < 3; r++) { int x = a[r][0]; a[r][0] = a[r][2]; a[r][2] = x; }
            else {
                memcpy(b, a, sizeof a);
                for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) a[r][c] = b[2 - c][r];
            }
        }
        n++;
    }
    if (n != 102) { fprintf(stderr, "orbits %d\n", n); exit(1); }
}
static inline int idx(int x, int y) {
    int i = 0;
    for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) i = 2 * i + g[y + dy + 2][x + dx + 2];
    return i;
}
static int missing(void) { int m = 0; for (int k = 0; k < 102; k++) m += cnt[k] == 0; return m; }
static void around(int x, int y, int d) {
    for (int dy = -1; dy <= 1; dy++) for (int dx = -1; dx <= 1; dx++) cnt[cls[idx(x + dx, y + dy)]] += d;
}
int main(int ac, char **av) {
    if (ac < 5) { fprintf(stderr, "usage: sa W H seconds seed\n"); return 2; }
    W = atoi(av[1]); H = atoi(av[2]); double secs = atof(av[3]); st = 88172645463325252ULL ^ (strtoull(av[4], 0, 10) * 0x9E3779B97F4A7C15ULL);
    if (W > MAXS || H > MAXS) return 2;
    orbits(); bgc = cls[0];
    int best = 1 << 30, bestm = 1000; const double lam = 3.0; time_t t0 = time(0);
    while (difftime(time(0), t0) < secs) {
        memset(g, 0, sizeof g); memset(cnt, 0, sizeof cnt);
        int pop = 0;
        for (int y = 0; y < H; y++) for (int x = 0; x < W; x++) if (rnd() & 1) { g[y + 2][x + 2] = 1; pop++; }
        for (int y = -1; y <= H; y++) for (int x = -1; x <= W; x++) cnt[cls[idx(x, y)]]++;
        cnt[bgc]++;
        int miss = missing();
        const long N = 40000000; const double T0 = 1.2, T1 = 0.08;
        for (long i = 0; i < N; i++) {
            double T = T0 * pow(T1 / T0, (double)i / N);
            int x = rnd() % W, y = rnd() % H, d = g[y + 2][x + 2] ? -1 : 1;
            around(x, y, -1); g[y + 2][x + 2] ^= 1; around(x, y, 1);
            int m2 = missing();
            double de = d + lam * (m2 - miss);
            if (de <= 0 || urnd() < exp(-de / T)) {
                pop += d; miss = m2;
                if (miss == 0 && pop < best) {
                    best = pop; printf("pop %d box %dx%d\n", pop, W, H);
                    for (int r = 0; r < H; r++) { for (int c = 0; c < W; c++) putchar(g[r + 2][c + 2] ? 'o' : '.'); putchar('\n'); }
                    fflush(stdout);
                }
                if (miss < bestm) bestm = miss;
            } else { around(x, y, -1); g[y + 2][x + 2] ^= 1; around(x, y, 1); }
        }
    }
    fprintf(stderr, "%dx%d seed %s: best pop %d, fewest missing seen %d\n", W, H, av[4], best == 1 << 30 ? -1 : best, bestm);
}
