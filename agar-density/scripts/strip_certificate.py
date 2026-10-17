#!/usr/bin/env python3
"""
Exact proof (or counterexample) for ALL oscillating agars of period p that repeat every W columns.

Idea
  Such an agar is an infinite vertical sequence of rows; each row is W cells x p phases (W*p bits, wrapping
  horizontally).  Life's rule constrains only consecutive triples of rows.  Build the graph whose nodes are pairs of
  consecutive rows and whose edges are valid triples, with edge weight  2*(live cell-phases of the middle row) - W*p.
  An agar = a cycle in this graph; its density exceeds 1/2 exactly when the cycle's total weight is positive.

  * If there is NO positive cycle, Bellman-Ford converges to an integer potential g with
        g(a,b) >= w(a,b,c) + g(b,c)      for every valid triple (a,b,c).
    Summing around any cycle gives total weight <= 0, i.e. density <= 1/2.  The script re-checks this inequality on
    every edge in exact integer arithmetic and saves g as the certificate.
  * If there IS a positive cycle, the script extracts it and prints the agar (rows for every phase): a counterexample.

Covers every vertical period at once (any number of rows), all periods dividing p, and all widths dividing W.

Requirements: pip install numpy
Usage
  python strip_certificate.py --width 6 --period 2 --jobs 16          # the family-1 band's natural class
  python strip_certificate.py --width 4 --period 2                     # quick test (seconds)
  python strip_certificate.py --width 4 --period 2 --rule B3/S012345678   # sanity check: must find a counterexample
Memory/time (rough): W*p = 12 bits -> 4096 row states, 16.7M row pairs; expect several GB and 10-60 min with 16 jobs.
"""
import argparse, time, os
import numpy as np
from multiprocessing import Pool

T0 = time.time()
def log(*a): print(f"[{time.time()-T0:8.0f}s]", *a, flush=True)

def parse_rule(s):
    b, sv = s.upper().split('/')
    return {int(ch) for ch in b[1:]}, {int(ch) for ch in sv[1:]}

def setup(W, p, rule):
    global G_W, G_P, G_BORN, G_SURV, G_BITS, G_PROF
    G_W, G_P = W, p; G_BORN, G_SURV = rule
    S = 1 << (W * p)
    x = np.arange(S, dtype=np.int64)
    G_BITS = np.stack([np.stack([(x >> (t*W + j)) & 1 for j in range(W)], 1) for t in range(p)], 1).astype(np.int8)  # S x p x W
    G_PROF = (np.roll(G_BITS, 1, 2) + G_BITS + np.roll(G_BITS, -1, 2)).astype(np.int8)                               # 3-cell sums

def valid_pairs_for(b):
    """All (a, c) such that every cell of middle row b obeys the rule in every phase, given rows a above, c below."""
    W, p = G_W, G_P; bb = G_BITS[b]
    S = G_BITS.shape[0]; ok = np.ones((S, S), bool)
    for t in range(p):
        for j in range(W):
            cur = int(bb[t, j]); nxt = int(bb[(t+1) % p, j])
            inner = int(bb[t, (j-1) % W]) + int(bb[t, (j+1) % W])
            T = np.zeros((4, 4), bool)
            for sa in range(4):
                for sc in range(4):
                    n = inner + sa + sc
                    T[sa, sc] = (1 if ((n in G_SURV) if cur else (n in G_BORN)) else 0) == nxt
            if T.all(): continue
            ok &= T[G_PROF[:, t, j][:, None], G_PROF[:, t, j][None, :]]
            if not ok.any(): break
    a, c = np.nonzero(ok)
    return b, a.astype(np.int32), c.astype(np.int32)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=6); ap.add_argument("--period", type=int, default=2)
    ap.add_argument("--rule", default="B3/S23"); ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--check-every", type=int, default=25)
    a = ap.parse_args(); W, p = a.width, a.period; rule = parse_rule(a.rule)
    S = 1 << (W * p); V = S * S
    log(f"width {W}, period {p}, rule {a.rule}: {S} row states, {V} row pairs")
    setup(W, p, rule)
    live = G_BITS.reshape(S, -1).sum(1).astype(np.int64)
    src, dst, wt = [], [], []
    with Pool(a.jobs, initializer=setup, initargs=(W, p, rule)) as pool:
        for k, (b, aa, cc) in enumerate(pool.imap_unordered(valid_pairs_for, range(S), chunksize=8)):
            if len(aa):                                              # node ids fit in 32 bits for W*p <= 15
                idt = np.int32 if V < 2**31 else np.int64
                src.append((aa.astype(np.int64) * S + b).astype(idt)); dst.append((np.int64(b) * S + cc.astype(np.int64)).astype(idt))
                wt.append(np.full(len(aa), 2*live[b] - W*p, np.int16))
            if k % 256 == 0: log(f"   enumerated middle rows {k}/{S}, edges so far {sum(len(x) for x in src)}")
    src = np.concatenate(src); dst = np.concatenate(dst); wt = np.concatenate(wt).astype(np.int64)
    order = np.argsort(src, kind='stable'); src, dst, wt = src[order], dst[order], wt[order]
    us, starts = np.unique(src, return_index=True)
    log(f"graph: {len(wt)} valid triples (edges); max edge weight {wt.max()}")
    # Bellman-Ford for the longest path with "stop anywhere" (g >= 0); converges iff no positive cycle.
    g = np.zeros(V, np.int64)
    for it in range(1, 10**9):
        vals = wt + g[dst]                                     # int64 arithmetic
        best = np.maximum.reduceat(vals, starts)
        new = g.copy(); new[us] = np.maximum(g[us], best)
        if np.array_equal(new, g):
            ok = bool(np.all(g[src] >= wt + g[dst]))           # exact integer re-check of every edge
            np.save(f"strip_certificate_W{W}_p{p}.npy", g)
            log(f"CONVERGED after {it} rounds. Certificate check on all {len(wt)} edges: {ok}")
            if ok: log(f"PROVEN: every period-{p} agar repeating every {W} columns has density <= 1/2")
            return
        g = new
        if it % a.check_every == 0:
            # follow each node's best successor; a cycle with positive weight is a counterexample
            seg_argmax = np.zeros(len(us), np.int64)
            ordv = np.lexsort((vals, np.repeat(np.arange(len(us)), np.diff(np.append(starts, len(vals))))))
            last = np.append(starts[1:], len(vals)) - 1
            seg_argmax = ordv[last]
            succ = np.full(V, -1, np.int64); wsucc = np.zeros(V, np.int64)
            succ[us] = dst[seg_argmax]; wsucc[us] = wt[seg_argmax]
            u = int(us[np.argmax(g[us])])
            seen = {}
            path = []
            while u != -1 and u not in seen:
                seen[u] = len(path); path.append(u); u = int(succ[u])
            if u != -1:
                cyc = path[seen[u]:]; tot = int(sum(wsucc[x] for x in cyc))
                if tot > 0:
                    rows = [x % S for x in cyc]
                    log(f"COUNTEREXAMPLE: positive cycle of {len(rows)} rows, weight {tot} "
                        f"(density {(sum(live[r] for r in rows))}/{len(rows)*W*p})")
                    for t in range(p):
                        print(f"phase {t}:"); [print("".join("o" if G_BITS[r, t, j] else "." for j in range(W))) for r in rows]
                    return
            log(f"   round {it}: still improving (max potential {int(g.max())})")

if __name__ == "__main__":
    main()
