"""Enumerate still-life agars on an a x b torus (rows x cols, up to translation), tile each into a
W x W window, and report forced cells and the maximal self-forcing subset of the window.
Prints each agar that beats the best forced count so far; stops at the first self-forcing subset.
usage: agars.py RULE a b [--W 24] [--max N] [--live 35 --dead 45]"""
import argparse, time, sys
import sfcore as S
from pysat.solvers import Solver

def canon_torus(g, a, b):
    best = None
    for dr in range(a):
        for dc in range(b):
            t = tuple(g[(r + dr) % a][(c + dc) % b] for r in range(a) for c in range(b))
            if best is None or t < best: best = t
    return best

def count_cubes(live_counts, dead_counts):
    """Forbidden-cube cover for 'live centre has count in live_counts, dead centre in dead_counts'."""
    ok = [0] * 512
    for i in range(512):
        n = bin(i & ~16).count("1")
        allowed = live_counts if i & 16 else dead_counts
        ok[i] = int(allowed is None or n in allowed)
    return S._cubes(ok, 1)

def enum_agars(rc, a, b, maxn, live_counts=None, dead_counts=None):
    v = lambda r, c: (r % a) * b + (c % b) + 1
    s = Solver(name="cadical153")
    if live_counts is not None or dead_counts is not None:
        for cube in count_cubes(live_counts, dead_counts):
            for r in range(a):
                for c in range(b):
                    s.add_clause([(-v(r + dr, c + dc) if bit else v(r + dr, c + dc))
                                  for k, bit in cube.items() for dr, dc in [S.OFFS[k]]])
    for r in range(a):
        for c in range(b):
            for val in (0, 1):
                # still life: if cell==val then next==val
                for cl in rc.clauses(v, r, c, val):
                    s.add_clause(cl + ([-v(r, c)] if val else [v(r, c)]))
    seen = set()
    n = 0
    while s.solve():
        m = s.get_model()
        g = [[int(m[v(r, c) - 1] > 0) for c in range(b)] for r in range(a)]
        # block all translates
        for dr in range(a):
            for dc in range(b):
                s.add_clause([(-v(r, c) if g[(r + dr) % a][(c + dc) % b] else v(r, c))
                              for r in range(a) for c in range(b)])
        k = canon_torus(g, a, b)
        if k in seen: continue
        seen.add(k)
        yield g
        n += 1
        if n >= maxn: break

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("rule"); ap.add_argument("a", type=int); ap.add_argument("b", type=int)
    ap.add_argument("--W", type=int, default=24); ap.add_argument("--max", type=int, default=100000)
    ap.add_argument("--live", help="allowed live-cell neighbour counts, e.g. 35 (default any)")
    ap.add_argument("--dead", help="allowed dead-cell neighbour counts, e.g. 45 (default any)")
    x = ap.parse_args()
    rc = S.RuleCNF(x.rule)
    t0 = time.time(); best = (-1, None); cnt = 0
    lc = None if x.live is None else {int(ch) for ch in x.live}
    dc = None if x.dead is None else {int(ch) for ch in x.dead}
    for g in enum_agars(rc, x.a, x.b, x.max, lc, dc):
        cnt += 1
        if sum(map(sum, g)) == 0: continue
        P = {(r, c): g[r % x.a][c % x.b] for r in range(x.W) for c in range(x.W)}
        F = S.forced_cells(rc, P, set(P))
        if len(F) > best[0]:
            best = (len(F), g)
            print(f"  agar #{cnt}: forced {len(F)}/{len(P)}  {''.join(''.join('o' if q else '.' for q in row)+'/' for row in g)}", flush=True)
        if F:
            D = S.max_self_forcing(rc, P, F, log=None)
            if D:
                print(f"  SELF-FORCING subset {len(D)} cells in agar #{cnt}", flush=True)
                import sfcheck; print(sfcheck.show(P, D))
                break
    print(f"{x.rule} {x.a}x{x.b} live={x.live} dead={x.dead}: {cnt} agars, best forced {best[0]}, {time.time()-t0:.1f}s")
