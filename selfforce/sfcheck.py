"""Maximal self-forcing subset of a still life (or any pattern) plus a dead margin.
usage: sfcheck.py file.rle [--rule R] [--margin M]"""
import argparse, time
import numpy as np
import sfcore as S

def grid_patch(g, margin):
    h, w = g.shape
    P = {}
    for r in range(-margin, h + margin):
        for c in range(-margin, w + margin):
            P[(r, c)] = int(g[r, c]) if 0 <= r < h and 0 <= c < w else 0
    return P

def show(P, D):
    rs = [r for r, _ in D]; cs = [c for _, c in D]
    out = []
    for r in range(min(rs), max(rs) + 1):
        out.append("".join(("o" if P[(r, c)] else ".") if (r, c) in D else " "
                           for c in range(min(cs), max(cs) + 1)))
    return "\n".join(out)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("rle"); ap.add_argument("--rule"); ap.add_argument("--margin", type=int, default=2)
    a = ap.parse_args()
    p = S.C.parse_rle(open(a.rle).read(), a.rule)
    rule = a.rule or p.rule.name
    g = np.asarray(p.cells)
    q = S.C.parse_rle(open(a.rle).read(), a.rule); q.step(1)
    print("rule", rule, "still life:", np.array_equal(np.asarray(q.cells), g) and q.bbox == p.bbox)
    P = grid_patch(g, a.margin)
    t = time.time()
    D = S.max_self_forcing(rule, P, set(P))
    print(f"max self-forcing subset: {len(D)} cells ({time.time()-t:.1f}s)")
    if D: print(show(P, D))
