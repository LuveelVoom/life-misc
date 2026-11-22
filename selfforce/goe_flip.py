"""Garden of Eden from a forcing agar: flip one output cell whose 3x3 stays forced.

If every predecessor of P on D minus {d} equals P on the 3x3 around d, then every predecessor
gives cell d the value f(P around d) = P_d (P is a still life). So P with d flipped has no
predecessor. Then: crop to the smallest window around d that is still an orphan, and greedily
unspecify cells while it stays an orphan.
usage: goe_flip.py RULE AGAR_ROWS W [--out FILE]"""
import argparse, random
import sfcore as S
from goe_sat import has_pred
from pysat.solvers import Solver

ap = argparse.ArgumentParser()
ap.add_argument("rule"); ap.add_argument("agar"); ap.add_argument("W", type=int)
ap.add_argument("--out"); ap.add_argument("--seed", type=int, default=0)
x = ap.parse_args()
rows = x.agar.strip("/").split("/"); a, b = len(rows), len(rows[0])
W = x.W
rc = S.RuleCNF(x.rule)
P = {(r, c): int(rows[r % a][c % b] == "o") for r in range(W) for c in range(W)}


def nbhd_forced(P, D, d):
    """All predecessors of P|D agree with P on the 3x3 around d?"""
    Dp = {(r + dr, c + dc) for (r, c) in D for dr, dc in S.OFFS}
    vid = {q: i + 1 for i, q in enumerate(sorted(Dp))}
    s = Solver(name="cadical153")
    for (r, c) in D:
        for cl in rc.clauses(lambda rr, cc: vid[(rr, cc)], r, c, P[(r, c)]):
            s.add_clause(cl)
    s.add_clause([(-vid[q] if P[q] else vid[q]) for q in
                  [(d[0] + dr, d[1] + dc) for dr, dc in S.OFFS]])
    res = s.solve(); s.delete()
    return not res


mid = W // 2
order = sorted(P, key=lambda q: max(abs(q[0] - mid), abs(q[1] - mid)))
d = None
for q in order[: 4 * a * b]:
    if nbhd_forced(P, set(P) - {q}, q):
        d = q; break
if d is None:
    raise SystemExit("no cell with a forced neighbourhood near the centre")
print("flip cell", d, "agar value", P[d], flush=True)
Q = dict(P); Q[d] ^= 1
assert has_pred(rc, Q) is None, "flipped window still has a predecessor?!"
print(f"{W}x{W} window with one flipped cell: orphan", flush=True)

# crop: smallest square-ish window around d that is still an orphan (shrink each side greedily)
r0, c0, r1, c1 = 0, 0, W, W
changed = True
while changed:
    changed = False
    for side in range(4):
        nr = [r0 + (side == 0), c0 + (side == 1), r1 - (side == 2), c1 - (side == 3)]
        if not (nr[0] <= d[0] < nr[2] and nr[1] <= d[1] < nr[3]): continue
        cells = [(r, c) for r in range(nr[0], nr[2]) for c in range(nr[1], nr[3])]
        if has_pred(rc, {q: Q[q] for q in cells}) is None:
            r0, c0, r1, c1 = nr; changed = True
box = [(r, c) for r in range(r0, r1) for c in range(c0, c1)]
print(f"cropped to {r1 - r0}x{c1 - c0}", flush=True)

# greedy unspecify (outermost first, random ties)
random.seed(x.seed)
spec = set(box)
for q in sorted(box, key=lambda q: (-max(abs(q[0] - d[0]), abs(q[1] - d[1])), random.random())):
    if q == d: continue
    trial = spec - {q}
    if has_pred(rc, {z: Q[z] for z in trial}) is None:
        spec = trial
print(f"orphan patch: {len(spec)} specified cells")
pic = "\n".join("".join(("o" if Q[(r, c)] else ".") if (r, c) in spec else " " for c in range(c0, c1))
                for r in range(r0, r1))
full = "\n".join("".join("o" if Q[(r, c)] else "." for c in range(c0, c1)) for r in range(r0, r1))
print("patch (blank = unspecified):\n" + pic)
print("full cropped window:\n" + full)
if x.out:
    with open(x.out, "w") as f:
        f.write(f"# goe_flip {x.rule} {x.agar} {W}; flip {d}; crop rows {r0}-{r1-1} cols {c0}-{c1-1}\n")
        f.write("# patch (blank = unspecified)\n" + pic + "\n# full window\n" + full + "\n")
