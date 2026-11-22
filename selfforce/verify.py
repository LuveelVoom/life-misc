"""Independent check of a self-forcing still life.
Input: picture file from sfcegar --out ('O'/'-' = live/dead in D, 'o'/'.' = live/dead outside D;
lines starting with '#' ignored). Checks:
  1. casim: the whole pattern is a still life.
  2. kissat: no y on D+ with f(y)|D = P|D and y != P somewhere on D. The encoding forbids every
     wrong 9-bit neighbourhood one clause each (no cube minimisation), table from casim.
usage: verify.py RULE FILE"""
import sys, subprocess, tempfile, os
import numpy as np
sys.path.insert(0, os.environ.get("CASIM", os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../casim")))
import casim as C

rule, fn = sys.argv[1], sys.argv[2]
lines = [l.rstrip("\n") for l in open(fn) if l.strip() and not l.startswith("#")]
H, W = len(lines), max(map(len, lines))
P, D = {}, set()
for r, l in enumerate(lines):
    for c, ch in enumerate(l.ljust(W, ".")):
        P[(r, c)] = int(ch in "Oo")
        if ch in "O-": D.add((r, c))
g = np.array([[P[(r, c)] for c in range(W)] for r in range(H)], dtype=np.uint8)
p = C.Pattern.from_coords([(c, r) for (r, c), v in P.items() if v], rule)
before = p.to_rle(); p.step(1)
print("still life (casim):", p.to_rle() == before, " population", int(g.sum()), " |D| =", len(D))
table = C.Rule(rule).table
offs = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
vid = {}
def v(q):
    if q not in vid: vid[q] = len(vid) + 1
    return vid[q]
cls = []
for (r, c) in D:
    for i in range(512):
        if table[i] != P[(r, c)]:
            cls.append([(-v((r + dr, c + dc)) if (i >> (8 - k)) & 1 else v((r + dr, c + dc)))
                        for k, (dr, dc) in enumerate(offs)])
cls.append([(-v(q) if P[q] else v(q)) for q in sorted(D)])
with tempfile.NamedTemporaryFile("w", suffix=".cnf", delete=False) as f:
    f.write(f"p cnf {len(vid)} {len(cls)}\n")
    for cl in cls: f.write(" ".join(map(str, cl)) + " 0\n")
    name = f.name
out = subprocess.run(["kissat", "-q", name], capture_output=True, text=True).stdout
os.unlink(name)
res = [l for l in out.splitlines() if l.startswith("s ")]
print("kissat on 'predecessor differing on D':", res[0] if res else out[:200])
print("SELF-FORCING" if res and "UNSATISFIABLE" in res[0] else "NOT VERIFIED")
