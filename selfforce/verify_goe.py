"""Independent orphan check: kissat on 'some y has f(y) = P on the specified cells', one clause per
wrong 9-bit neighbourhood (no cube cover), table from casim.
Input: text picture, 'o' live, '.' dead, space = unspecified; '#' lines ignored. Section selectable
by --section N (0-based picture block separated by '#' lines). --unflip r,c flips one cell back."""
import sys, subprocess, tempfile, os, argparse
sys.path.insert(0, os.environ.get("CASIM", os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../casim")))
import casim as C
ap = argparse.ArgumentParser(); ap.add_argument("rule"); ap.add_argument("file")
ap.add_argument("--section", type=int, default=0); ap.add_argument("--unflip")
x = ap.parse_args()
secs, cur = [], None
for l in open(x.file):
    if l.startswith("#"):
        if cur: secs.append(cur)
        cur = []
    else:
        cur = (cur or []) + [l.rstrip("\n")]
if cur: secs.append(cur)
lines = secs[x.section]
P = {(r, c): int(ch == "o") for r, l in enumerate(lines) for c, ch in enumerate(l) if ch in "o."}
if x.unflip:
    r, c = map(int, x.unflip.split(",")); P[(r, c)] ^= 1
T = C.Rule(x.rule).table
offs = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]
vid = {}
v = lambda q: vid.setdefault(q, len(vid) + 1)
cls = [[(-v((r + dr, c + dc)) if (i >> (8 - k)) & 1 else v((r + dr, c + dc))) for k, (dr, dc) in enumerate(offs)]
       for (r, c), val in P.items() for i in range(512) if T[i] != val]
with tempfile.NamedTemporaryFile("w", suffix=".cnf", delete=False) as f:
    f.write(f"p cnf {len(vid)} {len(cls)}\n" + "".join(" ".join(map(str, cl)) + " 0\n" for cl in cls))
out = subprocess.run(["kissat", "-q", f.name], capture_output=True, text=True).stdout; os.unlink(f.name)
s = [l for l in out.splitlines() if l.startswith("s ")][0]
print(f"{len(P)} specified cells, population {sum(P.values())}: {s} ->", "ORPHAN" if "UNSAT" in s else "has predecessor")
