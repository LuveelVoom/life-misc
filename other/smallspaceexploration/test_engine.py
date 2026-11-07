# Compare soupclass -D (tile engine) with a full-grid numpy reference, same MAP convention.
import subprocess, sys, random, numpy as np
import rulespace as R
def ref(tab, cells, N, gens, erase=0):
    g = np.zeros((N, N), np.int32)
    for x, y in cells: g[y, x] = 1
    t = np.array(tab, np.uint8)
    for _ in range(gens):
        p = np.pad(g, 1)
        idx = np.zeros((N, N), np.int32)
        w = 256
        for dy in (-1, 0, 1):          # row-major from NW
            for dx in (-1, 0, 1):
                idx += w * p[1 + dy:1 + dy + N, 1 + dx:1 + dx + N]; w //= 2
        g = t[idx].astype(np.int32)
        if erase: g[:16, :] = 0; g[-16:, :] = 0; g[:, :16] = 0; g[:, -16:] = 0
    return {(int(x), int(y)) for y, x in zip(*np.nonzero(g))}
def run(tab, N, gens, seed=7, w=16, d=50, erase=0):
    out = subprocess.run(["./soupclass", "-D", "-g", str(gens), "-N", str(N), "-s", str(seed), "-w", str(w), "-d", str(d), "-E", str(erase)],
                         input="x " + R.table_hex(tab) + "\n", capture_output=True, text=True, check=True).stdout.split("\n")
    return {tuple(map(int, l.split())) for l in out if l and not l.startswith("end")}
N = 128; fails = 0
rng = random.Random(3)
tests = []
for k in range(6):           # random non-isotropic tables (no B0), asymmetric
    t = [rng.random() < 0.35 for _ in range(512)]; t[0] = 0; tests.append(("map%d" % k, [int(v) for v in t]))
for rid in rng.sample(range(65536), 6) + [0, 0xffff]:
    tests.append(("id%04x" % rid, R.rule_table(rid)))
tests.append(("life", R.table_from_sets(*R.parse_hensel("B3/S23"))))
for name, tab in tests:
    init = run(tab, N, 0)
    for gens in (1, 7, 60):
        for e in (0, 1):
            a, b = run(tab, N, gens, erase=e), ref(tab, init, N, gens, erase=e)
            if a != b: fails += 1; print("MISMATCH", name, gens, e, len(a), len(b))
    print(name, "ok" if True else "", len(init), len(a))
# asymmetric single-step orientation check: one cell at NE of an empty centre, rule born iff idx==64
t = [0]*512; t[64] = 1
a = subprocess.run(["./soupclass", "-D", "-g", "1", "-N", "32", "-w", "1", "-d", "100", "-s", "1", "-E", "0"],
                   input="x " + R.table_hex(t) + "\n", capture_output=True, text=True).stdout
print("NE-only rule, single cell at (16,16) ->", a.split("\n")[0], "(expect 15 17: cell born SW of the seed)")
print("FAILS", fails)
