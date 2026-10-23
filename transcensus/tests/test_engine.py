#!/usr/bin/env python3
"""Engine check: transcensus-int --evolve against casim on random soups.

Covers Life, totalistic, INT and non-isotropic MAP rules (including a pure shift rule that pins the
orientation: next state = NE neighbour, so patterns move one cell SW per generation). Soups are placed
at x offsets that put them on and across the 64-bit word boundaries of the census grid. Generation
counts stay below the distance to the grid edge, so casim's unbounded plane gives the same answer.
With --generated it also builds and checks the fixed-rule (make RULE=...) path for every rule.
"""
import os, random, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "..", "..", "casim"))
import numpy as np
import casim as C

rng = random.Random(11)
shift_ne = np.array([(x >> 6) & 1 for x in range(512)], np.uint8)     # MAP weight 64 = NE
rand_map = np.array([rng.random() < 0.35 for _ in range(512)], np.uint8); rand_map[0] = 0
RULES = ["B3/S23", "B36/S23", "B2-a3/S23", "B3-cnqy4z/S23-a4i", "B2e3ai4n/S12k3-q4ent5",
         "B35y6ac/S2-in3-e4cy", "B2ci3ai4c8/S02ae3eijkq4iz5ar6i7e", C.map_string(shift_ne), C.map_string(rand_map)]

def rle_from_cells(cells):
    h = max(y for x, y in cells) + 1
    rows = []
    for y in range(h):
        xs = sorted(x for x, yy in cells if yy == y)
        row, cx = "", 0
        for x in xs:
            if x > cx: row += f"{x - cx}b"
            row += "o"; cx = x + 1
        rows.append(row)
    return "x = 0, y = 0\n" + "$".join(rows) + "!\n"

def run(exe, rule, cells, gens):
    p = subprocess.run([exe, "-r", rule, "--evolve", str(gens)], input=rle_from_cells(cells),
                       capture_output=True, text=True, check=True)
    return sorted((int(a), int(b)) for a, b in (l.split() for l in p.stdout.split("\n") if l))

def check(exe, label):
    bad = cases = 0
    for rule in RULES:
        for _ in range(12):
            dx = rng.choice([0, 20, 26, 31, 32, 40, 58])
            size = rng.choice([4, 8, 12])
            cells = [(dx + x, y) for y in range(size) for x in range(size) if rng.random() < 0.45]
            cells.append((dx, 0))           # fix the RLE origin so placement is exact
            cells = sorted(set(cells))
            gens = rng.choice([1, 2, 7, 25])
            if dx + size + gens >= 95:
                gens = 95 - dx - size - 1
            got = run(exe, rule, cells, gens)
            ref = C.Pattern.from_coords(cells, rule, background="literal")
            ref.step(gens)
            want = sorted(ref.coords())
            cases += 1
            if got != want:
                bad += 1
                print(f"[{label}] MISMATCH rule={rule[:30]} dx={dx} size={size} gens={gens}: {len(got)} vs {len(want)} cells")
    print(f"[{label}] {cases} cases, {bad} mismatches")
    return bad

bad = check(os.path.join(ROOT, "transcensus-int"), "generic")
if "--generated" in sys.argv:
    with tempfile.TemporaryDirectory() as d:
        for rule in RULES:
            h = os.path.join(d, "rule.h")
            with open(h, "w") as f:
                subprocess.run([os.path.join(ROOT, "transcensus-int"), "-r", rule, "--emit-c"], stdout=f, check=True)
            exe = os.path.join(d, "tc-gen")
            subprocess.run(["clang", "-std=c99", "-O3", "-mcpu=native", "-Wno-unused-function", f'-DRULE_GENERATED="{h}"',
                            os.path.join(ROOT, "transcensus-int.c"), "-o", exe, "-lm"], check=True)
            saved = RULES[:]
            RULES[:] = [rule]
            bad += check(exe, "generated " + rule[:24])
            RULES[:] = saved
sys.exit(1 if bad else 0)
