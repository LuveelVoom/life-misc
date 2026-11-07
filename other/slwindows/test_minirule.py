#!/usr/bin/env python3
"""minirule against casim (an independent Golly-compatible simulator). Optional: skipped
when casim isn't found (set CASIM=/path/to/casim-folder). Checks rule tables and canonical
names on random INT rules written in several spellings, one-generation steps on random
soups, periods, and RLE reading."""
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import minirule as M  # noqa: E402

CASIM = os.environ.get("CASIM", os.path.join(HERE, "..", "..", "..", "casim"))
sys.path.insert(0, CASIM)
try:
    import casim as C
except ImportError:
    print("casim not found (set CASIM); skipped")
    sys.exit(0)

rng = random.Random(1)


def random_rule():
    parts = []
    for c in "BS":
        s = c
        for k in range(9):
            if k in (0, 8) or not M.VALID.get(k):
                if rng.random() < 0.3 and not (c == "B" and k == 0):
                    s += str(k)
                continue
            r = rng.random()
            if r < 0.3:
                continue
            if r < 0.45:
                s += str(k)
                continue
            letters = [L for L in M.VALID[k] if rng.random() < 0.5]
            if not letters:
                continue
            s += str(k) + ("-" if rng.random() < 0.4 else "") + "".join(letters)
        parts.append(s)
    return "/".join(parts)


n = 0
for _ in range(3000):
    r = random_rule()
    for spelled in (r, r.lower(), r.replace("/", "")):
        a = M.parse(spelled)
        b = [int(v) for v in C.Rule(spelled).table]
        assert a == b, (spelled, "table")
        assert M.name(a) == C.Rule(spelled).name, (spelled, M.name(a), C.Rule(spelled).name)
        n += 1
print(f"rule tables and names: {n} spellings ok")

for t in range(300):
    r = random_rule() if t % 3 else "B3/S23"
    table = M.parse(r)
    if M.b0(table):
        continue
    cells = {(x, y) for x in range(12) for y in range(12) if rng.random() < 0.4}
    p = C.Pattern.from_coords(sorted(cells), r)
    p.step(1)
    assert M.step(cells, table) == {(int(x), int(y)) for x, y in p.coords()}, r
print("steps: ok")

for rle, rule in [("bo$2bo$3o!", "B3/S23"), ("2o$2o!", "B3/S23"), ("3o!", "B3/S23"),
                  ("bo$obo$bo!", "B3/S23-a5"), ("3o$o$bo!", "B36/S23")]:
    cells, _ = M.parse_rle(f"x = 3, y = 3, rule = {rule}\n{rle}")
    q = C.parse_rle(f"x = 3, y = 3, rule = {rule}\n{rle}")
    assert sorted(cells) == sorted((int(x), int(y)) for x, y in q.coords()), rle
    assert M.period(cells, M.parse(rule), 400) == q.period(400), (rle, rule)
print("rle and periods: ok")
