"""Check rule clauses and forced_cells against brute-force enumeration on tiny domains."""
import random
import numpy as np
import sfcore as S


def brute_forced(table, P, D):
    Dp = sorted({(r + dr, c + dc) for (r, c) in D for dr, dc in S.OFFS})
    n = len(Dp)
    idx = {x: i for i, x in enumerate(Dp)}
    ys = ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1).astype(np.int64)
    ok = np.ones(1 << n, bool)
    T = np.array(table)
    for (r, c) in D:
        k = np.zeros(1 << n, np.int64)
        for b, (dr, dc) in enumerate(S.OFFS):
            k |= ys[:, idx[(r + dr, c + dc)]] << (8 - b)
        ok &= T[k] == P[(r, c)]
    good = ys[ok]
    return {x for x in D if (good[:, idx[x]] == P[x]).all()}


random.seed(1)
rnd = [0] + [random.randint(0, 1) for _ in range(511)]
rnd2 = [0] + [int(random.random() < 0.15) for _ in range(511)]
rules = ["B3/S23", "B3/S23-a5", "B2-a3/S12-a", S.C.map_string(rnd), S.C.map_string(rnd2),
         "B/S012345678", "B1/S1", "B3/S012345678"]
shapes = [[(0, 0), (0, 1), (1, 0), (1, 1)], [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)],
          [(0, 0), (1, 1), (0, 2)], [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]]
for rule in rules:
    rc = S.RuleCNF(rule)
    for i in range(512):
        for v in (0, 1):
            viol = any(all(((i >> (8 - k)) & 1) == b for k, b in cube.items()) for cube in rc.cubes[v])
            assert viol == (rc.table[i] != v), (rule, i, v)
    bad = 0
    nforced = 0
    for t in range(60):
        D = random.choice(shapes)
        if t % 2:
            P = {x: random.randint(0, 1) for x in D}
        else:  # an image, so it has predecessors
            y = {(r, c): random.randint(0, 1) for r in range(-1, 3) for c in range(-1, 4)}
            P = {(r, c): rc.table[sum(y[(r + dr, c + dc)] << (8 - b) for b, (dr, dc) in enumerate(S.OFFS))]
                 for (r, c) in D}
        a = S.forced_cells(rc, P, set(D))
        b = brute_forced(rc.table, P, D)
        nforced += len(b)
        if a != b:
            bad += 1
            print("MISMATCH", rule, P, a, b)
    print(rule[:20], "ok" if not bad else f"{bad} bad", f"(forced cells seen: {nforced})")
