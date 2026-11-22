"""CEGAR search for a still life S containing a self-forcing patch D.

Box of H x W cells (rows x cols) with vacuum outside: one ring of fixed-dead cells that must stay
dead. Optional agar fixed on a core rectangle. D = box cells at Chebyshev distance >= inset from
the box border (inset 0 = the whole box).

Generator SAT: S is a still life (box + ring), agar core fixed, plus learned clauses.
Checker SAT: y on D+ with f(y)|D = S|D and y != S somewhere on D. UNSAT => D is self-forcing
and S is written to --out. Otherwise the checker looks for the y with the fewest differences
(at most K differing cells, K = 1, 2, 4, ..., up to --kmax). Its differences are split into
clusters (cells within distance 2 joined); the smallest cluster touching D is a counterexample
on its own and stays one for every S' that agrees with S on the cluster dilated by 2. That local
configuration is blocked at every translate, and under the 8 symmetries for isotropic rules,
where each involved cell has the same status (outside box / box minus D / D).

--rigid 1 forbids predecessors that differ from S in one cell, --rigid 2 also in two cells within
distance 2, encoded up front (the loop would find and block exactly these).
--live/--dead restrict neighbour counts of live/dead box cells away from the box edge.
--minpop excludes near-empty still lifes. --wrap: columns wrap with period W (an infinite band).
Generator UNSAT: no still life in this box has D self-forcing.

usage: sfcegar.py RULE H W [--inset k] [--rigid 2] [--minpop N] [--agar ROWS --core r0 c0 r1 c1]
                  [--live 35 --dead 45] [--wrap] [--out FILE]
"""
import argparse, time, random, sys
import numpy as np
import sfcore as S
from pysat.solvers import Solver
from agars import count_cubes
from pysat.card import ITotalizer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rule"); ap.add_argument("H", type=int); ap.add_argument("W", type=int)
    ap.add_argument("--agar"); ap.add_argument("--core", type=int, nargs=4)
    ap.add_argument("--agar-off", type=int, nargs=2, default=[0, 0], help="agar phase (row, col)")
    ap.add_argument("--inset", type=int, default=0)
    ap.add_argument("--wrap", action="store_true")
    ap.add_argument("--live"); ap.add_argument("--dead")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--maxiter", type=int, default=10**9)
    ap.add_argument("--minpop", type=int, default=1)
    ap.add_argument("--out")
    ap.add_argument("--rigid", type=int, default=0, help="encode 1-flip (1) or 1+2-flip (2) rigidity up front")
    ap.add_argument("--nosym", action="store_true", help="block only the counterexample's own position")
    ap.add_argument("--kmax", type=int, default=64, help="largest difference bound tried before taking any counterexample")
    x = ap.parse_args()
    H, W = x.H, x.W
    rc = S.RuleCNF(x.rule)
    random.seed(x.seed)
    wrap = x.wrap
    nm = (lambda q: (q[0], q[1] % W)) if wrap else (lambda q: q)
    box = [(r, c) for r in range(H) for c in range(W)]
    inbox = lambda q: 0 <= q[0] < H and (wrap or 0 <= q[1] < W)

    # generator variables: box cells; outside = constant 0
    gv = {q: i + 1 for i, q in enumerate(box)}
    TRUE = len(gv) + 1
    def g(r, c):
        q = nm((r, c))
        return gv[q] if inbox(q) else -TRUE
    gen = Solver(name="cadical153")
    gen.add_clause([TRUE])
    ring = set()
    for (r, c) in box:
        for dr, dc in S.OFFS:
            q = nm((r + dr, c + dc))
            if not inbox(q): ring.add(q)
    for q in box + sorted(ring):
        r, c = q
        for val in (0, 1):
            if q in ring and val == 1: continue
            for cl in rc.clauses(g, r, c, val):
                cl2 = []
                sat = False
                for l in cl:
                    if l == TRUE: sat = True; break
                    if l == -TRUE: continue
                    cl2.append(l)
                if sat: continue
                gen.add_clause(cl2 + ([-g(r, c)] if val else [g(r, c)]) if q not in ring else cl2)
    fixed = set()
    if x.agar:
        rows = x.agar.strip("/").split("/"); a, b = len(rows), len(rows[0])
        r0, c0, r1, c1 = x.core
        for r in range(r0, r1):
            for c in range(c0, c1):
                bit = rows[(r + x.agar_off[0]) % a][(c + x.agar_off[1]) % b] == 'o'
                gen.add_clause([gv[(r, c)] if bit else -gv[(r, c)]]); fixed.add((r, c))
    if x.live:
        cubes = count_cubes({int(ch) for ch in x.live}, None if x.dead is None else {int(ch) for ch in x.dead})
        for (r, c) in box:
            if (r, c) in fixed: continue
            if min(r, H - 1 - r, *( [] if wrap else [c, W - 1 - c])) < 1: continue
            for cube in cubes:
                cl = [g(r + S.OFFS[k][0], c + S.OFFS[k][1]) * (-1 if bit else 1) for k, bit in cube.items()]
                if any(l == TRUE for l in cl): continue
                gen.add_clause([l for l in cl if l != -TRUE])
    if x.minpop > 1:
        from pysat.card import CardEnc, EncType
        enc = CardEnc.atleast([gv[q] for q in box], bound=x.minpop, top_id=TRUE, encoding=EncType.seqcounter)
        for cl in enc.clauses: gen.add_clause(cl)
    dist = lambda q: min(q[0], H - 1 - q[0]) if wrap else min(q[0], H - 1 - q[0], q[1], W - 1 - q[1])
    D = [q for q in box if dist(q) >= x.inset]
    Dset = set(D)
    Dp = sorted({nm((r + dr, c + dc)) for (r, c) in D for dr, dc in S.OFFS})
    free = [q for q in box if q not in fixed]
    print(f"box {H}x{W} wrap={wrap} |D|={len(D)} free={len(free)} fixed={len(fixed)}", flush=True)

    # Up-front rigidity: no predecessor differs from S in exactly 1 (or 2 nearby) cells.
    # For flip set F (touching D), some d in D near F must change output when F is flipped:
    # aux a_{F,d} -> neighbourhood of d lies in Sens(mask) = {n : T(n ^ mask) != centre(n)}.
    # Sound: these are exactly the counterexamples with |diff| = |F|, which the loop would find
    # and block anyway (pairs farther apart than 2 act independently, so singles cover them).
    if x.rigid:
        nv = max(TRUE, gen.nof_vars() if hasattr(gen, "nof_vars") else TRUE)
        if x.minpop > 1: nv = max(nv, enc.nv)
        cube_cache = {}
        def insens_cubes(mask):
            if mask not in cube_cache:
                sens = [int(rc.table[n ^ mask] != ((n >> 4) & 1)) for n in range(512)]
                cube_cache[mask] = S._cubes(sens, 1)
            return cube_cache[mask]
        flipsets = [(z,) for z in D]
        if x.rigid >= 2:
            seenp = set()
            for z1 in D:
                for dr in range(-2, 3):
                    for dc in range(-2, 3):
                        z2 = nm((z1[0] + dr, z1[1] + dc))
                        if z2 == z1: continue
                        key = frozenset((z1, z2))
                        if key in seenp: continue
                        seenp.add(key); flipsets.append((z1, z2))
        nrc = 0; aux_cache = {}
        for F in flipsets:
            ds = {nm((f[0] + dr, f[1] + dc)) for f in F for dr, dc in S.OFFS}
            ors = []
            for d in ds:
                if d not in Dset: continue
                mask = 0
                for k, (dr, dc) in enumerate(S.OFFS):
                    if nm((d[0] + dr, d[1] + dc)) in F: mask |= 1 << (8 - k)
                if (d, mask) in aux_cache:
                    ors.append(aux_cache[(d, mask)]); continue
                nv += 1; a = nv; ors.append(a); aux_cache[(d, mask)] = a
                for cube in insens_cubes(mask):
                    cl = [-a]; sat = False
                    for k, bit in cube.items():
                        l = g(d[0] + S.OFFS[k][0], d[1] + S.OFFS[k][1])
                        l = -l if bit else l
                        if l == TRUE: sat = True; break
                        if l != -TRUE: cl.append(l)
                    if not sat: gen.add_clause(cl); nrc += 1
            gen.add_clause(ors); nrc += 1
        print(f"rigidity {x.rigid}: {len(flipsets)} flip sets, {nrc} clauses", flush=True)

    def status(q):
        q = nm(q)
        return 0 if not inbox(q) else (2 if q in Dset else 1)
    syms = [lambda d: (d[0], d[1]), lambda d: (-d[1], d[0]), lambda d: (-d[0], -d[1]), lambda d: (d[1], -d[0]),
            lambda d: (d[0], -d[1]), lambda d: (-d[0], d[1]), lambda d: (d[1], d[0]), lambda d: (-d[1], -d[0])]
    if not x.nosym and not S.C.is_isotropic(np.array(rc.table, dtype=np.uint8)):
        syms = syms[:1]
    if x.nosym: syms = []
    seen_cl = set()
    def add_block(nb, val):
        """nb: cells (normalised), val: S value per cell. Add clause and all its valid images."""
        nb = list(nb)
        base = [(q, status(q), val[q]) for q in nb]
        added = 0
        if any(q in fixed for q in nb) or not syms:
            cl = [(-gv[q] if val[q] else gv[q]) for q in nb if inbox(q) and q not in fixed]
            gen.add_clause(cl); return 1
        a0 = nb[0]
        rels = [((q[0] - a0[0], q[1] - a0[1]), st, v) for q, st, v in base]
        for sy in syms:
            img = [(sy(d), st, v) for d, st, v in rels]
            for tr in range(-3, H + 3):
                for tc in range(-3, W + 3):
                    ok = True; cl = []
                    for (dr, dc), st, v in img:
                        q = nm((tr + dr, tc + dc))
                        if status(q) != st or q in fixed: ok = False; break
                        if st: cl.append(-gv[q] if v else gv[q])
                    if ok:
                        key = frozenset(cl)
                        if key not in seen_cl:
                            seen_cl.add(key); gen.add_clause(cl); added += 1
        return added

    t0 = time.time(); it = 0; blk_sizes = []; nclauses = 0
    while it < x.maxiter:
        gen.set_phases([random.choice((-1, 1)) * gv[q] for q in free[:: max(1, len(free) // 64)]])
        if not gen.solve():
            print(f"GENERATOR UNSAT after {it} iterations: no such still life ({time.time()-t0:.0f}s)")
            return
        m = set(l for l in gen.get_model() if l > 0)
        Sv = lambda q: (1 if gv[nm(q)] in m else 0) if inbox(nm(q)) else 0
        # checker
        vid = {q: i + 1 for i, q in enumerate(Dp)}
        yv = lambda r, c: vid[nm((r, c))]
        ch = Solver(name="cadical153")
        for (r, c) in D:
            for cl in rc.clauses(yv, r, c, Sv((r, c))):
                ch.add_clause(cl)
        ch.add_clause([(-vid[q] if Sv(q) else vid[q]) for q in D])
        it += 1
        if not ch.solve():
            ch.delete()
            P = {q: Sv(q) for q in box}
            print(f"FOUND after {it} iterations ({time.time()-t0:.0f}s): pop {sum(P.values())}")
            pic = "\n".join("".join(("O" if P[(r, c)] else "-") if (r, c) in Dset else ("o" if P[(r, c)] else ".")
                                    for c in range(W)) for r in range(H))
            print(pic, flush=True)
            if x.out:
                with open(x.out, "w") as f:
                    f.write(f"# sfcegar {' '.join(sys.argv[1:])}\n# D: inset {x.inset}\n{pic}\n")
            return
        # smallest-difference counterexample: at most K differing cells, K = 1, 2, 4, ...
        ym = set(l for l in ch.get_model() if l > 0)
        dl = [(-vid[q] if Sv(q) else vid[q]) for q in Dp]
        nd = sum(1 for q in Dp if (vid[q] in ym) != bool(Sv(q)))
        tot = ITotalizer(lits=dl, ubound=min(nd, x.kmax), top_id=len(vid) + 1)
        for cl in tot.cnf.clauses: ch.add_clause(cl)
        K = 1
        while K < min(nd, x.kmax + 1):
            if ch.solve(assumptions=[-tot.rhs[K]]):
                ym = set(l for l in ch.get_model() if l > 0); break
            K *= 2
        ch.delete()
        diff = [q for q in Dp if (vid[q] in ym) != bool(Sv(q))]
        # cluster: cells within Chebyshev distance 2 joined
        dset = set(diff); comp = {}
        for q in diff:
            if q in comp: continue
            stack = [q]; comp[q] = q
            while stack:
                u = stack.pop()
                for dr in range(-2, 3):
                    for dc in range(-2, 3):
                        w = nm((u[0] + dr, u[1] + dc))
                        if w in dset and w not in comp:
                            comp[w] = q; stack.append(w)
        groups = {}
        for q in diff: groups.setdefault(comp[q], []).append(q)
        cands = [gr for gr in groups.values() if any(q in Dset for q in gr)]
        gr = min(cands, key=len)
        nb = {nm((r + dr, c + dc)) for (r, c) in gr for dr in range(-2, 3) for dc in range(-2, 3)}
        cl = [(-gv[q] if gv[q] in m else gv[q]) for q in nb if inbox(q) and q not in fixed]
        blk_sizes.append(len(cl))
        if not cl:
            print(f"local counterexample touches only fixed cells: the fixed part is not rigid ({len(gr)} diff cells)")
            return
        nclauses += add_block(nb, {q: Sv(q) for q in nb})
        if it % 200 == 0 or it < 5:
            bs = blk_sizes[-200:]
            print(f"it {it}: last diff cluster {len(gr)} cells, block {len(cl)} lits, "
                  f"avg block {sum(bs)/len(bs):.1f}, clauses {nclauses}, {time.time()-t0:.0f}s", flush=True)
    print("maxiter reached")


if __name__ == "__main__":
    main()
