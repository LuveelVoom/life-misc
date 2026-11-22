"""Core: rule clauses and maximal self-forcing subset of a patch.

Definition (Torma-Salo): patch P on finite domain D is self-forcing if every y with
f(y)|D = P has y|D = P. f(y)|D depends only on y on D+ (D dilated by the Moore
neighbourhood), so this is a finite UNSAT check over y on D+.

Cells are (r, c) = (row, column); row grows down. MAP index: NW 256 ... SE 1, centre 16.
"""
import os, sys
sys.path.insert(0, os.environ.get("CASIM", os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../casim")))
import casim as C
from pysat.solvers import Solver

OFFS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]  # bit 8..0


def rule_table(rule):
    return list(C.Rule(rule).table)


def _cubes(table, v):
    """Cubes (dict pos->bit) covering exactly the 9-bit inputs whose output != v.
    Greedy expansion; each cube is a forbidden partial assignment."""
    bad = {i for i in range(512) if table[i] != v}
    left = set(bad)
    cubes = []
    while left:
        m = min(left)
        cube = {k: (m >> (8 - k)) & 1 for k in range(9)}
        for k in range(9):
            trial = {kk: b for kk, b in cube.items() if kk != k}
            free = [kk for kk in range(9) if kk not in trial]
            ok = True
            for s in range(1 << len(free)):
                idx = 0
                for kk, b in trial.items():
                    idx |= b << (8 - kk)
                for j, kk in enumerate(free):
                    idx |= ((s >> j) & 1) << (8 - kk)
                if idx not in bad:
                    ok = False
                    break
            if ok:
                cube = trial
        cubes.append(cube)
        free = [kk for kk in range(9) if kk not in cube]
        for s in range(1 << len(free)):
            idx = 0
            for kk, b in cube.items():
                idx |= b << (8 - kk)
            for j, kk in enumerate(free):
                idx |= ((s >> j) & 1) << (8 - kk)
            left.discard(idx)
    return cubes


class RuleCNF:
    def __init__(self, rule):
        self.rule = rule
        self.table = rule_table(rule)
        if self.table[0]:
            raise ValueError("B0 rules not supported")
        self.cubes = {0: _cubes(self.table, 0), 1: _cubes(self.table, 1)}

    def clauses(self, var, r, c, v):
        """Clauses forcing next-gen cell (r,c) == v; var(r,c) -> literal of y."""
        out = []
        for cube in self.cubes[v]:
            cl = []
            for k, b in cube.items():
                dr, dc = OFFS[k]
                x = var(r + dr, c + dc)
                cl.append(-x if b else x)
            out.append(cl)
        return out


def forced_cells(rc, P, D, solver="cadical153", log=None, norm=None):
    """Cells of D forced (every predecessor y of P|D agrees with P there).
    P: dict cell->0/1 covering D. Returns set.
    norm: optional map on cells for wraparound (e.g. column mod L); D must be normalised."""
    nm = norm or (lambda x: x)
    Dp = set()
    for (r, c) in D:
        for dr, dc in OFFS:
            Dp.add(nm((r + dr, c + dc)))
    vid = {cell: i + 1 for i, cell in enumerate(sorted(Dp))}
    var = lambda r, c: vid[nm((r, c))]
    s = Solver(name=solver)
    for (r, c) in D:
        for cl in rc.clauses(var, r, c, P[(r, c)]):
            s.add_clause(cl)
    cand = set(D)
    nxt = len(vid) + 1
    while cand:
        act = nxt; nxt += 1
        # some candidate differs from P
        s.add_clause([-act] + [(-vid[x] if P[x] else vid[x]) for x in cand])
        if not s.solve(assumptions=[act]):
            s.add_clause([-act])
            break
        m = s.get_model()
        mv = set(l for l in m if l > 0)
        diff = {x for x in cand if (vid[x] in mv) != bool(P[x])}
        cand -= diff
        s.add_clause([-act])
    s.delete()
    return cand


def max_self_forcing(rule_or_rc, P, D, log=print, norm=None):
    rc = rule_or_rc if isinstance(rule_or_rc, RuleCNF) else RuleCNF(rule_or_rc)
    D = set(D)
    it = 0
    while D:
        F = forced_cells(rc, P, D, norm=norm)
        it += 1
        if log:
            log(f"  iter {it}: |D|={len(D)} forced={len(F & D)}")
        if F >= D:
            return D
        D = D & F
    return D


def is_self_forcing(rc, P, D):
    return forced_cells(rc, P, set(D)) >= set(D)
