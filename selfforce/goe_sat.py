"""Random-pattern orphan test: SAT-check whether random n x n patterns have a predecessor.
usage: goe_sat.py RULE n trials [--density 0.5] [--seed 0]"""
import argparse, random, time
import sfcore as S
from pysat.solvers import Solver

def has_pred(rc, P, cells=None):
    """P: dict cell->0/1 (any shape). Returns predecessor dict or None."""
    cells = cells or list(P)
    vid = {}
    def v(r, c):
        if (r, c) not in vid: vid[(r, c)] = len(vid) + 1
        return vid[(r, c)]
    s = Solver(name="cadical153")
    for (r, c) in cells:
        for cl in rc.clauses(v, r, c, P[(r, c)]):
            s.add_clause(cl)
    ok = s.solve()
    y = None
    if ok:
        m = set(l for l in s.get_model() if l > 0)
        y = {k: int(i in m) for k, i in vid.items()}
    s.delete()
    return y if ok else None

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("rule"); ap.add_argument("n", type=int); ap.add_argument("trials", type=int)
    ap.add_argument("--density", type=float, default=0.5); ap.add_argument("--seed", type=int, default=0)
    x = ap.parse_args()
    rc = S.RuleCNF(x.rule); random.seed(x.seed)
    orph = 0; t = time.time()
    for i in range(x.trials):
        P = {(r, c): int(random.random() < x.density) for r in range(x.n) for c in range(x.n)}
        if has_pred(rc, P) is None:
            orph += 1
            if orph == 1:
                print("first orphan:"); print("\n".join("".join("o" if P[(r, c)] else "." for c in range(x.n)) for r in range(x.n)))
    print(f"{x.rule} n={x.n} d={x.density}: {orph}/{x.trials} orphans, {time.time()-t:.1f}s")
