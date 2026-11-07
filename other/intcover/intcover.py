#!/usr/bin/env python3
"""Smallest pattern containing every INT transition (all 102) at least once.

A transition is a D4 orbit of (centre, 8 neighbours) -- 51 with a dead centre,
51 with a live one. A pattern contains it if some cell of the infinite plane
(dead outside the pattern) has that neighbourhood. B0 is therefore free.

Exact CP-SAT model in a W x H box, dead outside it; the 1-cell ring around the
box is checked too. Default objective: minimum population in that box. OPTIMAL
means optimal for that box only. --feasible just asks whether the box can work.

  python3 intcover.py 12 12 --time 300 --out best.rle
  python3 intcover.py 10 10 --feasible

Index convention: Golly MAP, NW 256, N 128, NE 64, W 32, C 16, E 8, SW 4, S 2,
SE 1; x right, y down.
"""
import argparse, itertools, os, sys, time

# Optional independent check with casim (a separate simulator); skipped if it isn't found.
CASIM = os.environ.get("CASIM", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "casim"))
OFFS = [(dx, dy) for dy in (-1, 0, 1) for dx in (-1, 0, 1)]   # MSB first
LETTERS = {0: "", 1: "ce", 2: "cekain", 3: "cekainyqjr", 4: "cekainyqjrtwz",
           5: "cekainyqjr", 6: "cekain", 7: "ce", 8: ""}


def orbit_classes():
    """class id per 512 indices, built from the square symmetries."""
    def grid(i): return [[(i >> (8 - 3 * r - c)) & 1 for c in range(3)] for r in range(3)]
    def idx(g): return sum(g[r][c] << (8 - 3 * r - c) for r in range(3) for c in range(3))
    def rot(g): return [[g[2 - c][r] for c in range(3)] for r in range(3)]
    def flip(g): return [row[::-1] for row in g]
    cls, n = [-1] * 512, 0
    for i in range(512):
        if cls[i] < 0:
            g = grid(i)
            for _ in range(4):
                g = rot(g)
                cls[idx(g)] = n
                cls[idx(flip(g))] = n
            n += 1
    assert n == 102
    return cls


def covered(cells, cls):
    live = set(cells)
    cand = {(x + dx, y + dy) for x, y in live for dx, dy in OFFS}
    seen = {cls[0]}                       # background
    for x, y in cand:
        i = 0
        for dx, dy in OFFS:
            i = 2 * i + ((x + dx, y + dy) in live)
        seen.add(cls[i])
    return seen


def to_rle(cells, w, h):
    live = set(cells)
    rows = []
    for y in range(h):
        row = "".join("o" if (x, y) in live else "b" for x in range(w)).rstrip("b")
        runs = [(k, len(list(g))) for k, g in itertools.groupby(row)]
        rows.append("".join(f"{n if n > 1 else ''}{k}" for k, n in runs))
    return f"x = {w}, y = {h}, rule = B3/S23\n" + "$".join(rows) + "!\n"


def casim_check(cells):
    """Independent check: Golly letter semantics via casim. For each transition
    T, the rule containing only T leaves something alive iff T occurs."""
    sys.path.insert(0, CASIM)
    try:
        import casim as C
    except ImportError:
        return None, None
    missing = []
    names = [("B", n, l) for n in range(1, 9) for l in (LETTERS[n] or " ")]
    names += [("S", n, l) for n in range(9) for l in (LETTERS[n] or " ")]
    for t, n, l in names:
        tr = f"{n}{l.strip()}"
        rule = f"B{tr}/S" if t == "B" else f"B/S{tr}"
        p = C.Pattern.from_coords(cells, rule)
        p.step(1)
        if p.population == 0:
            missing.append(t + tr)
    return len(names), missing


def solve(w, h, a, cls):
    from unittest import mock             # ortools imports pandas but never needs it
    sys.modules["pandas"] = mock.MagicMock()   # here, and the anaconda pandas is broken
    from ortools.sat.python import cp_model
    m = cp_model.CpModel()
    v = {(x, y): m.NewBoolVar(f"c{x}_{y}") for y in range(h) for x in range(w)}
    by_class = [[] for _ in range(102)]
    for y in range(-1, h + 1):
        for x in range(-1, w + 1):
            lits = [v.get((x + dx, y + dy)) for dx, dy in OFFS]
            one = []
            for i in range(512):
                bits = [(i >> (8 - k)) & 1 for k in range(9)]
                if any(b and l is None for b, l in zip(bits, lits)):
                    continue
                z = m.NewBoolVar("")
                one.append(z)
                by_class[cls[i]].append(z)
                for b, l in zip(bits, lits):
                    if l is not None:
                        m.AddImplication(z, l if b else l.Not())
            m.AddExactlyOne(one)
    for k in range(102):
        if k != cls[0]:
            m.AddBoolOr(by_class[k])
    pop = sum(v.values())
    if a.feasible:
        pass
    else:
        m.Add(pop >= 55)                  # window-count LP bound, see README
        m.Minimize(pop)
    if a.hint:                            # rows of . and o, e.g. the last block sa printed
        rows = [r.strip() for r in open(a.hint) if r.strip() and not r.startswith("pop")][-h:]
        for (x, y), b in v.items():
            m.AddHint(b, rows[y][x] == "o")
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = a.time
    s.parameters.num_workers = a.threads
    t0 = time.time()

    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            cells = [c for c, b in v.items() if self.Value(b)]
            print(f"[{time.time() - t0:7.1f}s] pop {len(cells)} bound {self.BestObjectiveBound():.0f}",
                  flush=True)
            if a.out:
                open(a.out, "w").write(to_rle(cells, w, h))

    st = s.Solve(m, CB())
    name = s.StatusName(st)
    cells = [c for c, b in v.items() if s.Value(b)] if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    return name, cells, s.BestObjectiveBound() if not a.feasible else None, time.time() - t0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("width", type=int)
    p.add_argument("height", type=int)
    p.add_argument("--time", type=float, default=120, help="seconds (default 120)")
    p.add_argument("--threads", type=int, default=10)
    p.add_argument("--feasible", action="store_true", help="stop at any solution, no objective")
    p.add_argument("--hint", help="warm start: file of ./o rows (sa output)")
    p.add_argument("--out", help="write best RLE here on each improvement")
    a = p.parse_args()
    cls = orbit_classes()
    status, cells, bound, t = solve(a.width, a.height, a, cls)
    print(f"{a.width}x{a.height}: {status} in {t:.1f}s" + (f", bound {bound:.0f}" if bound is not None else ""))
    if cells:
        seen = covered(cells, cls)
        total, missing = casim_check(cells)
        extra = (f"; casim {total - len(missing)}/{total} (B0 free)" + (f" MISSING {missing}" if missing else "")
                 if total is not None else " (casim not found, independent check skipped)")
        print(f"pop {len(cells)}; own check {len(seen)}/102" + extra)
        print(to_rle(cells, a.width, a.height), end="")


if __name__ == "__main__":
    main()
