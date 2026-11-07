#!/usr/bin/env python3
"""
Coverage-lemma search, version 2.

New options (after the first genuine counterexample to the radius-1 rule):
  --history H   require the time-t pattern to have H generations of past (every oscillator phase has an
                infinite past, so this is always valid; it removes Garden-of-Eden-style counterexamples)
  --radius R    charge-splitting radius (1 = original rule; 2 = negative windows offer to positives within
                distance 2, splitting equally among them)
  --type N      run only positive type #N (e.g. rerun the counterexample with new options)

Soundness: for any fixed (R, H), if every positive type has worst-case margin >= 0, then every positive window in
every oscillator is covered, every window ends <= 0, and density <= 1/2 for every period.  Offers are floored and
r(p) is ceiled, so the check is conservative; negative margins are re-checked exactly.

Usage examples
  python coverage_search2.py --type 0 --history 2            # does the counterexample survive 2 steps of history?
  python coverage_search2.py --type 0 --radius 2 --history 1
  python coverage_search2.py --radius 2 --history 2 --workers 8 --time 900
"""
# v3: tables and positive-matchers are pre-filtered using the centre window's fixed cells (large memory saving).
import argparse, pickle, time
from fractions import Fraction
import numpy as np
from ortools.sat.python import cp_model
import life_agar_density as L

SCALE = 10**6

def charges(cert_path):
    W = L.reachable_windows(); A, b, c = L.build_lp(W)
    obj = pickle.load(open(cert_path, "rb")); u = obj["dual"] if isinstance(obj, dict) else obj
    At = A[:-1].T.tocsr(); r = []
    for w in range(len(W)):
        s = Fraction(int(c[w])) - Fraction(1, 2)
        for k in range(At.indptr[w], At.indptr[w+1]): s -= int(At.data[k]) * u[At.indices[k]]
        r.append(s)
    return W, r

def win_bits(w):
    bits = []
    for s in w: bits += [int(s) & 1, (int(s) >> 1) & 1]
    return bits

def life_link(m, src, dst, rad):
    """dst[(i,j)] = Life(src) for |i|,|j| <= rad (dst vars must already exist)"""
    for i in range(-rad, rad+1):
        for j in range(-rad, rad+1):
            n = m.NewIntVar(0, 8, "")
            m.Add(n == sum(src[(i+a, j+b)] for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)))
            b2, b3, z = m.NewBoolVar(""), m.NewBoolVar(""), m.NewBoolVar(""); x = src[(i, j)]; y = dst[(i, j)]
            m.Add(n == 2).OnlyEnforceIf(b2); m.Add(n != 2).OnlyEnforceIf(b2.Not())
            m.Add(n == 3).OnlyEnforceIf(b3); m.Add(n != 3).OnlyEnforceIf(b3.Not())
            m.AddBoolAnd([x, b2]).OnlyEnforceIf(z); m.AddBoolOr([x.Not(), b2.Not()]).OnlyEnforceIf(z.Not())
            m.AddBoolOr([b3, z]).OnlyEnforceIf(y); m.AddBoolAnd([b3.Not(), z.Not()]).OnlyEnforceIf(y.Not())

def filtered(rows, a, b, wp):
    """keep only table rows compatible with the centre window's fixed cells on the overlap"""
    mask = np.ones(len(rows), bool)
    for q in range(9):
        i, j = a + q // 3 - 1, b + q % 3 - 1
        if abs(i) <= 1 and abs(j) <= 1:
            s = int(wp[(i+1)*3 + (j+1)])
            mask &= (rows[:, 2*q] == (s & 1)) & (rows[:, 2*q+1] == ((s >> 1) & 1))
    return rows[mask]

def solve_type(p_idx, W, r, table, posrows, R, H, workers, tlimit, decide=False):
    m = cp_model.CpModel()
    C = 2*R + 1                                   # cells needed at t and t+1
    layers = []                                   # layers[h] = time t-h, radius C+1+h
    for h in range(H+1):
        rad = C + 1 + h
        layers.append({(i, j): m.NewBoolVar("") for i in range(-rad, rad+1) for j in range(-rad, rad+1)})
    X = layers[0]
    Y = {(i, j): m.NewBoolVar("") for i in range(-C, C+1) for j in range(-C, C+1)}
    life_link(m, X, Y, C)
    for h in range(1, H+1): life_link(m, layers[h], layers[h-1], C + h)
    bits_at = lambda a, b: sum(([X[(a+di-1, b+dj-1)], Y[(a+di-1, b+dj-1)]] for di in range(3) for dj in range(3)), [])
    for lit, val in zip(bits_at(0, 0), win_bits(W[p_idx])): m.Add(lit == val)
    maxslack = int(table[:, 18].max())
    slack, pos = {}, {(0, 0): 1}
    for a in range(-2*R, 2*R+1):
        for b in range(-2*R, 2*R+1):
            if (a, b) == (0, 0): continue
            bits = bits_at(a, b); pv = m.NewBoolVar("")
            if max(abs(a), abs(b)) <= R:                      # needs its slack: full table
                s = m.NewIntVar(0, maxslack, ""); m.AddAllowedAssignments(bits + [s, pv], filtered(table, a, b, W[p_idx]).tolist()); slack[(a, b)] = s
            else:                                             # only needs "is positive"
                ms = []
                for row in filtered(posrows, a, b, W[p_idx]).tolist():
                    mk = m.NewBoolVar(""); lits = [l if v else l.Not() for l, v in zip(bits, row)]
                    m.AddBoolAnd(lits).OnlyEnforceIf(mk); m.AddBoolOr([l.Not() for l in lits]).OnlyEnforceIf(mk.Not())
                    ms.append(mk)
                m.AddBoolOr(ms).OnlyEnforceIf(pv); m.AddBoolAnd([x.Not() for x in ms]).OnlyEnforceIf(pv.Not())
            pos[(a, b)] = pv
    offers = []
    for (a, b), s in slack.items():
        kq = m.NewIntVar(1, (2*R+1)**2, "")
        m.Add(kq == sum(pos[(a+da, b+db)] for da in range(-R, R+1) for db in range(-R, R+1)))
        off = m.NewIntVar(0, maxslack, ""); m.AddDivisionEquality(off, s, kq); offers.append(off)
    rp = -((-r[p_idx] * SCALE) // 1)
    margin = m.NewIntVar(-10**9, 10**9, ""); m.Add(margin == sum(offers) - int(rp))
    if decide: m.Add(margin <= -1)          # decision mode: only ask whether a failure exists
    else: m.Minimize(margin)
    solver = cp_model.CpSolver(); solver.parameters.num_search_workers = workers; solver.parameters.max_time_in_seconds = tlimit
    st = solver.Solve(m); pats = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        pats = [[[solver.Value(lay[(i, j)]) for j in range(-(C+1+h), C+2+h)] for i in range(-(C+1+h), C+2+h)]
                for h, lay in enumerate(layers)]
    return st, solver, margin, pats

def exact_recheck(Xt, W, r, R):
    P = np.array(Xt); c0 = P.shape[0] // 2; C = 2*R + 1
    def nxt(i, j):
        n = sum(P[i+a, j+b] for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0))
        return 1 if (n == 3 or (P[i, j] and n == 2)) else 0
    idx = {int(k): n for n, k in enumerate((W * (4 ** np.arange(9))).sum(1))}
    def rw(a, b):
        code = 0
        for q in range(9):
            di, dj = divmod(q, 3); i, j = c0+a+di-1, c0+b+dj-1
            code += (int(P[i, j]) | (nxt(i, j) << 1)) * 4 ** q
        return r[idx[code]]
    Rv = {(a, b): rw(a, b) for a in range(-2*R, 2*R+1) for b in range(-2*R, 2*R+1)}
    got = Fraction(0)
    for a in range(-R, R+1):
        for b in range(-R, R+1):
            if (a, b) == (0, 0) or Rv[(a, b)] >= 0: continue
            k = sum(1 for da in range(-R, R+1) for db in range(-R, R+1) if Rv[(a+da, b+db)] > 0)
            got += (-Rv[(a, b)]) / k
    return got, Rv[(0, 0)]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cert", default="certificate.pkl"); ap.add_argument("--radius", type=int, default=1)
    ap.add_argument("--history", type=int, default=0); ap.add_argument("--type", type=int, default=None)
    ap.add_argument("--start", type=int, default=0); ap.add_argument("--end", type=int, default=10**9)
    ap.add_argument("--decide", action="store_true", help="faster yes/no search (INFEASIBLE = covered)")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--time", type=float, default=900)
    a = ap.parse_args(); t0 = time.time()
    W, r = charges(a.cert)
    positive = [i for i, rv in enumerate(r) if rv > 0]
    table = [win_bits(w) + [int(((-rv) * SCALE) // 1) if rv < 0 else 0, 1 if rv > 0 else 0] for w, rv in zip(W, r)]
    posrows = np.array([win_bits(W[i]) for i in positive], dtype=np.int64)
    table = np.array(table, dtype=np.int64)
    print(f"[{time.time()-t0:.0f}s] positive types {len(positive)}; radius {a.radius}; history {a.history}", flush=True)
    todo = [(a.type, positive[a.type])] if a.type is not None else list(enumerate(positive))[a.start:a.end]
    out = open(f"coverage_R{a.radius}_H{a.history}.txt", "a")
    for n, p in todo:
        st, solver, margin, pats = solve_type(p, W, r, table, posrows, a.radius, a.history, a.workers, a.time, a.decide)
        if a.decide and st == cp_model.INFEASIBLE:
            msg = f"type #{n} (window {p}): COVERED (no failing configuration exists)"
        elif (not a.decide) and st == cp_model.OPTIMAL and solver.Value(margin) >= 0:
            msg = f"type #{n} (window {p}): COVERED, worst-case margin {solver.Value(margin)}/{SCALE}"
        elif pats is not None and solver.Value(margin) < 0:
            got, rp = exact_recheck(pats[0], W, r, a.radius)
            bad = got < rp
            msg = (f"type #{n} (window {p}): margin {solver.Value(margin)}/{SCALE} [{solver.StatusName(st)}] -> exact received "
                   f"{float(got):.6f} vs r {float(rp):.6f}: " + ("GENUINE COUNTEREXAMPLE" if bad else "rounding artefact"))
            if bad:
                with open(f"cex_R{a.radius}_H{a.history}_type{n}.txt", "w") as f:
                    for h, lay in enumerate(pats):
                        f.write(f"time t-{h}:\n" + "\n".join("".join("o" if v else "." for v in row) for row in lay) + "\n\n")
        else:
            msg = f"type #{n} (window {p}): UNDECIDED ({solver.StatusName(st)}; bound {solver.BestObjectiveBound()}/{SCALE})"
        print(f"[{time.time()-t0:.0f}s] " + msg, flush=True); out.write(msg + "\n"); out.flush()
