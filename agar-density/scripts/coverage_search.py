#!/usr/bin/env python3
"""
Search for counterexamples to the COVERAGE LEMMA (the last step of the oscillating-agar density proof).

Setting
  Every two-generation 3x3 window w (cells at time t and t+1) has an exact charge r(w) from the LP certificate:
      r(w) = [centre alive at t] - sum_rows A[row,w]*u[row] - 1/2          (u = certificate dual)
  For every oscillator, density - 1/2 = average of r over all windows and phases.  Allowed windows have r <= 0;
  some "family" windows have r > 0.

Rule (radius 1)
  Every window q with r(q) < 0 splits |r(q)| equally among the positive windows within king distance 1 of q
  (same generation pair).  Negative windows never give more than they have.

Coverage lemma
  In every two-generation configuration, every positive window p receives >= r(p).
  => every window ends <= 0 => density <= 1/2 for EVERY period.

What this script checks
  For each positive window type p (placed at the centre), it asks CP-SAT to MINIMISE
        sum_q floor(slack(q) / k(q))  -  ceil(r(p))          (all charges scaled by SCALE)
  over ALL 9x9 patterns at time t (time t+1 follows from Life).  k(q) = number of positive windows within
  distance 1 of q.  Floors on offers and the ceiling on r(p) make this conservative: if the minimum is >= 0 for
  every type, the lemma is PROVEN (modulo code correctness).  A negative minimum is a candidate counterexample;
  it is then re-checked in exact rational arithmetic.

Requirements:  pip install ortools numpy scipy
               plus life_agar_density.py in the same folder, and its outputs from `certify`
               (reach.npy, certificate.pkl, excluded_windows.pkl).

Usage
  python coverage_search.py                         # all positive types
  python coverage_search.py --start 0 --end 100     # a slice (run several slices in parallel)
  python coverage_search.py --workers 8 --time 600  # CP-SAT workers / time limit per type (seconds)
Results are appended to coverage_results.txt; counterexamples are saved as counterexample_<type>.txt.
"""
import argparse, pickle, time, sys
from fractions import Fraction
import numpy as np
from ortools.sat.python import cp_model
import life_agar_density as L

SCALE = 10**6

def charges(cert_path):
    W = L.reachable_windows()
    A, b, c = L.build_lp(W)
    obj = pickle.load(open(cert_path, "rb"))
    u = obj["dual"] if isinstance(obj, dict) else obj          # list of Fractions, one per LP row
    nr = A.shape[0]
    At = A[:-1].T.tocsr()                                        # exclude normalisation row
    r = []
    for w in range(len(W)):
        s = Fraction(int(c[w])) - Fraction(1, 2)
        for k in range(At.indptr[w], At.indptr[w+1]):
            s -= int(At.data[k]) * u[At.indices[k]]
        r.append(s)
    return W, r

def build_table(W, r):
    """rows: 18 bits (x,y per cell, row-major) + slack_int + positive flag"""
    rows = []
    for w, rw in zip(W, r):
        bits = []
        for s in w: bits += [int(s) & 1, (int(s) >> 1) & 1]
        slack = int(((-rw) * SCALE) // 1) if rw < 0 else 0      # floor(|r|*SCALE): conservative
        rows.append(bits + [slack, 1 if rw > 0 else 0])
    return rows

def solve_type(p_idx, W, r, table, workers, tlimit):
    m = cp_model.CpModel()
    X = {(i, j): m.NewBoolVar(f"x{i}_{j}") for i in range(-4, 5) for j in range(-4, 5)}      # time t, 9x9
    Y = {}
    for i in range(-3, 4):                                                                   # time t+1, 7x7
        for j in range(-3, 4):
            n = m.NewIntVar(0, 8, f"n{i}_{j}")
            m.Add(n == sum(X[(i+a, j+b)] for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)))
            b2, b3, z, y = (m.NewBoolVar(f"{t}{i}_{j}") for t in ("b2", "b3", "z", "y"))
            m.Add(n == 2).OnlyEnforceIf(b2); m.Add(n != 2).OnlyEnforceIf(b2.Not())
            m.Add(n == 3).OnlyEnforceIf(b3); m.Add(n != 3).OnlyEnforceIf(b3.Not())
            m.AddBoolAnd([X[(i, j)], b2]).OnlyEnforceIf(z); m.AddBoolOr([X[(i, j)].Not(), b2.Not()]).OnlyEnforceIf(z.Not())
            m.AddBoolOr([b3, z]).OnlyEnforceIf(y); m.AddBoolAnd([b3.Not(), z.Not()]).OnlyEnforceIf(y.Not())
            Y[(i, j)] = y
    # centre window fixed to type p
    wp = W[p_idx]
    for q in range(9):
        di, dj = divmod(q, 3); di -= 1; dj -= 1
        m.Add(X[(di, dj)] == (int(wp[q]) & 1)); m.Add(Y[(di, dj)] == ((int(wp[q]) >> 1) & 1))
    maxslack = max(row[18] for row in table)
    slack, pos = {}, {}
    for a in range(-2, 3):
        for b in range(-2, 3):
            if (a, b) == (0, 0): pos[(a, b)] = 1; continue
            bits = []
            for q in range(9):
                di, dj = divmod(q, 3); di += a - 1; dj += b - 1
                bits += [X[(di, dj)], Y[(di, dj)]]
            s = m.NewIntVar(0, maxslack, f"s{a}_{b}"); pv = m.NewBoolVar(f"p{a}_{b}")
            m.AddAllowedAssignments(bits + [s, pv], table)
            slack[(a, b)] = s; pos[(a, b)] = pv
    offers = []
    for a in range(-1, 2):
        for b in range(-1, 2):
            if (a, b) == (0, 0): continue
            kq = m.NewIntVar(1, 9, f"k{a}_{b}")
            m.Add(kq == sum(pos[(a+da, b+db)] for da in (-1, 0, 1) for db in (-1, 0, 1)))
            off = m.NewIntVar(0, maxslack, f"o{a}_{b}")
            m.AddDivisionEquality(off, slack[(a, b)], kq)            # floor division: conservative
            offers.append(off)
    rp = -((-r[p_idx] * SCALE) // 1)                                  # ceil(r(p)*SCALE)
    margin = m.NewIntVar(-10**9, 10**9, "margin")
    m.Add(margin == sum(offers) - int(rp))
    m.Minimize(margin)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.max_time_in_seconds = tlimit
    st = solver.Solve(m)
    pattern = None
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        pattern = [[solver.Value(X[(i, j)]) for j in range(-4, 5)] for i in range(-4, 5)]
    return st, solver, margin, pattern

def exact_recheck(pattern, W, r):
    """Exact received charge at the centre for a 9x9 pattern (list of lists, time t)."""
    P = np.array(pattern)
    def nxt(i, j):
        n = sum(P[i+a, j+b] for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0))
        return 1 if (n == 3 or (P[i, j] and n == 2)) else 0
    Yv = {(i, j): nxt(i, j) for i in range(1, 8) for j in range(1, 8)}
    idx = {int(k): n for n, k in enumerate((W * (4 ** np.arange(9))).sum(1))}
    def win(ci, cj):
        code = 0
        for q in range(9):
            di, dj = divmod(q, 3)
            i, j = ci + di - 1, cj + dj - 1
            code += (int(P[i, j]) | (Yv[(i, j)] << 1)) * 4 ** q
        return idx[code]
    R = {(a, b): r[win(4+a, 4+b)] for a in range(-2, 3) for b in range(-2, 3)}
    got = Fraction(0)
    for a in range(-1, 2):
        for b in range(-1, 2):
            if (a, b) == (0, 0) or R[(a, b)] >= 0: continue
            k = sum(1 for da in (-1, 0, 1) for db in (-1, 0, 1) if R[(a+da, b+db)] > 0)
            got += (-R[(a, b)]) / k
    return got, R[(0, 0)]

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cert", default="certificate.pkl")
    ap.add_argument("--start", type=int, default=0); ap.add_argument("--end", type=int, default=10**9)
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--time", type=float, default=600)
    a = ap.parse_args()
    t0 = time.time()
    W, r = charges(a.cert)
    positive = [i for i, rv in enumerate(r) if rv > 0]
    print(f"[{time.time()-t0:.0f}s] windows {len(W)}, positive types {len(positive)}", flush=True)
    table = build_table(W, r)
    out = open("coverage_results.txt", "a")
    for n, p in enumerate(positive[a.start:a.end], start=a.start):
        st, solver, margin, pattern = solve_type(p, W, r, table, a.workers, a.time)
        name = solver.StatusName(st)
        if st == cp_model.OPTIMAL and solver.Value(margin) >= 0:
            msg = f"type #{n} (window {p}): COVERED, worst-case margin {solver.Value(margin)}/{SCALE}"
        elif st in (cp_model.OPTIMAL, cp_model.FEASIBLE) and solver.Value(margin) < 0:
            got, rp = exact_recheck(pattern, W, r)
            verdict = "GENUINE COUNTEREXAMPLE (rule fails here)" if got < rp else "rounding artefact (exactly covered)"
            msg = (f"type #{n} (window {p}): margin {solver.Value(margin)}/{SCALE} [{name}] -> exact: received {float(got):.6f} "
                   f"vs r {float(rp):.6f}: {verdict}")
            with open(f"counterexample_{n}.txt", "w") as f:
                f.write("\n".join("".join("o" if v else "." for v in row) for row in pattern) + "\n")
        else:
            msg = f"type #{n} (window {p}): UNDECIDED ({name}; best bound {solver.BestObjectiveBound()}/{SCALE}) - rerun with more time"
        print(f"[{time.time()-t0:.0f}s] " + msg, flush=True); out.write(msg + "\n"); out.flush()
