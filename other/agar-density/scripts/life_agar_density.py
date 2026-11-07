#!/usr/bin/env python3
"""
Oscillating-agar density in Conway's Game of Life -- reproduction + remaining computations.

Requirements:  pip install numpy scipy highspy python-sat python-flint

Subcommands
  certify   Rebuild the all-period LP, run the "peeling" (find the culprit window families), and produce an EXACT
            rational certificate that every oscillator avoiding the excluded window types has density <= 1/2.
  embed     Exact maximum density of REAL oscillators (period p, torus width x L rows) that contain a given band
            pattern in two consecutive phases.  Use this to test the two remaining families.
  patch     Largest rectangles that can be filled using only the windows of the LP's current extremal solution.

Examples
  python life_agar_density.py certify
  python life_agar_density.py embed --family 1 --period 2 --width 6  --rows 6-30
  python life_agar_density.py embed --family 2 --period 2 --width 4  --rows 6-30
  python life_agar_density.py embed --family 2 --period 3 --width 6  --rows 6-24

Conventions: each cell of a two-generation window has state  0 '.' dead->dead, 1 'v' alive->dead,
2 '^' dead->alive, 3 'o' alive->alive.  Windows are 3x3, row-major.
"""
import argparse, itertools, os, pickle, sys, time
from fractions import Fraction
import numpy as np

T0 = time.time()
def log(*a): print(f"[{time.time()-T0:7.0f}s]", *a, flush=True)
def life(c, n): return 1 if (n == 3 or (c and n == 2)) else 0

# ----------------------------------------------------------------------------------------------------------------
# 1. Reachable two-generation 3x3 windows (3x3 at t, 3x3 at t+1) produced by some 5x5 pattern at t
# ----------------------------------------------------------------------------------------------------------------
def reachable_windows(cache="reach.npy"):
    if os.path.exists(cache): reach = np.load(cache)
    else:
        log("enumerating all 2^25 5x5 patterns ...")
        reach = np.zeros(1 << 18, bool)
        for s in range(0, 1 << 25, 1 << 22):
            P = np.arange(s, s + (1 << 22), dtype=np.int64)
            core = np.zeros_like(P); nxt = np.zeros_like(P)
            for k in range(9):
                r, c = divmod(k, 3); r += 1; c += 1
                cc = (P >> (r*5+c)) & 1
                n = sum((P >> ((r+a)*5+(c+b))) & 1 for a in (-1,0,1) for b in (-1,0,1) if (a,b) != (0,0))
                core |= cc << k; nxt |= ((n == 3) | ((cc == 1) & (n == 2))).astype(np.int64) << k
            reach[core | (nxt << 9)] = True
        np.save(cache, reach)
    codes = np.nonzero(reach)[0]
    core, nxt = codes & 511, codes >> 9
    W = np.stack([((core >> k) & 1) | (((nxt >> k) & 1) << 1) for k in range(9)], 1)
    log("reachable windows:", len(W))
    return W

# ----------------------------------------------------------------------------------------------------------------
# 2. The all-period LP:  spatial consistency (3x2 / 2x3 overlaps), stationarity (layer-0 3x3 marginal = layer-1),
#    normalisation.  Objective: frequency of a live centre at time t.
# ----------------------------------------------------------------------------------------------------------------
def build_lp(W):
    from scipy.sparse import coo_matrix
    N = len(W)
    def key(cols, layer=None):
        v = np.zeros(N, np.int64)
        for k in cols: v = v*4 + (W[:, k] if layer is None else (W[:, k] >> layer) & 1)
        return v
    pairs = [(key([0,1,3,4,6,7]), key([1,2,4,5,7,8])), (key([0,1,2,3,4,5]), key([3,4,5,6,7,8])),
             (key(range(9), 0), key(range(9), 1))]
    rows, cols, vals, off = [], [], [], 0
    for a, b in pairs:
        u, inv = np.unique(np.concatenate([a, b]), return_inverse=True)
        rows += [off+inv[:N], off+inv[N:]]; cols += [np.arange(N)]*2; vals += [np.ones(N), -np.ones(N)]; off += len(u)
    rows.append(np.full(N, off)); cols.append(np.arange(N)); vals.append(np.ones(N))
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(off+1, N)).tocsr()
    b = np.zeros(off+1); b[-1] = 1
    c = (W[:, 4] & 1).astype(float)
    return A, b, c

# culprit windows: centre is a cell of a vanishing 2x2 block (>=3 'v', rest '.', centre 'v') and has a surviving
# ('o') orthogonal neighbour outside that block
def culprit_mask(W):
    outer = {(0,1,3,4): (5,7), (1,2,4,5): (3,7), (3,4,6,7): (1,5), (4,5,7,8): (1,3)}
    m = np.zeros(len(W), bool)
    for blk, out in outer.items():
        st = W[:, list(blk)]
        van = np.isin(st, [0,1]).all(1) & ((st == 1).sum(1) >= 3) & (W[:, 4] == 1)
        m |= van & ((W[:, out[0]] == 3) | (W[:, out[1]] == 3))
    return m

def solve_lp(A, b, c, keep, tlim=600):
    from scipy.optimize import linprog
    res = linprog(-c[keep], A_eq=A[:, keep], b_eq=b, bounds=(0, None), method='highs-ipm',
                  options={'time_limit': tlim})
    return res

# ----------------------------------------------------------------------------------------------------------------
# 3. Exact certificate: HiGHS IPM + crossover -> basis -> exact rational dual (FLINT) -> exact check of all windows
# ----------------------------------------------------------------------------------------------------------------
def exact_certificate(A, b, c, keep, out="certificate.pkl"):
    import highspy, flint
    A2 = A[:, keep].tocsc(); c2 = c[keep]; nw, nr = A2.shape[1], A2.shape[0]
    h = highspy.Highs(); h.setOptionValue('output_flag', False)
    h.setOptionValue('solver', 'ipm'); h.setOptionValue('run_crossover', 'on')
    lp = highspy.HighsLp(); lp.num_col_ = nw; lp.num_row_ = nr
    inf = highspy.kHighsInf
    lp.col_cost_ = (-c2).astype(float); lp.col_lower_ = np.zeros(nw); lp.col_upper_ = np.full(nw, inf)
    lp.row_lower_ = b.astype(float); lp.row_upper_ = b.astype(float)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A2.indptr; lp.a_matrix_.index_ = A2.indices; lp.a_matrix_.value_ = A2.data.astype(float)
    h.passModel(lp); h.run()
    basis = h.getBasis()
    Bc = np.array([s == highspy.HighsBasisStatus.kBasic for s in basis.col_status])
    Br = np.array([s == highspy.HighsBasisStatus.kBasic for s in basis.row_status])
    colsB, rowsN = np.nonzero(Bc)[0], np.nonzero(~Br)[0]
    log(f"basis valid={basis.valid}; exact system {len(colsB)} x {len(rowsN)}")
    sub = A2[rowsN][:, colsB].T.toarray().astype(np.int64)
    M = flint.fmpq_mat(len(colsB), len(rowsN), [int(v) for v in sub.ravel()])
    rhs = flint.fmpq_mat(len(colsB), 1, [-int(c2[w]) for w in colsB])
    Y = M.solve(rhs); log("solved exactly")
    y = [Fraction(0)]*nr
    for i, r in enumerate(rowsN): e = Y[i, 0]; y[r] = Fraction(int(e.p), int(e.q))
    u = [-v for v in y]
    worst = Fraction(-10**9)
    for w in range(nw):                       # exact check:  c_w - sum_{consistency rows} A_rw u_r  <= bound
        s = Fraction(int(c2[w]))
        for k in range(A2.indptr[w], A2.indptr[w+1]):
            r = A2.indices[k]
            if r != nr-1: s -= int(A2.data[k]) * u[r]
        worst = max(worst, s)
    pickle.dump({'dual': u, 'keep': keep}, open(out, 'wb'))
    log(f"EXACT BOUND = {worst}  (certificate saved to {out})")
    return worst

# ----------------------------------------------------------------------------------------------------------------
# 4. Patch filling with a given window set (SAT)
# ----------------------------------------------------------------------------------------------------------------
def rect_fillable(S, h, w, budget=400000):
    from pysat.solvers import Cadical153
    from pysat.formula import IDPool
    pool = IDPool(); v = lambda bt, i, j: pool.id((bt, i, j)); cls = []
    for i in range(1, h-1):
        for j in range(1, w-1):
            sel = [pool.id(('s', i, j, k)) for k in range(len(S))]; cls.append(sel)
            for k, wd in enumerate(S):
                for q in range(9):
                    r, cc = divmod(q, 3)
                    for bit in (0, 1):
                        lit = v(bit, i+r-1, j+cc-1); cls.append([-sel[k], lit if (wd[q] >> bit) & 1 else -lit])
    s = Cadical153(bootstrap_with=cls); s.conf_budget(budget); return s.solve_limited()

# ----------------------------------------------------------------------------------------------------------------
# 5. Embedding a band into REAL oscillators and maximising density
# ----------------------------------------------------------------------------------------------------------------
FAMILY = {1: ["vvo^^o", "vv.^^.", ".oo.oo", "^^.vv.", "^^ovvo"],        # blinking-block band (period 6)
          2: ["^.", "^^", "vv", "vv", "^^", "^."]}                        # moving double-stripe variant (period 2)

def embed_max(band, p, wd, L, budget=600):
    from pysat.solvers import Cadical153
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    TPL = [(bt, life(bt[0], sum(bt[1:]))) for bt in itertools.product((0,1), repeat=9)]
    pool = IDPool(); v = lambda t, i, j: pool.id((t % p, i % L, j % wd)); cls = []
    for t in range(p):
        for i in range(L):
            for j in range(wd):
                lits = [v(t,i,j)] + [v(t,i+a,j+b) for a in (-1,0,1) for b in (-1,0,1) if (a,b) != (0,0)]
                for bits, r in TPL:
                    cls.append([-l if bb else l for l, bb in zip(lits, bits)] + [v(t+1,i,j) if r else -v(t+1,i,j)])
    per = len(band[0])
    for i, row in enumerate(band):
        for j in range(wd):
            s = '.v^o'.index(row[j % per])
            cls += [[v(0,i,j)] if s & 1 else [-v(0,i,j)], [v(1,i,j)] if s & 2 else [-v(1,i,j)]]
    allv = [v(t,i,j) for t in range(p) for i in range(L) for j in range(wd)]
    k, best, t0 = 0, None, time.time()
    while time.time()-t0 < budget:
        S = Cadical153(bootstrap_with=cls + CardEnc.atleast(lits=allv, bound=k, vpool=pool,
                                                            encoding=EncType.seqcounter).clauses)
        r = S.solve()
        if not r: return best, "optimal"
        mod = set(x for x in S.get_model() if x > 0); best = sum(1 for x in allv if x in mod); k = best + 1
    return best, "time limit (lower bound only)"

# ----------------------------------------------------------------------------------------------------------------
def cmd_certify(args):
    W = reachable_windows(); A, b, c = build_lp(W); N = len(W)
    cul = culprit_mask(W); removed = np.zeros(N, bool)
    for it in range(20):
        keep = ~removed; res = solve_lp(A, b, c, keep)
        u = -res.eqlin.marginals; rc = (A.T @ u) - c; rc[removed] = np.inf
        tight = np.nonzero(rc < 1e-7)[0]; tc = tight[cul[tight]]
        log(f"peel round {it}: LP bound {-res.fun:.6f}; tight culprit windows {len(tc)}")
        if -res.fun < 0.5 + 1e-9 or len(tc) == 0: break
        removed[tc] = True
    pickle.dump(removed, open("excluded_windows.pkl", "wb"))
    log(f"excluded {int(removed.sum())} window types (saved to excluded_windows.pkl); certifying exactly ...")
    exact_certificate(A, b, c, ~removed)

def cmd_embed(args):
    band = FAMILY[args.family]; lo, hi = map(int, args.rows.split('-'))
    for L in range(max(lo, len(band)+1), hi+1):
        best, note = embed_max(band, args.period, args.width, L, args.budget)
        tot = args.period*args.width*L
        msg = "no oscillator contains it" if best is None else f"best density {best}/{tot} = {best/tot:.4f}"
        log(f"family {args.family}, period {args.period}, width {args.width}, rows {L}: {msg} ({note})")
        if best is not None and best*2 > tot: log("   *** DENSITY ABOVE 1/2 -- CHECK THIS CASE ***")

def cmd_patch(args):
    W = reachable_windows(); A, b, c = build_lp(W); keep = np.ones(len(W), bool)
    if os.path.exists("excluded_windows.pkl") and args.peeled:
        keep = ~pickle.load(open("excluded_windows.pkl", "rb"))
    res = solve_lp(A, b, c, keep); xs = np.zeros(len(W)); xs[keep] = res.x
    S = W[xs > 1e-9]; log(f"LP bound {-res.fun:.6f}; support {len(S)} windows")
    for h in range(3, 10):
        log(f"height {h}: " + ' '.join(f"w{w}:{ {True:'Y',False:'-',None:'?'}[rect_fillable(S,h,w)] }"
                                          for w in (6, 10, 14, 20)))

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); sp = ap.add_subparsers(dest="cmd", required=True)
    sp.add_parser("certify")
    e = sp.add_parser("embed"); e.add_argument("--family", type=int, default=1); e.add_argument("--period", type=int, default=2)
    e.add_argument("--width", type=int, default=6); e.add_argument("--rows", default="6-24"); e.add_argument("--budget", type=int, default=600)
    p = sp.add_parser("patch"); p.add_argument("--peeled", action="store_true")
    a = ap.parse_args()
    {"certify": cmd_certify, "embed": cmd_embed, "patch": cmd_patch}[a.cmd](a)
