#!/usr/bin/env python3
"""
Larger-window all-period LP for the oscillating-agar density problem (Conway's Life).

Pipeline for an h x w window (two generations per cell):
  1. enumerate   all genuinely occurring two-generation h x w windows ((h+2) x (w+2) patterns at time t)
  2. solve       the all-period LP (spatial consistency + stationarity + normalisation); prints a RIGOROUS bound
                 and saves the duals, so certify never has to re-solve
     --sym       solve the reflection-symmetric LP instead (one column per mirror-orbit, up to 4x smaller, same value)
  3. peel        repeatedly exclude culprit windows that are tight at the optimum, until the bound reaches 1/2;
                 the final round's duals are saved and certified automatically
  4. certify     rigorous bound from SAVED duals (no re-solve), exact integer check -> "density <= 1/2 + eps"
  5. export      write the (peeled) LP as MPS for an exact rational solver (SoPlex) -> exactly 1/2

Requirements: pip install numpy scipy highspy
Usage:
  python big_window_lp.py enumerate --h 3 --w 4
  python big_window_lp.py solve     --h 3 --w 4 --threads 16
  python big_window_lp.py peel      --h 3 --w 4 --threads 16
  python big_window_lp.py certify   --h 3 --w 4 --peeled
  python big_window_lp.py export    --h 3 --w 4 --peeled          # -> lp_3x4_peeled.mps
Cell states: 0 dead->dead, 1 alive->dead, 2 dead->alive, 3 alive->alive.  Windows are row-major.
"""
import argparse, os, time, pickle
import numpy as np
from scipy.sparse import coo_matrix

T0 = time.time()
def log(*a): print(f"[{time.time()-T0:8.0f}s]", *a, flush=True)

# ---------------------------------------------------------------------------------------------------- 1. enumerate
def _chunk(args):
    s, e, h, w = args; Wd = w + 2; nc = h * w
    P = np.arange(s, e, dtype=np.int64); core = np.zeros_like(P); nxt = np.zeros_like(P)
    for k in range(nc):
        r, c = divmod(k, w); r += 1; c += 1
        cc = (P >> (r*Wd + c)) & 1
        n = sum((P >> ((r+a)*Wd + (c+b))) & 1 for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0))
        core |= cc << k; nxt |= ((n == 3) | ((cc == 1) & (n == 2))).astype(np.int64) << k
    return np.unique(core | (nxt << nc))

def enumerate_windows(h, w, chunk_bits=22, jobs=1):
    out = f"windows_{h}x{w}.npy"
    if os.path.exists(out): return np.load(out)
    nbits = (h + 2) * (w + 2); nc = h * w
    tasks = [(s, min(s + (1 << chunk_bits), 1 << nbits), h, w) for s in range(0, 1 << nbits, 1 << chunk_bits)]
    codes = np.zeros(0, np.int64); pending = []
    from multiprocessing import Pool
    with Pool(jobs) as pool:                                   # parallel over chunks
        for i, part in enumerate(pool.imap_unordered(_chunk, tasks, chunksize=4)):
            pending.append(part)
            if len(pending) >= 32:                             # merge in batches (one sort per 32 chunks)
                codes = np.unique(np.concatenate([codes] + pending)); pending = []
            if i % 64 == 0: log(f"enumerate {h}x{w}: {i/len(tasks):6.1%} done, {len(codes)} windows so far")
    codes = np.unique(np.concatenate([codes] + pending))
    core, nxt = codes & ((1 << nc) - 1), codes >> nc
    Wn = np.stack([((core >> k) & 1) | (((nxt >> k) & 1) << 1) for k in range(nc)], 1).astype(np.int8)
    np.save(out, Wn); log(f"{len(Wn)} windows saved to {out}"); return Wn

# ---------------------------------------------------------------------------------------------------- 2. LP
def build_lp(Wn, h, w):
    N = len(Wn); W = Wn.astype(np.int64)
    def key(cells, layer=None):
        v = np.zeros(N, np.int64)
        for k in cells: v = v*4 + (W[:, k] if layer is None else (W[:, k] >> layer) & 1)
        return v
    cell = lambda r, c: r*w + c
    left  = [cell(r, c) for r in range(h) for c in range(w-1)]; right = [cell(r, c) for r in range(h) for c in range(1, w)]
    top   = [cell(r, c) for r in range(h-1) for c in range(w)]; bot   = [cell(r, c) for r in range(1, h) for c in range(w)]
    allc  = list(range(h*w))
    pairs = [(key(left), key(right)), (key(top), key(bot)), (key(allc, 0), key(allc, 1))]
    rows, cols, vals, off = [], [], [], 0
    for a, b in pairs:
        u, inv = np.unique(np.concatenate([a, b]), return_inverse=True)
        rows += [off+inv[:N], off+inv[N:]]; cols += [np.arange(N)]*2; vals += [np.ones(N), -np.ones(N)]; off += len(u)
    rows.append(np.full(N, off)); cols.append(np.arange(N)); vals.append(np.ones(N))
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(off+1, N)).tocsc()
    b = np.zeros(off+1); b[-1] = 1
    rc, cc = h // 2, w // 2 if w % 2 else w // 2 - 1          # an interior cell with a full 3x3 around it
    c = (W[:, cell(rc, cc)] & 1).astype(float)
    return A, b, c, (rc, cc)

def culprit_mask(Wn, h, w, centre):
    """3x3 culprit test (survivor-bordered vanishing block) on the sub-window around `centre`."""
    rc, cc = centre
    sub = np.stack([Wn[:, (rc+dr)*w + (cc+dc)] for dr in (-1, 0, 1) for dc in (-1, 0, 1)], 1).astype(np.int64)
    outer = {(0,1,3,4): (5,7), (1,2,4,5): (3,7), (3,4,6,7): (1,5), (4,5,7,8): (1,3)}
    m = np.zeros(len(Wn), bool)
    for blk, out in outer.items():
        st = sub[:, list(blk)]
        van = np.isin(st, [0,1]).all(1) & ((st == 1).sum(1) >= 3) & (sub[:, 4] == 1)
        m |= van & ((sub[:, out[0]] == 3) | (sub[:, out[1]] == 3))
    return m

# ---------------------------------------------------------------------------------------------------- 3. symmetry + solve
def reflection_images(Wn, h, w):
    """Index of each window's image under the 4 reflections (identity, left-right, up-down, both).
    Reachability and every constraint family are invariant under these, so the LP is too."""
    N = len(Wn); W = Wn.astype(np.int64)
    def keys(M): return (M * (4 ** np.arange(M.shape[1], dtype=np.int64))).sum(1)
    base = keys(W); order = np.argsort(base); sb = base[order]
    idx = lambda M: order[np.searchsorted(sb, keys(M))]
    g = W.reshape(N, h, w)
    return [np.arange(N), idx(g[:, :, ::-1].reshape(N, -1)), idx(g[:, ::-1, :].reshape(N, -1)), idx(g[:, ::-1, ::-1].reshape(N, -1))]

def symmetric_objective(Wn, h, w, centre):
    """Density read off the centre cell AND its mirror images (equal for every feasible x, by consistency)."""
    rc, cc = centre
    cells = {(rc, cc), (rc, w-1-cc), (h-1-rc, cc), (h-1-rc, w-1-cc)}
    return np.mean([(Wn[:, r*w + q] & 1) for r, q in cells], axis=0).astype(float)

def ipm_solve(A, b, cost_max, threads):
    """maximise cost_max . x  s.t. A x = b, x >= 0  (HiGHS interior point).  Returns (value, row_dual)."""
    import highspy
    A = A.tocsc(); nw, nr = A.shape[1], A.shape[0]
    h = highspy.Highs(); h.setOptionValue('output_flag', True)
    h.setOptionValue('solver', 'ipm'); h.setOptionValue('run_crossover', 'off'); h.setOptionValue('threads', threads)
    lp = highspy.HighsLp(); lp.num_col_ = nw; lp.num_row_ = nr; inf = highspy.kHighsInf
    lp.col_cost_ = (-cost_max).astype(float); lp.col_lower_ = np.zeros(nw); lp.col_upper_ = np.full(nw, inf)
    lp.row_lower_ = b.astype(float); lp.row_upper_ = b.astype(float)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data.astype(float)
    h.passModel(lp); h.run()
    return -h.getInfo().objective_function_value, np.array(h.getSolution().row_dual)

def solve_lp(A, b, c, removed, threads, images=None):
    """Full LP, or (images given) the LP restricted to reflection-symmetric solutions: one column per orbit.
    A symmetric optimum always exists, so both give the same value; the symmetric LP is up to 4x smaller."""
    if images is None:
        keep = ~removed
        val, y = ipm_solve(A[:, keep], b, c[keep], threads)
        return val, y
    N = A.shape[1]
    orb = np.minimum.reduce(images); _, orb = np.unique(orb, return_inverse=True); K = orb.max() + 1
    bad = np.zeros(K, bool); np.logical_or.at(bad, orb, removed)
    M = coo_matrix((np.ones(N), (np.arange(N), orb)), shape=(N, K)).tocsc()[:, ~bad]
    log(f"symmetric LP: {M.shape[1]} orbit columns (instead of {int((~removed).sum())} windows)")
    val, y = ipm_solve(A @ M, b, M.T @ c, threads)
    return val, y

def rigorous_bound(A, c, removed, row_dual, images=None, D=10**8):
    """Exact-integer certificate.  u = rounded max-form dual.  Without symmetry: density <= max_w (c_w - (A^T u)_w).
    With symmetry the dual is averaged over the 4 reflections (the reflections act on constraint rows as signed
    permutations), which amounts to averaging (A^T u) over each window's 4 images."""
    u = np.rint(-row_dual[:-1] * D).astype(np.int64)
    r = A[:-1].T.tocsr().astype(np.int64) @ u
    keep = ~removed
    if images is None:
        z = int(((np.rint(c * D)).astype(np.int64) - r)[keep].max()); return z, D
    G = len(images); cs = np.rint(c * D * G).astype(np.int64)          # c has quarters at most; D*G keeps it integral
    rbar = sum(r[img] for img in images)
    z = int((cs - rbar)[keep].max()); return z, D * G

def export_mps(A, b, c, keep, path):
    import highspy
    A2 = A[:, keep].tocsc(); c2 = c[keep]; nw, nr = A2.shape[1], A2.shape[0]
    h = highspy.Highs(); lp = highspy.HighsLp(); lp.num_col_ = nw; lp.num_row_ = nr; inf = highspy.kHighsInf
    lp.col_cost_ = (-c2).astype(float); lp.col_lower_ = np.zeros(nw); lp.col_upper_ = np.full(nw, inf)
    lp.row_lower_ = b.astype(float); lp.row_upper_ = b.astype(float)
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A2.indptr; lp.a_matrix_.index_ = A2.indices; lp.a_matrix_.value_ = A2.data.astype(float)
    h.passModel(lp); h.writeModel(path); log("wrote", path, "(minimise -density; optimum -1/2 means density <= 1/2)")

# ---------------------------------------------------------------------------------------------------- CLI
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["enumerate", "solve", "peel", "certify", "export"])
    ap.add_argument("--h", type=int, default=3); ap.add_argument("--w", type=int, default=4)
    ap.add_argument("--threads", type=int, default=8); ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--peeled", action="store_true")
    ap.add_argument("--sym", action="store_true", help="solve the reflection-symmetric LP (up to 4x smaller, same value)")
    a = ap.parse_args(); tag = f"{a.h}x{a.w}"
    Wn = enumerate_windows(a.h, a.w, jobs=a.jobs)
    if a.cmd == "enumerate": raise SystemExit
    A, b, c, centre = build_lp(Wn, a.h, a.w); N = len(Wn)
    images = reflection_images(Wn, a.h, a.w) if a.sym else None
    if a.sym: c = symmetric_objective(Wn, a.h, a.w, centre)
    log(f"LP {tag}: {N} windows, {A.shape[0]} rows")
    excl_file = f"excluded_{tag}.npy"
    removed = np.load(excl_file) if (a.peeled and os.path.exists(excl_file)) else np.zeros(N, bool)
    dual_file = f"duals_{tag}{'_peeled' if a.peeled else ''}{'_sym' if a.sym else ''}.npy"
    def finish(val, y, rem, name):
        np.save(name, y); z, D = rigorous_bound(A, c, rem, y, images)
        log(f"LP value {val:.9f}; RIGOROUS bound: density <= {z}/{D} = {z/D:.9f}   (duals saved to {name})")
    if a.cmd == "solve":
        val, y = solve_lp(A, b, c, removed, a.threads, images); finish(val, y, removed, dual_file)
    elif a.cmd == "peel":
        cul = culprit_mask(Wn, a.h, a.w, centre); AcT = A[:-1].T.tocsr()     # transpose once, reuse every round
        for rnd in range(50):
            val, y = solve_lp(A, b, c, removed, a.threads, images)
            u = -y[:-1]; d = (AcT @ u) - c + (-y[-1])                         # reduced cost of every window (>= 0)
            if images is not None: d = sum(d[img] for img in images) / len(images)
            d[removed] = np.inf
            tight = np.nonzero(d < 1e-7)[0]; tc = tight[cul[tight]]
            if images is not None and len(tc):                            # remove whole orbits to keep symmetry
                tc = np.unique(np.concatenate([img[tc] for img in images]))
            log(f"peel round {rnd}: bound {val:.9f}; culprit windows removed this round {len(tc)}")
            if val < 0.5 + 1e-9 or len(tc) == 0: break
            removed[tc] = True; np.save(excl_file, removed)
        log(f"excluded {int(removed.sum())} window types -> {excl_file}")
        finish(val, y, removed, f"duals_{tag}_peeled{'_sym' if a.sym else ''}.npy")
    elif a.cmd == "certify":
        if os.path.exists(dual_file):
            log(f"using saved duals {dual_file} (no re-solve)"); y = np.load(dual_file); val = float('nan')
        else:
            val, y = solve_lp(A, b, c, removed, a.threads, images)
        finish(val, y, removed, dual_file)
    elif a.cmd == "export":
        export_mps(A, b, c, ~removed, f"lp_{tag}{'_peeled' if a.peeled else ''}.mps")
