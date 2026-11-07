#!/usr/bin/env python3
"""All-period LP bound for oscillating agars in Conway's Life (3x3 two-generation windows).
Requirements: pip install numpy scipy highspy        Usage: python bound_16_31.py
Cell states: 0 dead->dead, 1 alive->dead, 2 dead->alive, 3 alive->alive; windows are row-major 3x3."""
import numpy as np, highspy
from scipy.sparse import coo_matrix

# 1. Reachable windows: run every 5x5 pattern at time t, record (3x3 core at t, 3x3 core at t+1).
reach = np.zeros(1 << 18, bool)
for s in range(0, 1 << 25, 1 << 22):
    P = np.arange(s, s + (1 << 22), dtype=np.int64)
    core = np.zeros_like(P); nxt = np.zeros_like(P)
    for k in range(9):
        r, c = divmod(k, 3); r += 1; c += 1
        cc = (P >> (r*5 + c)) & 1
        n = sum((P >> ((r+a)*5 + (c+b))) & 1 for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0))
        core |= cc << k; nxt |= ((n == 3) | ((cc == 1) & (n == 2))).astype(np.int64) << k
    reach[core | (nxt << 9)] = True
codes = np.nonzero(reach)[0]
W = np.stack([((codes >> k) & 1) | (((codes >> (9+k)) & 1) << 1) for k in range(9)], 1)
N = len(W); print("reachable windows:", N)

# 2. Constraints: left/right 3x2 overlaps, top/bottom 2x3 overlaps, stationarity (t layer = t+1 layer), normalisation.
def key(cells, layer=None):
    v = np.zeros(N, np.int64)
    for k in cells: v = v*4 + (W[:, k] if layer is None else (W[:, k] >> layer) & 1)
    return v
pairs = [(key([0,1,3,4,6,7]), key([1,2,4,5,7,8])), (key([0,1,2,3,4,5]), key([3,4,5,6,7,8])),
         (key(range(9), 0), key(range(9), 1))]
rows, cols, vals, off = [], [], [], 0
for a, b in pairs:
    u, inv = np.unique(np.concatenate([a, b]), return_inverse=True)
    rows += [off+inv[:N], off+inv[N:]]; cols += [np.arange(N)]*2; vals += [np.ones(N), -np.ones(N)]; off += len(u)
rows.append(np.full(N, off)); cols.append(np.arange(N)); vals.append(np.ones(N))
A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(off+1, N)).tocsc()
b = np.zeros(off+1); b[-1] = 1
c = (W[:, 4] & 1).astype(float)                      # centre alive at time t
print("rows:", A.shape[0])

# 3. Solve: minimise -density.
h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('solver', 'ipm')
lp = highspy.HighsLp(); lp.num_col_ = N; lp.num_row_ = A.shape[0]; inf = highspy.kHighsInf
lp.col_cost_ = -c; lp.col_lower_ = np.zeros(N); lp.col_upper_ = np.full(N, inf)
lp.row_lower_ = b; lp.row_upper_ = b
lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
lp.a_matrix_.start_ = A.indptr; lp.a_matrix_.index_ = A.indices; lp.a_matrix_.value_ = A.data.astype(float)
h.passModel(lp); h.run()
print(f"LP optimum: {-h.getInfo().objective_function_value:.9f}   (16/31 = {16/31:.9f})")

# 4. Certificate: round the dual to integers over D and take max_w (c_w D - (A^T u)_w) in exact integer arithmetic.
D = 10**8
u = np.rint(-np.array(h.getSolution().row_dual)[:-1] * D).astype(np.int64)
z = int(((c * D).astype(np.int64) - A[:-1].T.tocsr().astype(np.int64) @ u).max())
print(f"RIGOROUS bound: density <= {z}/{D} = {z/D:.9f}")
