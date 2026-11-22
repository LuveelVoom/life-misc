"""Forced cells and maximal self-forcing subset of a W x W window of one agar.
usage: agar1.py RULE AGAR_ROWS W   (rows separated by '/', 'o' live)"""
import sys, time, sfcore as S, sfcheck
rule, rows, W = sys.argv[1], sys.argv[2].strip('/').split('/'), int(sys.argv[3])
a, b = len(rows), len(rows[0])
P = {(r, c): int(rows[r % a][c % b] == 'o') for r in range(W) for c in range(W)}
rc = S.RuleCNF(rule); t = time.time()
F = S.forced_cells(rc, P, set(P))
print("first-pass forced", len(F)); print(sfcheck.show(P, F))
D = S.max_self_forcing(rc, P, F)
print("max self-forcing", len(D), f"{time.time()-t:.1f}s")
if D: print(sfcheck.show(P, D))
