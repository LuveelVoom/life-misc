#!/bin/sh
# Known-answer tests, about a minute. Needs python3 with pysat and numpy, casim, kissat.
cd "$(dirname "$0")"
set -e
fail() { echo "FAIL: $*"; exit 1; }
A=.oooo../.oo..o./o.o..oo/...oooo/o.ooo../oo..o../oo...oo/

echo "== clause cover and forced_cells vs brute force"
python3 test_core.py | tee /dev/stderr | grep -q "bad" && fail test_core

echo "== 7x7 agar, 20x20 window: 131 forced cells"
python3 -c "
import sfcore as S
rows='$A'.strip('/').split('/')
P={(r,c):int(rows[r%7][c%7]=='o') for r in range(20) for c in range(20)}
n=len(S.forced_cells(S.RuleCNF('B3/S23-a5'),P,set(P))); print(n); assert n==131
" || fail agar

echo "== GoE: 20x22 window and 389-cell patch have no predecessor, unflipped control has one"
python3 verify_goe.py B3/S23-a5 results/goe/a5_flip.txt --section 1 | grep -q ORPHAN || fail goe_full
python3 verify_goe.py B3/S23-a5 results/goe/a5_flip.txt --section 0 | grep -q ORPHAN || fail goe_patch
python3 verify_goe.py B3/S23-a5 results/goe/a5_flip.txt --section 1 --unflip 10,11 | grep -q "has predecessor" || fail control

echo "== sfcegar finds nothing quickly and doesn't crash (12x12, 50 iterations)"
python3 sfcegar.py B3/S23-a5 12 12 --inset 1 --minpop 48 --rigid 1 --maxiter 50 | tail -1 | grep -q "maxiter" || fail sfcegar

echo "all passed"
