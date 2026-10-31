# slwindows

Small still lifes containing every 3x3 window a still life can contain, for Life-like
and INT rules (`--rule`, default B3/S23). B0 and non-isotropic rules are rejected.
Not useful for anything.

```
python3 slwindows.py classify --rule B3/S23-a5   # writes windows_<rule>.json, needed by find
python3 slwindows.py find 13 13 --mode orbit --incremental --maxpop 85 --out x.rle
python3 slwindows.py find 31 31 --mode orbit --sym --incremental --out x.rle
```

`--mode oriented` needs every window, `--mode orbit` one per rotation/reflection
class. `--sym` forces D4 symmetry, so covering the classes covers every orientation.
Patterns are stable but not forced to be strict or connected. Each solution is
checked with casim plus a direct window scan, independent of the CNF.

Windows that are neither proved possible nor refuted stay legal but aren't targets;
`classify` lists them as open (raise `--radius`).

Known answers for `classify`: B3/S23 259 windows (51 classes), B3/S23-a5 336 (66),
B2-a3/S12-a 81 (16); none open. `find 8 8 --mode orbit` on B3/S23 is UNSAT.

## Results for B3/S23 (2026-10-31)

- Windows: 259 realizable (51 classes). Of the rest, 228 fail at the centre and 25
  are refuted by requiring all 9 window cells stable in a free 5x5. None open.
- Orbit: best 67 cells, 13x13 (`best_orbit_13x13_p67.rle`). 8x8 and 10x10 UNSAT
  (kissat, no proof checked). 11x11 and 12x12 undecided (11x11 got 22/50 targets in
  280 s, 12x12 42/50 in 400 s).
- Oriented: best 31x31, 420 cells and 32x32, 408 cells, both D4-symmetric. With D4:
  25x25 and smaller tested sizes UNSAT, 28-30 undecided.

Rejected: plain kissat on 12x12+ orbit and 14x14-28x28 oriented (no answer in
60-190 s); adding oriented targets one by one without symmetry (24x24 reached
56/258 and 32x32 102/258 in 200 s).
