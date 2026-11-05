# intcover

Hunts for the smallest pattern containing all 102 INT transitions (51 B, 51 S) at
least once, on a dead infinite plane. B0 is covered by the background.

```
clang -O3 -mcpu=native -o sa sa.c -lm
./sa 15 15 120 1 > run.txt                 # annealer: W H seconds seed
python3 intcover.py 12 12 --feasible       # CP-SAT: is this box possible?
python3 intcover.py 15 15 --hint run.txt --time 300 --out best.rle
```

`intcover.py` needs ortools and stubs out pandas, whose anaconda install is broken.
It checks each solution itself and with casim: for each transition T, the rule
containing only T must leave a live cell after one step.

## Results (2026-10-31, 10 threads/processes)

| Box | Result | Method |
|---|---|---|
| 10x10 | impossible | CP-SAT, 15 s; solver verdict, no certificate |
| 11x11, 11x12, 10x12 | undecided | CP-SAT 90 s; annealer stuck at 5 and 2 missing on the first two |
| 12x12 | population 73 | annealer |
| 13x13, 14x14 | population 70 | annealer |
| 15x15 (bounding box) | **population 68**, in `results/best.rle` | annealer, run in a 16x16 box |

Bounds on population: at most 68 (found), at least 55.

- 51 by counting: the 51 live-centre transitions need 51 different live cells
  (equally, 102 transitions total 459 cells and each live cell is in 9 windows).
  CP-SAT never raised this bound in the box runs above.
- 55 by an LP on window counts n_w, valid in any box: minimise population, with
  n_w >= 1 summed over each transition class, and each nonzero 3x2 block appearing
  equally often as the left and right part of a window (2x3 likewise, top and
  bottom). Optimum 54.5, for 3x3 windows and for 4x4 tiles; the LP duals were checked
  in exact fractions, with no external proof checker. `results/best.rle` satisfies every
  constraint exactly, objective 68.
- Evidence only: the integer version reached a bound of 56.42 (so 57) in 10 s of
  floating-point SCIP, stopped early. Not a proof.

Rejected: warm-starting CP-SAT from the 69-cell pattern in 15x15; after 240 s it
was only at 78, because the hint fixes cell variables only.
