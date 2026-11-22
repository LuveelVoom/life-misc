# selfforce

Self-forcing patches and Gardens of Eden in range-1 Moore rules (Life-like, INT, MAP; B0 rules
are rejected), with SAT (pysat/CaDiCaL). Found so far: a 20×22 Garden of Eden in B3/S23-a5.

A patch P on a finite domain D is self-forcing (Törmä–Salo) if every y with f(y)|D = P has
y|D = P, checked as UNSAT over y on D dilated by one cell. A still life containing one can't
arise from anything without it.

Needs python3 with `python-sat` and NumPy, and casim for rule tables (looked for at
`../../casim`, or set `CASIM`). `verify.py` and `verify_goe.py` also need kissat. Cells are
(row, column) with rows growing downwards; neighbourhood bits are in MAP order (NW 256 … SE 1).

```
sh test.sh                                                      # known answers, ~1.5 min
python3 goe_flip.py B3/S23-a5 <agar rows> 24 --out goe.txt      # GoE from a forcing agar
python3 verify_goe.py B3/S23-a5 goe.txt --section 1             # kissat check of a GoE
python3 agars.py B3/S23-a5 7 7 --W 20 --live 35 --dead 45       # forcing agars on a 7x7 torus
python3 sfcegar.py B3/S23-a5 16 16 --inset 1 --rigid 2 --minpop 85 --out hit.txt
```

## Files

| File | What it does |
|---|---|
| `sfcore.py` | Rule → CNF (greedy cube cover per output value), `forced_cells`, `max_self_forcing` (iterates D ← forced(D); gives the maximal self-forcing subset). Optional wraparound. |
| `sfcheck.py` | Maximal self-forcing subset of an RLE plus a dead margin. |
| `agars.py` | Enumerates agars on an a×b torus, optionally with live/dead neighbour counts restricted, and measures forcing in a tiled window. |
| `agar1.py`, `edges.py` | One agar: forced set and erosion; edge stability on strips that wrap in the other direction. |
| `sfcegar.py` | CEGAR search for a still life with a self-forcing patch: a generator SAT proposes still lifes, a checker SAT finds the predecessor with the fewest differences, and that local configuration is blocked at every equivalent position. See its docstring. |
| `sweep.sh` | The 12-configuration B3/S23-a5 sweep of `sfcegar.py` (background jobs). |
| `verify.py` | Independent check of a self-forcing still life: casim stability and kissat on a separate, unminimised encoding. |
| `goe_flip.py` | Garden of Eden from a forcing agar: flips an output cell whose 3×3 stays forced without it, crops, then unspecifies cells greedily. |
| `verify_goe.py` | Independent orphan check: kissat on an unminimised encoding. |
| `goe_sat.py` | Random-pattern orphan sampling; `has_pred`. |
| `goe.c` | Row-by-row beam search for orphans (the state is the set of predecessor row pairs). `clang -O3 -o goe goe.c`. |
| `test_core.py`, `test.sh` | Clause cover and `forced_cells` against brute force on 8 rules (including random MAP); known answers. |

## Garden of Eden, B3/S23-a5

`results/goe/a5_goe_20x22.rle`, population 237: a window of the 7×7 agar
`.oooo../.oo..o./o.o..oo/...oooo/o.ooo../oo..o../oo...oo` with row 10, column 11 (0-based)
switched on. Every predecessor of the surrounding agar has the agar's 3×3 at that cell, which
keeps it dead. `results/goe/a5_flip.txt` also has a 389-cell sub-patch (population 217) that is
already an orphan.

Both are UNSAT under `verify_goe.py` and in LLS (`results/goe/lls_*.txt`: generation 0 all `*`
on the pattern plus a one-cell ring, generation 1 the pattern with the ring and unspecified cells
`*'`). With the cell flipped back, both find a predecessor, and LLS's predecessor steps to the
window in casim. No UNSAT proof was checked.

## Self-forcing still lifes, B3/S23-a5: not found

- No cell is forced in the 97-cell still life from `../other/slwindows` (nor in its 67-cell
  B3/S23 one).
- Agars (`results/agars/`): no forcing on tori from 2×2 to 5×5 (333 agars). With live counts
  {3,5} and dead {4}, only stripe agars, none forcing. With live {3,5} and dead {4,5}, the 7×7
  agar above forces 131/400 cells of a 20×20 window and 986/1296 of a 36×36 one, but erodes to
  nothing (986 → 609 → 238 → 0). The 9×10 agar
  `.o..o.ooo./.o.oo..oo./o.oooo..oo/oo.o..oo.o/o.o.o.o.../o..o.ooo.o/.oooo.ooo./..oo.o..oo/o.o..o..oo`
  forces 150/400. The sweep was stopped after 77 of the 448 jobs in `jobs.txt`.
- Straight row and column cuts of the 7×7 agar erode completely (strips of height 24, every
  phase). 20 random wall completions (period 14, 4-row walls) also eroded completely.
- `sfcegar.py`, boxes 12–24 with insets 0–2, two of them with the 7×7 agar as core
  (`results/cegar/`). Sweep 1 (no `--rigid`): 1,200–19,000 iterations per job in about 16
  minutes, with counterexamples stuck at 1–3 differing cells. Sweep 2 (`--rigid 2`): about 2
  minutes, 200–1,800 iterations, counterexamples of at least 3 cells. Neither found anything or
  reached UNSAT; both were stopped by hand. Rejected: blocking a whole counterexample instead of
  its smallest cluster, since every clause then covered the whole box.

## Orphans by other methods, B3/S23-a5

- Random patterns (`goe_sat.py`, density 0.5): 0/200 at 6×6–12×12, 0/50 at 20×20, 0/20 at
  30×30; also 0/100 at 12×12 with density 0.3, 0.7 or 0.85. For B3/S23: 6/50 at 20×20 and 20/20
  at 30×30.
- `goe.c`, width 6, beam 3000: plateaus at 1282 predecessor row pairs by row 24
  (`results/goebeam/`). Width 7 is too slow at that beam size.

## Not tested

- The symmetric images in `sfcegar.py` against a brute-force case; `--wrap` mode never run.
- `goe_flip.py` on other agars or rules. Its crop and unspecify steps are greedy, so the GoE is
  not minimal.
- 3-flip rigidity; larger agars.
