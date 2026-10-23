# Findings (2026-10-22)

## A rule where random soups emit ships in all three direction classes
Found by `explore` (examples/diag_ship_explore.cfg, `reward oblique 20`), then improved with `climb`
(findings/three_ship_climb.cfg, 400 iterations).

    B2ce3akqy4enjytw5ckq6ek8/S12en3ejq4ainjwz5jqry6ain7e8        (climbed; rates from 20000 fresh soups)

| object | cells | per 16x16 soup |
|---|---|---|
| c/4 diagonal ship `2bo$obo!` | 3 | ~3 (capped at 4 per soup) |
| c/2 orthogonal ship `bo$obo2$obo!` | 5 | ~3 (capped) — uncapped ~5.7 in the parent rule |
| **(2,1)c/9 knightship** `obo$2o$2bo!` | 5 | **0.80** |
| 3c/10 orthogonal ship `2o2$o2$2o!` | 5 | ~0.09 |

bad soups 3.2%, mean cost 1.16M cell-updates. In 43 of 20000 soups the soup also contains a c/2-ship emitter
whose population keeps growing with ship deletion off (gun candidate; soups are logged as RLE when found).
Parent rule (straight from explore, knightship 0.64/soup, bad 1.9%):

    B2ce3akqy4enjytw5cknq6ak8/S12en3ejq4ajwz5-eikn6-an7e8

## Common small knightships (each in its own rule; all checked with the original engine too)
| ship | cells | rate/soup | rule |
|---|---|---|---|
| (3,2)c/19 `2o$bo$o2$o!` | 5 | 0.10 | B2ce3akqy4ejytw5cknq6ak8/S12in3ajq4ejywz5cajqr6-an7e8 (the start rule itself) |
| (3,1)c/12 `bo$o$o$o!` | 4 | 0.34 | B2ce3aikqy4ejytw5ceknq6ak8/S13aqy4eywz5cajqr6-an7e8 |
| (5,4)c/25 `2bo$2o$2bo$o!` | 5 | 0.11 | B2ce3akqy4einjytw5caknq6ak8/S12in3eajq4eqywz5cajqr6-an7e8 |
| (4,1)c/16 `2bobo$3bo$3bo$bo$o!` | 6 | 0.005 | B2ce3eakq4eknytw5-eijr6cak8/S12i3kjqy4ejywz5canjq6-n7e8 |

The 12/10/17-cell (3,2)c/19 ships in discoveries_obl.txt are much rarer (~2e-4 per soup); climbing on them got
no signal even with 20000 soups per evaluation.

## Variety runs (2000 explore iterations each, same start)
C1 soups: 179 new behaviours + 74 new soup behaviours (puffers/rakes/guns). D4 soups: 160 + 88.
Only 99 behaviours were found by both; 151 only with C1 and 149 only with D4.
Logs: findings/explore_variety_C1_discoveries.txt, findings/explore_variety_D4_discoveries.txt.

## Probe (catalysts / signals / fossils), 2026-10-22
- Life sanity check (`probe B3/S23`, 20000 soups, `probe_stills=40`): eater 1 and its 9-cell variant eat gliders
  and LWSSes; the diagonal line of cells is a speed-c fuse.
- Three-ship rule above: the only common still life is the 2-cell domino, and it eats the (2,1)c/9 knightship.
- 300-iteration explore with `probe 1` from the diagonal-ship start rule (findings/explore_probe_discoveries.txt):
  catalytic setups found — 132 eater, 81 converter (e.g. a 9-cell c/2 ship turned into a 4- or 6-cell c/2 ship by
  an 8-cell still life), 23 splitter, 25 pass-through (a c/4 ship crossing straight-line still lifes of 13–16
  cells and coming out unchanged), 117 other catalysed reactions; fuses at 12 different speeds on orthogonal and
  diagonal wicks; no reflector and no non-destructive signal yet.
