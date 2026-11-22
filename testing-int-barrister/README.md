# testing-int-barrister

Barrister (mvr/barrister, upstream commit 9fb5e0c) ported to isotropic non-totalistic (INT) and MAP
rules, plus the catalyst searches run with it.

| Folder | What it is |
|---|---|
| `barrister/` | Unmodified upstream clone, the behaviour oracle for B3/S23. Not in the repository: clone [mvr/barrister](https://github.com/mvr/barrister) at 9fb5e0c with `--recursive`. |
| `barrister-int/` | First port: rule-aware step, still using Life's still-life solver. Superseded. |
| `barrister-int2/` | Current version: generic still-life solver driven by the rule's own 3×3 table. Any INT or MAP rule without B0 works; B0 rules are rejected. See `NOTES-INT.md` for changes, limits and fixes. |
| `work/` | Searches and post-processing (below). |

Build: `make` in `barrister-int2/` (Apple clang++; toml11 is included as a plain copy in `toml/`, so the
`.gitmodules` entry is unused).
New TOML key: `rule = "B3/S23-a5"` (Hensel or MAP). Each completed result gets an
`INT check:` line (stable in the rule, catalyst region restored).
Tests: `barrister-int2/tests-int/` (`intA.toml`, `intB.toml` find the 19-cell catalyst in
`int_example.rle` in B2i3/S23-c; `check.py <output>` re-simulates the completions with casim).
With the current source, `intA` gives 6 completions, all still lifes, in about 80 s and `intB`
gives 1 at once. `work/reg_intA.out` is from the 7 Nov build, which also reported a 26-cell
completion that isn't a still life.

Not exhaustive: unknown cells are only those inside the box marked in the input (state 2). Every
other cell is fixed OFF, and completions don't reach beyond it, so a catalyst that needs a cell
outside the box is never found.

## work/

| Folder | Search |
|---|---|
| `pi/` | B36/S23 pi catalysts (first port). |
| `q4i/`, `q4i2/`, `q4i3/` | B3/S23-q4i: catalysts for a glider + block reaction; `q4i2` is where the stale-frontier pruning bug was found. |
| `pipeline/` | Lone-object catalyst searches: B3/S23-a5 `b3o$o2bo$o2bo` (near and `far` regions) and B3/S23-q4i pi (`piq4i_*`). |
| `v2/` | B3/S23-a5 reruns comparing `barrister-int` and `barrister-int2`; `miss/` holds the cases checked for missed results. |
| `ships/` | B3/S23-a5 c/4 ships: `tAhead`, `tSide` (Tship `bo$2o$bo`), `gAhead` (glider). |
| `conduit/` | B3/S23-a5 conduits for `2bo$2bo$bobo$o2bo$b2o`. |
| `scratch/` | Ad-hoc inputs: B34k/S23 and B3/S23-q4i searches. |
| `obj4/` | B3/S23-a5 catalysts for `2bo$bobo$2ob2o$2bo`. |
| `gblock/` | B3/S23-q4i glider + block: catalysts that restore the block, kept ≥ 3 cells from everything on before generation 8 (`t8`) or 16 (`t16`) and from the glider's lane (`gen.py`, `check.py`). |
| `selfcopy/` | B2i3/S23-c: catalysts for the head `2bo$b3o$2o2bo` that stay intact ≥ 30 generations (`headcats.py`), and reactions that leave an isolated copy of the start pattern (`copies.py`). Summaries in `*_cats.txt`, `copies_now.txt`, `hw6_noneaters.txt`. |

Scripts (all need casim, found relative to this folder):
- `pipeline/genlone.py name rule rle x0 x1 y0 y1 [maxact]`: writes a search TOML whose region is a box
  minus the cells within 2 of the object.
- `pipeline/reactions.py rule base.rle out_prefix result.out ...`: keeps catalysts stable in the rule
  that recover for 10 generations, groups them by reaction (first-hit generation, leftover at gen
  300, first 8 generations after the hit), writes the smallest catalyst of each group. `--exact`
  groups by the whole reaction. Each group is labelled with its leftover (spaceships with period,
  displacement and direction, still lifes, oscillators, catalyst cells lost later);
  `--min-leftover 1` drops plain eaters. "Leftover" counts everything that isn't the catalyst at gen 300, so a
  ship escaping from a reflector shows up as leftover.
- `pipeline/lonecheck.py`: older grouping (by interacting catalyst cells and leftover).
- `conduit/conduitcheck.py rule object out ...`: results where an isolated copy of the object appears
  after the hit.
- `align.py`: finds the object's offset in a result from the Winner block. `regrid.py`: re-lays out
  summary RLEs.

Raw outputs over 1 MB are not in the repository (`.gitignore`): `gblock/t16.out` (11 MB),
`selfcopy/hw6.out` (19 MB), `selfcopy/h25.out`, `selfcopy/hcat.out`, `pipeline/piq4i_far6.out`
(7.9 MB), `pipeline/piq4i_far6.snap`, `pipeline/piq4i_N.out`, `pipeline/a5_farW.out`, `pi/R6.out`,
`pi/R6.fin`. Their summaries are.
