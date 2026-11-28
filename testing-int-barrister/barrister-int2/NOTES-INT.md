# Barrister INT port, generic still-life solver (barrister-int2)

New TOML keys:
- `rule = "B2i3/S23-c"` (Hensel INT or MAP; default B3/S23). B0 rejected.
- `unsafe-allow-unrepresented-survival = true` (default false): accept rules with survival
  outside counts 2-3 (e.g. S4i). The still-life solver cannot represent those cells, so
  still lifes using them are never found.

Changes vs upstream (mvr/barrister 9fb5e0c):
- `IntRule.hpp` (+ `hensel_*.inc` from ntCatForce-plus): rule tables, exact step,
  three-valued step table over all 3^9 neighbourhoods, neighbour-count shadows.
- `LifeUnknownState.hpp`: all three `StepMaintaining*` replaced by a per-cell version
  (exact where the 3x3 is known; stable value where the 3x3 matches the background;
  three-valued lookup otherwise). The Life circuits in `bitslicing/` are unused.
- `Barrister.cpp`: `OptionsFor` uses the rule's count shadows; branching is on the
  value of an unknown stable cell next to the frontier cell instead of on the
  frontier cell's transition; interaction start is set in `TryAdvance`; each
  completed result gets an exact `INT check:` line.
- `LifeAPI.h`: `Step`/`StepFor` use the rule table for non-Life rules.

Limits:
- Still-life solver and completion are still Life's. Accepted rules: survival only
  at counts 2-3, B3 complete. Completions may be Life-stable but not stable in the
  rule (flagged `still-life=NO`); results whose only INT-stable completion was
  missed are lost. Not exhaustive.
- The `catalyst-region-restored` flag looks at one generation only; repeated or later
  hits make it say NO for valid catalysts.
- Slower than upstream: no count reasoning in the step, and the frontier is
  recomputed in full after every branch.

Tests (`tests-int/`): `intB.toml` (top of the catalyst given) and `intA.toml` (whole
catalyst unknown) both find the catalyst from `int_example.rle`; `check.py <output>`
re-simulates completed results with casim.

Fixes:
- 2026-11-07: branches were pruned when the chosen frontier cell had no unknown stable
  neighbour. The cause was a stale `frontier.next` (the neighbourhood had become known after
  it was computed). Now the frontier is recomputed instead. Found because a valid boat
  catalyst (B3/S23-q4i test in work/q4i2) was missed in a large region but found in a small
  one. After the fix: minimal.toml 16 results (was 7), searches run longer.
  `INTDBG=1` in the environment logs the remaining prune points to stderr.
- Remaining crude prune: a frontier cell whose 3x3 holds active (non-stable) unknowns
  and no stable unknown is still pruned (`PRUNE active-unknown-neighbour`).

## barrister-int2: generic still-life solver (2026-11-07)

`barrister-int2` replaces upstream's Life-specific still-life solver (per-cell neighbour-count
options live2/live3/dead0..dead6 and the espresso circuits in `bitslicing/`) with one that uses
the rule itself, so any INT or MAP rule without B0 works and `unsafe-allow-unrepresented-survival`
is ignored. The old file is kept as `LifeStableState.life.hpp.orig`.

- `IntRule.hpp`: table `still[(unk9<<9)|on9]` over all on/off/unknown 3x3 windows: whether some
  completion is stable at the centre, and which unknown cells are forced.
- `LifeStableState.hpp` (rewritten, same interface): `Propagate` applies the table to every window
  around a known ON cell, incrementally (only windows near cells changed since the last check;
  `fullCheck` after direct assignments). Windows without a known ON cell can always be completed
  with OFF cells (no B0), so no inconsistency is missed; forced cells there (B1/B2 rules) are not
  deduced, which only weakens pruning. `TestUnknowns` tries both values of unknown cells next to
  ON cells. `CompleteStable` branches on unknowns next to ON cells, then checks "rest OFF" exactly.
  `dead0` now means "ON or next to a known ON cell" (the catalyst's zone, which is how Barrister
  uses it).
- `LifeUnknownState.hpp`: since the background is a real still life of the rule, every cell whose
  3x3 equals the background keeps its background value without a lookup.
- `Barrister.cpp`: `SetForced` is generic: for a frontier cell, each allowed transition is checked
  by enumerating the unknown stable cells of its 3x3; impossible transitions are removed and cells
  equal in every remaining case are forced. The count-based `OptionsFor` is gone.

Measured (single thread, M4 Pro; other searches were running, so timings are rough):
- B3/S23 `minimal.toml`: 17 results, all stable; search nodes 97k vs 83k for barrister-int
  (Life counts), 3 s vs 1 s (completion with minimise-results is the difference).
- B2i3/S23-c tests A/B: still find the 19-cell catalyst; all completions stable (the old build
  produced unstable ones).
- B3/S23-a5: no catalyst unstable in the rule (old build: 48 of 267). The search tree is larger
  (4.7M vs 0.75M nodes on one test) because still lifes with S5 cells are now allowed; the old
  solver excluded them.

Rejected idea: keeping each frontier transition as a lasting window constraint on the stable
cells (like upstream's per-cell options): no fewer nodes on the a5 test (30 s vs 26 s), removed.

## Results file and progress (2026-11-27)
- Results go to `output-file`, default `<input>-results.rle` (next to the .toml). Normal mode rewrites
  it after every new solution (via `<file>.tmp` + rename) and once more at the end; the `#C` first line
  says whether the search was still running. The file holds all solutions (the "All solutions" list);
  the trimmed and filtered lists are still only printed. Metasearch writes it after each round, as before.
- stderr gets a line every 10 s: `progress: P% (rough), path ..., nodes, solutions, seconds`. P counts
  each node as 2^-depth of the tree and both halves of a branch as equal, so it reaches exactly 100%
  at the end but can sit near 0 for most of a run (intA: 0.0000% until 90 s of 92 s). `path` is the
  current node's first 32 branch choices (0 = cell set off first, 1 = the second branch).
- Search output on stdout unchanged (intA, intB byte-identical to the previous build).
  The previous binary is `Barrister.pre-progress`.
