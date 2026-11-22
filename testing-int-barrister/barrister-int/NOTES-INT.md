# Barrister INT port (crude prototype)

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
