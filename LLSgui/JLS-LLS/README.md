# JLS-LLS

JavaLifeSearch 1.7 with one addition: **Search → Run in LLS (INT rule)...** (Ctrl+L). It takes the
field you are editing, runs it in Logic Life Search with the rule you type, and shows the output.

Start it with `Run JLS.command`, which builds on first use. Requires a JDK (9+), a C++ compiler,
`python3`, and kissat or cadical (`lls/solvers/` links to the Homebrew copies).

The JLS 1.7 source has no license text, so it isn't in this repository. The build takes it from
the JLS download, `JLS-1.7.zip` from the ConwayLife.com JLS thread: put it next to this folder
or set `JLS_ZIP`.

## The dialog

- **INT rule**: any rule LLS parses, including partial rules (`pB3/S23`). Blank uses the JLS rule.
- **Solutions to find**: LLS `-n`.
- **Generation 0 must not be empty**: `--population >0`. JLS has no such option.
- **Exclude lower periods**: `--force_change 0 P/q` for each prime q dividing P. Skipped when there
  is a shift or translation. JLS would accept lower periods, so untick it to match JLS.
- **Extra LLS arguments**: split on spaces, no quoting.
- **SAT solver**: kissat, cadical, or `llsbt (backtracker)`. The backtracker is a lifesrc-style search
  in C++ with its own threads. It's erratic: often much faster than SAT, sometimes far slower. It needs
  a fully specified rule, and refuses B0 rules and active cells per layer limits.
- **Parallel SAT copies, or backtracker threads**: with SAT, this many copies of LLS run with different
  seeds, and every other copy prefers dead cells. The first to answer wins and the rest are killed.
  Any copy's "Unsatisfiable" is a full proof. With `--parameters` in the extra arguments, seeds are lost.
- **`-o FILE` in the extra arguments**: only the winning copy's solutions are appended.
- **Stop** kills LLS and the solver. So does closing JLS.

The first line of the log is the LLS command. Solutions are printed as RLE, with a one-cell border
around the field. Reset a JLS search before running: only the edited field is sent.

## What's where

| Path | What it is |
|---|---|
| `src/jls/LlsRunDialog.java` | The new dialog. |
| `ContentPane.patch` | The menu item and 3 small accessors, applied to the JLS original when building. |
| `lls/` | LLS from `logic-life-search-master(1).zip`, the newer of the two zips, with the patches below. |
| `backtrack.py`, `native/llsbt.cpp` | The backtracker, from `LLS-GUI`. |
| `builder.py` | `LLS-GUI`'s LLS-as-a-library layer, used by `backtrack.py`. It sets the default background to vacuum. |
| `test/run.sh` | Known-answer test, a few seconds. |

## Changes to LLS's JLS import

Each is marked `JLS-LLS` in `lls/lls` or `lls/src/formatting.py`. Without them, LLS runs a
different search with no warning. The JLS source (`CellArray.getCellPtrEx`) is the reference.

- **Shift and translation** were ignored, so a spaceship search became an oscillator search. Now
  generation P is free and tied to generation 0 by `--symmetry pP xR yD <transform>`.
- **Negative shifts** (`x-1`) are now parsed.
- **4-Fold-Diagonal** was missing and gave no symmetry at all. Now `D4X`.
- **`--symmetry` on the command line** replaced the file's symmetry. Now they combine.
- **`layers_active_cells_variables_only=No`** never took effect (a `== None` test on a `False` default).
- **Cell values** the importer doesn't know now stop with an error. They used to become free cells
  in generation 0. JLS's error flag (32) is masked off.

- **`src/primer.py`** used `type` statements (Python 3.12+). Started from Finder, `python3` is the
  system 3.9, so they are now plain assignments. Run under 3.9 and 3.12.

Refused with an error, as LLS can't express them: horizontal or vertical tiling, outer space unset.

Warned about and ignored: the "variables only" option of the generation 0 cell limit, and the live
cells per layer limit.

## Conventions

- JLS: `cell(P, c, r) = cell(0, T(c - right, r - down))`; off the field is background.
- LLS: `cell(0, x, y) = cell(P, f(x, y))`, so f is T's inverse followed by the shift. That is why
  JLS Rotate-90 maps to `RO3` and Rotate-270 to `RO1`.
- x is the column, y is the row, and y increases downwards in both programs.

## Not tested

- Against JLS's own search: no comparison of solution counts on the same field.
- Rotate translations on a real search. The mapping is checked cell by cell, for all 8 translations
  (`test/wrap_agreement.py`).
- B0 rules. Solver runs with subperiods, frozen, unchecked or unset cells.
- Parallel SAT with cadical's `--phase=0`, and with several solutions (`-n` above 1).
- Whether the backtracker and SAT agree on solution counts.
- Whether JLS can save cell values outside the importer's list, such as frozen ON.
- `layers_active` constraints through this route.
