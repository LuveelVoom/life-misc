# Handoff: density of oscillating agars in Conway's Life

Start here. This file is the complete state of a long research conversation, written so a new
conversation (on any account) can continue it. Upload this whole zip, or at least this file plus the
PDF and `scripts/`, then paste `CONTINUE_PROMPT.md` as your first message.

Packaged 2026-10-17.

## 1. The problem

An **oscillating agar** is a pattern that fills the plane periodically (equivalently, lives on a finite
torus) and returns to itself after p generations under Life (B3/S23). **Conjecture:** its average density
over a period is at most 1/2. Still lifes: proved by Elkies (1998, arXiv:math/9905194). Oscillators: open.
Best published bound found: 1176/2087 ≈ 0.563 (Keith Amling, 2024, per LifeWiki "Density" page).

## 2. Results obtained (all still need independent verification)

Computer-assisted (single code base, `scripts/`):

| | Result | Checked by |
|---|---|---|
| A | Every oscillating agar, any period: density ≤ 205847955/400000000 ≈ 0.5146199 (3×4 windows, reflection-symmetric LP) | rounded dual, exact integer arithmetic |
| B | Same with 3×3 windows: ≤ 51612905/100000000 ≈ 0.5161291 | rounded dual, exact integer arithmetic |
| C | **Exclusion theorem:** if none of 769 listed 3×3 two-generation window types occurs, density ≤ 1/2 | exact rational dual (HiGHS basis + FLINT), z(u) = 1/2 exactly |
| D | If none of 3,350 listed 3×4 window types occurs, density ≤ 0.500000010 | rounded dual, exact integer arithmetic. **Reproduced on the user's own machine** (log in `your_files/`) |
| E | Density ≤ 1/2 for agars with period h in one direction: p=1: h≤6; p=2: h≤5; p=3: h≤3; p=4,5: h≤2; p=6: h=1 | Bellman–Ford integer certificates |

By hand:

| | Result | Note |
|---|---|---|
| F | Any periodic set where every member has ≤ 3 neighbours in the set has density ≤ 1/2 | **This is Elkies' theorem.** Our discharging proof is ours; rules resemble Chu–Stuckey–de la Banda (2009) area rules |
| G | Period-2 agar with no two adjacent X cells and no two adjacent Y cells: density ≤ 1/2 | correct but narrow |
| H | Density α ≤ (k+4)/12, k = average number of in-set neighbours of set cells | standard eigenvalue (Hoffman-type) bound; tight for stripes and full plane; not useful for the main problem |
| I | Band-forcing lemma: the family-1 band forces `YY.XX.` above and `XX.YY.` below; band + forced rows has density 20/42 < 1/2 | period 2; SAT-checked over 30 columns without assuming periodicity |
| J | The `^.^.` edge rows of family 2 cannot occur (strips repeating every 2 columns) | what remains is the plain double stripe, density exactly 1/2 |

**The LP.** Variables: frequency x_w of each reachable two-generation window (87,851 for 3×3;
2,175,813 for 3×4). Constraints: left/right and top/bottom overlap consistency, stationarity (time-t layer
frequencies = time-(t+1) layer frequencies; this is what makes it period-independent), normalisation.
Objective: centre alive at t. Weak duality: any multipliers u give density ≤ max_w (c_w − (Aᵀu)_w).
LP optima (floating point, NOT proved exact): ≈0.5161290 (matches 16/31) for 3×3; ≈0.5146198 (matches 88/171) for 3×4.

**Exclusion ("peeling").** Remove culprit windows that are tight (zero reduced cost) at the current
optimum, re-solve, repeat. Culprit window = centred on a cell of a *vanishing block* (2×2 block, ≥3 live at t,
none at t+1) with a surviving edge-neighbour outside the block. Removing all 2,480 culprit types gives exactly 1/2.
Peeling: round 0 ≈0.516129 → remove 405 (family 1); round 1 0.505279 → remove 364 (family 2); round 2 = 1/2.
**This is a heuristic**: different runs gave 458, 464 or 769 windows; each run ending in an exact 1/2 proves
the theorem for its own list. No systematic (e.g. minimal) exclusion search has been done.

**The two families.**
- Family 1, blinking-block band (5 rows, period 6 columns, density 8/15):
  ```
  XXSYYS
  XX.YY.
  .SS.SS
  YY.XX.
  YYSXXS
  ```
  (S always alive, `.` always dead, X alive phase A only, Y alive phase B only.) Real agar containing it:
  6×14 torus, density 10/21: `x = 6, y = 14, rule = B3/S23:T6,14` then
  `3o2bo$2o$b2ob2o$3b2o$2b4o$2o$2o$2b4o$3b2o$2ob2o$2o$3o2bo$3b2o$3b2o!`
- Family 2: moving double stripe with half-filled `^.^.` edge rows (edge rows impossible, result J).

**Structure of the LP optimum (SAT, CaDiCaL via PySAT).** The 98 windows in the support of one 3×3 optimum:
8×8 patch exists, 9×9 and 10×10 do not (re-run 2026-10: genuinely UNSAT, not budget timeouts). Strips 3–5 rows
fit every width tested (4–30); 6 rows: widths 4–10 yes, 14/20/30 no; 7–9 rows: widths 4–8 yes, 10/14/20/30 no.
Using all 6,851 tight windows: a culprit can centre a 9×9 patch but not 11×11 (re-run: UNSAT); the resulting
cut only moves the LP value 0.516129 → 0.516119. No DRAT proofs were produced.

**Other band checks.** Real 6-wide period-2 agars containing the band: none at heights 7, 8, 10, 12, 15–18;
densities 0.476, 0.456, 0.467 at heights 14, 19, 20; heights 9, 11, 13 untested. Band piece m columns long
(+1 row each side), max excess of live cell-phases over half: +8 (m=6), −1 (m=12), −3 (m=18); only these m checked.

## 3. Approaches that failed (don't repeat without a new idea)

Fixed-window LPs/charge rules (blocked at ≈0.5161 / ≈0.5146); 4×4 windows only around dying triples (≥0.5149);
requiring window history (all 98 optimal windows arise after 4 generations; deeper not checked);
slicing / one phase per cell (refuted by a triangle agar, 8×10 torus, phase A
`x = 10, y = 8, rule = B3/S23:T10,8` `b2o$5bob3o$4ob2obo$3b2o$2b2ob5o$obobo2bo$obobo$o2bo5bo!`);
pushing charge one generation forward (fails on ~12% of 5×5 patterns); radius-1 charge sharing / coverage lemma
(CP-SAT found real counterexamples — the `cex_R1_H1_*.txt` files in `your_files/`); local surgery; S-to-D
matching along lines (7/17 random period-2 agars fail); local counting around one blinking block (isolated
block's ring can have S − D = 4; ring identity c + 2e = 12); splicing identical rows of a smallest
counterexample (no bound); reduction to bounded period (no proof); spectral/pairwise positivity.

## 4. The critique of the write-up and what was done

An outside reader sent this critique of the first PDF (verbatim summary of their points): R1 misworded (it is a
statement about a single configuration) and not new; R2, R3 rewordings of R1; R4 a one-line argument, not a lemma;
R5 "AI tryharding", probably incorrect; R6 "total bs" though such a bound would be nice; R7 useless (no machinery
uses it); R8 appeared skipped; R9 proves nothing; the exclusion idea is "cool" but the search was unsystematic;
the 9×9 impossibility is the most interesting finding but its method wasn't explained; appendix useless.

Assessment (after two full double-checks): critique fair on R1–R4, R7, R9, R8 numbering, unsystematic peeling,
missing 9×9 method, appendix. **Wrong on R5 and R6**: both are correct (re-derived; R6 also tested on 20,000
random tori). Double-checks additionally found and fixed: R1 is exactly Elkies' theorem; untested heights
(9, 11, 13) had been claimed; width/length claims overstated; 16/31 and 88/171 presented as exact when they are
numerical matches; "8 digits" should be ~7; wrong cross-reference; `.SS.SS` row misdescribed; Section 6 wording
didn't match the script's `--width` meaning; history-depth claim overstated. Confirmed: degree-lemma rules (exact
check over all relevant windows), reflection-symmetry property (3×3, 2×3, 3×2, 2×4; **not checked on 3×4 itself**),
exact 1/2 reproduction, 8×8 patch validity.

`writeup/oscillating_agar_density_writeup.pdf` is the **final revised version** (20 pages) with all of this applied.

## 5. Open problems / next steps

1. The conjecture is equivalent to: the 769 excluded window types never push density above 1/2 in any oscillator.
2. Arrangements: what prevents short/irregular band pieces from being packed densely?
3. Systematic exclusion search (smallest / most structured exclusion set; larger windows).
4. Verification: independent reimplementation of enumeration, constraints and certificate checks; DRAT proofs
   for the SAT results; check the reflection property on 3×4 directly.
5. Bounded-period reduction: if a counterexample exists, is there one with bounded spatial or temporal period?

## 6. What is in this zip

| Path | Contents |
|---|---|
| `writeup/` | final PDF |
| `scripts/` | all scripts + README (canonical). `coverage_search.py` is the older v1 of `coverage_search2.py` |
| `patterns/` | RLE files for the family-1 band, family-2 strip, and the 6×14 real agar |
| `data/reach.npy` | boolean mask over 2^18 window codes, 87,851 reachable 3×3 windows (regenerated by scripts anyway) |
| `your_files/` | files the user uploaded: coverage-lemma counterexamples (`cex_R1_H1_type*.txt`, output of `coverage_search2.py --radius 1 --history 1`); the 3×4 peel log from the user's Mac; an ortools/NumPy 1.x-vs-2.x import error log (fix: use the fresh `agar` conda env from the README, not base anaconda) |

**Not included (lost when the workspace was reset; all regenerable or non-essential):** the full conversation
transcript; the PDF's source script; `certificate.pkl` / `excluded_windows.pkl` (regenerate:
`python life_agar_density.py certify`, ~5 min); the saved LP optimum used for the patch/SAT checks and those
small exploratory scripts (regenerate by solving the 3×3 all-period LP, taking the windows with x_w > 1e-9, and
running a SAT tiling encoding: one selector per window per interior position, each selector fixes the 9 cells'
two-time states, at-least-one selector per position). Two Claude Docs write-ups on the old account are earlier,
superseded versions of the PDF.

## 7. Reproduction quick start

```
conda create -n agar python=3.12 -y && conda activate agar
pip install numpy scipy highspy python-sat python-flint ortools
cd scripts
python life_agar_density.py certify                         # result C, ~5 min, ends "EXACT BOUND = 1/2"
python big_window_lp.py certify --h 3 --w 3                 # result B, ~1 min
python big_window_lp.py solve --h 3 --w 4 --sym --threads 16   # result A, ~10 min on 16 threads
python big_window_lp.py peel  --h 3 --w 4 --sym --threads 16   # result D, ~45 min
python strip_certificate.py --width 5 --period 2            # result E (one case), ends PROVEN
python band_forcing_check.py                                # result I solver check
python check_windows.py                                     # where excluded windows occur (after certify)
```
