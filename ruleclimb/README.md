# rulesearch

Constrained search over range-1 Moore rules (Life-like, isotropic non-totalistic, or MAP).
Constraints are *exact*: patterns must evolve exactly as in a reference rule. Within those
constraints it either climbs toward making target objects common (`climb`) or looks for
new objects while keeping targets common (`explore`).

    make                                      # OpenMP build (macOS: picks Homebrew g++-N automatically)
    make nomp CXX=clang++                     # single-threaded build with any C++17 compiler
    ./rulesearch CONFIG check [RULE...]       # forced transitions, pattern periods, test rules
    ./rulesearch CONFIG census RULE [N] [SEED]
    ./rulesearch CONFIG climb ITERS
    ./rulesearch CONFIG explore ITERS [POPCAP]
    ./rulesearch CONFIG dump RULE SEED        # soup as RLE
    ./rulesearch CONFIG probe RULE [N] [SEED] # catalysts, signals/fuses on wicks, fossils (see below)
    ./rulesearch CONFIG climb 500 seed=3 soups=400 soup=16,16,0.4   # key=value overrides config lines
    OMP_NUM_THREADS=n sets threads; results do not depend on thread count.

Ctrl-C during `climb`/`explore` finishes the current iteration and prints the best rule / final population
as usual (a second Ctrl-C quits at once). Example configs are in `examples/`.

## Conventions (followed Golly source `gollybase/liferules.cpp`)
- 512-entry table index: row-major from NW, NW = 256 ... SE = 1, centre = 16 (MAP order).
- MAP: 86 base64 chars, 6 bits per char MSB first, bit k = table[k].
- Hensel letters: Golly's neighbourhood tables (current letter set, not the older v/y/r variant).
- RLE: x right, y down. Output RLE wraps only between runs.

## Config (one directive per line, `#` comments)
| key | meaning |
|---|---|
| `rule R` | reference rule (Hensel `B../S..` or `MAP...`). Constraint evolutions come from it. |
| `symmetry totalistic\|isotropic\|none` | search variables = symmetry classes (18 / 102 / 512). Reference rule must respect it. |
| `pattern NAME RLE [split]` | constraint: evolve until the state repeats (translation allowed) or dies; every transition used is frozen. RLE may also be `@file.rle` (an RLE file; header and `#` lines are skipped). `split` makes one constraint per component (components = groups 3+ cells apart), named NAME.0, NAME.1 … in order of top-left corner; it aborts if components interact. |
| `target NAME [W]` | objective term `W*log(rate)`, rate = mean per-soup count (capped). Target object = the cycle the pattern reaches, in every orientation allowed by the symmetry. |
| `min_rate NAME r` | feasibility threshold for a target (per soup) |
| `fix CLASS v` | extra forced transition: `B2a 0`, `S3 1` (totalistic), `i273 1` (MAP index) |
| `soup W H D` | soup size and density (default 16 16 0.5) |
| `soups N` | soups per evaluation (default 1000) |
| `grid G` `gmax` `pmax` | universe G×G (256), generation cap (2000), population cap = explosion (2500) |
| `cap c` | per-soup cap on target counts (4) |
| `max_bad f` | max fraction of exploded / unstabilised / too-slow soups |
| `max_cost M` | max mean cost (million cell-updates per soup) — keeps soups short |
| `cost_weight w` `cost_ref M` | soft penalty `w*log(cost/cost_ref)` above cost_ref (default: start rule's cost) |
| `soup_cost_limit M` | abort a single soup past M million cell-updates (counted bad). Default min(20, 10*max_cost) |
| `min_osc_period p` | oscillators below this period are not "interesting" (default 3) |
| `maxp P` | max period searched when classifying (256) |
| `start RULE` | starting rule(s); forced classes are overwritten to satisfy constraints |
| `lahc L` `reseed R` `seed S` | climb history length, new soup set every R iterations, RNG seed |
| `log FILE` `archive FILE` | discoveries log (RLE) and persistent object archive (one id per line); `none` (or `log=` on the command line) turns either off |
| `best FILE` | climb: best rule so far, rewritten on every improvement; explore: population, rewritten every 25 iterations |
| `reward KIND W` | explore: multiply novelty credit of discoveries of KIND (`ship`, `oblique`, `diagonal`, `orthogonal`, `osc`, `growth`, `gun`, `catalyst`, `signal`) by W; several allowed, they multiply |
| `soup_symmetry S` | `C1` (default, original soups), `D2` (mirror, 2W×H), `C2` (180°, W×2H), `C4`, `D4`, `D8` (2n×2n, n=max(W,H)). Symmetric soups favour symmetric oscillators/ships that C1 soups rarely make |
| `growth 0\|1` | detect puffers/rakes/guns in soups (default 1). `growth 0` reproduces the old engine exactly |
| `probe 0\|1` | explore: run the probe on every feasible child and give novelty credit for new catalysts / signals / fossils (slower: ~1–2 s per rule) |
| `probe_ships N` `probe_stills N` | probe library size: most common ships / still lifes (4 / 16) |
| `probe_gens G` `probe_grid G` | generation limit of one collision (600), grid for probe runs (512) |
| `wick_period K` `wick_cells N` `wick_max N` | wick tiles: period ≤ K (6), ≤ N cells per tile (8), at most N wicks probed (400) |
| `fossil_cells N` `fossil_rate r` `fossil_max N` | fossils: ≥ N cells (12), in ≤ r·soups (0.002), report at most N (12) |
| `gun_ships N` | a soup that never settles and emitted ≥ N copies of one ship is reported as a gun candidate (10) |

## Algorithm
1. **Exact constraints.** Each pattern is run until its state repeats; all neighbourhoods seen (bbox+1)
   are frozen. Repeat ⇒ behaviour holds forever. Free variables = remaining symmetry classes.
2. **Common random numbers.** Parent and mutant are scored on the same soups, so differences
   are not swamped by soup-to-soup noise. Soup sets are rotated every `reseed` iterations.
3. **Racing.** Mutants run on a 32-soup probe, then batches of N/8. Paired per-soup differences of a
   linearised objective are tested after each batch; a mutant clearly worse (mean + 2.5 SE < 0),
   too slow, or too explosive is dropped early. Typically 50–70% of mutants stop early.
4. **Usage-weighted mutation.** Transition-class usage is counted in the parent's soups; classes are
   picked with probability ∝ sqrt(usage) + small ε, and every mutant flips at least one class that
   actually fired. (A class that never fires in those soups cannot change their result.)
5. **climb**: late-acceptance hill climbing (history `lahc`) on the objective
   Σ w·log(rate) − penalties (thresholds, bad soups, cost).
6. **explore**: constrained novelty search. Population of feasible rules (cap POPCAP); parent chosen
   ∝ (1+credit)/sqrt(1+visits). Each feasible child's census is scanned for non-target ships and
   oscillators (period ≥ min_osc_period); pseudo-objects (independent sub-groups) are removed;
   objects get a canonical id (min over phases and allowed orientations). Credit: new behaviour
   (type, period, velocity) 5 for ships, 1+log2(P)/2 for oscillators; a new shape of a known
   behaviour 0.1. Ids persist in the archive file, so reruns only report new things.

## Soup behaviours (puffers, rakes, guns)
Soups that never settle used to be just "bad". Now, with `growth 1` (default):
- **linear growth**: population rises by exactly +d every P generations (P ≤ 300) over the last 4 periods, checked
  when the soup ends (explosion, gmax, cost limit) or first touches the grid border. Puffers, rakes, blocked guns.
- **emitter**: the soup state repeats while ships were deleted inside the repeating window (deleted escaping ships
  keep coming back) — a gun. Previously such a soup was counted as *stable*. Also: a soup that never settles and
  emitted ≥ `gun_ships` copies of one ship. Each new emitter report is re-run without ship deletion: a gun's
  population keeps growing (reported), while a larger moving object whose front part is the ship does not
  (reported as probably not a gun; such objects also inflate ship counts in the census).
These are logged as `# NEW SOUP BEHAVIOUR ...` with the whole soup as RLE (open it in Golly), archived by
(period, growth) / (ship velocity), listed by `census`, and earn explore credit (3 per new growth period, 10 per
new gun ship type, times `reward growth` / `reward gun`). In explore they are recorded for every fully evaluated
rule, feasible or not (rules with common puffers are usually over `max_bad`). Census/climb results are unchanged.

## Probe: catalysts, signals, fossils
`probe RULE` runs a census (collecting the final ash), then three active searches. In explore, `probe 1` does the
same for every feasible child. New items go to the log/archive like other discoveries, each with an RLE of the
setup to open in Golly.
- **Catalysts.** Library = the most common real ships (pseudo-objects removed) and still lifes. Every ship is fired
  at every still life in each of its orientations, along every lane that comes within 2 cells (one ship phase is
  enough, since a still life is time-invariant). A collision counts only if the still life is disturbed (the region
  around it differs from still life + ship evolving separately) and is later back exactly, with everything else
  identified. Outcomes: **eater** (nothing else left), **reflector N deg** (the same ship leaves in another
  direction), **converter** (a different ship leaves), **splitter** (several ships), **pass-through** (same ship, same
  direction, after a real interaction), **catalysed reaction** (still life restored, other debris made).
  Explore credit: reflector 8, converter 5, splitter 5, pass-through 4, catalysed 1.5, eater 1 (× `reward catalyst`).
  Check: in Life it finds eater 1 and its 9-cell variant eating gliders and LWSSes.
- **Wicks, signals, fuses.** Wicks are stable periodic lines made of a small tile repeated along (k,0) or (k,k), whose
  pieces are not still lifes on their own (rows of independent still lifes are rejected). Each cell of the tile in
  the middle of a long wick is toggled and the damage (cells differing from the intact wick) is followed for 64
  gens inside the light cone of the wick's ends. **signal** = bounded damage whose front moves with an exact period
  and changes wick cells; **fuse** = the front moves the same way but the damaged part keeps growing (the wick is
  consumed). Credit: signal 6, fuse 2 (× `reward signal`). Check: in Life the diagonal line of cells is found as a
  speed-c fuse.
- **Fossils** (my reading of "index fossils"): large isolated objects in the final ash (≥ `fossil_cells`) that occur
  in few soups — leftovers that mark rare reactions. Pseudo-objects are excluded.
Not done yet (ideas): oscillators as catalysts (needs phase enumeration), ship–ship collisions (syntheses),
sparkers, signals on oscillating wicks, catalysts for soup reactions other than ship impacts.

## Census details
- Soups 16×16 at density 1/2 by default (bit-identical to the earlier census.cpp seeds).
- Ships are counted when isolated and nothing lies ahead in a cone along their direction; then removed.
- Oscillators are counted at the end if nothing lies within 2 cells of their full-period envelope.
- Soup ends when the hash repeats with period ≤ 60 twice, at `gmax`, or at `pmax` population.
- Classification of unknown final clusters goes to `maxp`; larger periods are missed.

## Known limits
- Escaping-ship test is a heuristic: a ship counted as escaping can in principle be caught by a faster one.
- Pseudo-object test: some partition (≤ 7 8-connected components) into groups that each evolve independently AND are each periodic; components that are mere sparks do not make an object pseudo.
- `split` checks for interaction only up to the longest component's repeat time.
- B0 rules are not supported (B0 is always forced off).

## Performance notes (2026-10)
All of these give bit-identical results to the previous engine (same config + seed, any thread count);
checked on census, climb and explore runs.
- Cluster classification runs on a local byte grid instead of `std::set`/`vector` cell lists, and stops as soon as
  the evolution enters a cycle that does not contain the starting state (it then can never recur).
- Before the final scan only ships matter, so a cluster is not classified at all when every one of the 8 possible
  directions of travel is already blocked by other debris (the escape-cone test is the same one used afterwards).
- On ARM (NEON) the soup and classification sweeps evaluate 16 cells per step: the 512-entry table is packed into a
  64-byte bit table, one `tbl` lookup on the NW..E bits, one per-lane shift on the SW,S,SE bits.
- Transition usage (mutation weights) is computed lazily: only for rules that are actually mutated from, by
  re-running their first 64 soups with counting on. Before, every evaluation ran its first 64 soups on the slow
  counting path, which in `climb` was the whole 32-soup race probe.
- `check`/constraint setup records the table entries a pattern uses in a 512-entry array (it used to store one
  entry per cell per generation, which took >16 GB on a chaotic non-repeating pattern), and stops with an error
  when a pattern grows past `pmax` cells.
- `climb` reports the best *feasible* rule (it used to report the highest f even when infeasible).
- Discoveries are merged in soup order and processed in key order, so the reported phase / first-seen seed of an
  object no longer depends on thread scheduling.
Measured on an M4 Pro (rule of `examples/diag_ship_explore.cfg`): census of 2000 soups 2.45 s -> 0.42 s single
thread; 60 climb iterations 167 s -> 21 s (1 thread), 34 s -> 3.4 s (14 threads). See `speedups.txt` for what was
tried and rejected.
