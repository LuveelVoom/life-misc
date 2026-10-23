# rulesearch

Constrained search over range-1 Moore rules (Life-like, isotropic non-totalistic, or MAP).
Constraints are *exact*: patterns must evolve exactly as in a reference rule. Within those
constraints it either climbs toward making target objects common (`climb`) or looks for
new objects while keeping targets common (`explore`).

    g++ -O2 -march=native -fopenmp -std=c++17 -o rulesearch rulesearch.cpp
    ./rulesearch CONFIG check [RULE...]       # forced transitions, pattern periods, test rules
    ./rulesearch CONFIG census RULE [N] [SEED]
    ./rulesearch CONFIG climb ITERS
    ./rulesearch CONFIG explore ITERS [POPCAP]
    ./rulesearch CONFIG dump RULE SEED        # soup as RLE
    OMP_NUM_THREADS=n sets threads; results do not depend on thread count.

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
| `pattern NAME RLE [split]` | constraint: evolve until the state repeats (translation allowed) or dies; every transition used is frozen. `split` makes one constraint per component (components = groups 3+ cells apart), named NAME.0, NAME.1 … in order of top-left corner; it aborts if components interact. |
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
| `log FILE` `archive FILE` | discoveries log (RLE) and persistent object archive (one id per line) |

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
