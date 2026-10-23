# transcensus-int

A soup census of a rule's small *precursors*: patterns of at most 12 cells that turn up in soups and then grow
into an active reaction. Results are grouped by the reaction they grow into, most common first. It is a port of `transcensus.c`
(a 2016 census for B3/S23; the original source is in `../tc.zip`, unpacked unmodified in `orig/`) to Life-like,
isotropic non-totalistic (Hensel/INT) and isotropic MAP rules. For B3/S23 the output is identical to the
original's, and it can use all your cores.

## What the census does (unchanged from the original)

Each soup is a random 8×8 block on a 192×192 grid with a dead border (cells outside the grid are always dead).
The population is followed for up to 256 generations:

1. Once the population is ≤ 12 (`-a`), the smallest-population state seen is kept as the **min pattern**.
2. When it then grows to ≥ 20 (`-b`), it counts as active again. If it falls back to ≤ 12 first, it starts over.
3. Once active, when 4 (`-c`) different populations below 96 have been seen, the current pattern is the
   **key pattern** and the soup counts as *settled* (used).

Key patterns are merged under translation and the 8 symmetries (by a 64-bit hash, as in the original).
Each one keeps up to 16 different min patterns (precursors). After 10M settled soups (`-n`), every key pattern
seen ≥ 200 times (`-t`) is printed with its best precursor: the precursor with the highest count × size
weight (weight 128 for ≤ 4 cells down to 3 for 12 cells).

## Build and run

```
make                           # transcensus-int (Apple clang + Homebrew libomp)
./transcensus-int -r B2-a3/S23 -j 10 > out.rle 2> log.txt
make fixed RULE=B2-a3/S23      # transcensus-int-B2-a3_S23: about 3x faster, only runs that rule
./transcensus-int --help
```

Options: `-r` rule, `-n` settled soups, `-t` report threshold, `-g` generations, `-a/-b/-c` the thresholds
above, `-k` key pattern store size, `-s` seed, `-j` threads, `-l FILE` plain list.
Defaults are the original's constants.

Output on stdout is the original's: a `Key pattern N, count = C:` line per reported pattern, then one
LifeHistory RLE with the precursors in state 1 and their counts written above them in state 4.
For any rule other than Life, `rule = LifeHistory` there is only for display: running the RLE in Golly would
run Life. Use `-l FILE` for a list with one line per reported pattern (count, tab, the precursor as one-line
RLE body), which you can paste into Golly under the census rule.

## Speed (M4 Pro)

| run | time |
|---|---|
| original, B3/S23, 200k settled, 1 thread | 13.5 s |
| generic build, B3/S23, 1 thread | 44.9 s |
| generic build, B3/S23, `-j 10` | 6.6 s |
| `make fixed`, B3/S23, 1 thread | 14.8 s |
| generic / fixed, B3-cnqy4z/S23-a4i, 20k settled, 1 thread | 3.2 s / 0.9 s |

Single runs, so expect a few percent of noise. `make fixed` builds a version for one rule only, about 3× faster
than the generic build. The full default census (10M settled) takes about 11 minutes
single-threaded for Life with the original or the fixed build.

**Watch the settle rate** in the progress lines (`Testing soup N, M used`). In rules whose soups rarely drop
to ≤ 12 cells within 256 generations, almost nothing settles. B2-a3/S23 settled 1,325 of 18.6M soups in
22 minutes on 10 threads. For such rules, change `-a/-b/-g` or pick another rule; the defaults are tuned for Life.

## Rules

- Golly syntax (`B3/S23`, `b3s23`, `23/3`, `B2-a3/S23`, `B2-a3_S23`, `MAP...`). The parser is a C port of
  casim's, which follows Golly's `liferules.cpp`. Checked equal to casim on 225 rule strings (`tests/test_rules.py`).
- **B0 rules are rejected.** The grid has a fixed dead background and no strobing.
- **Non-isotropic MAP rules are rejected**, because the census merges all 8 orientations.
- Bounded-grid suffixes (`:T...`), other neighbourhoods and multistate rules are rejected.
- Threads don't change the output: `-j 10` gives exactly what `-j 1` gives.
- In chaotic rules the key pattern store can fill up. The run then stops early but still reports; raise `-k`
  (each slot costs about 1.2 KB, only when used).

## Tests

`tests/` checks rule parsing and the simulation against casim, and `tests/kat.sh` checks known census outputs
(including byte-identical output to the original for B3/S23). Implementation notes and test details are in
[DEVELOPMENT.md](DEVELOPMENT.md).
