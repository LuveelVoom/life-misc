# transcensus-int

Port of `transcensus.c` (a 2016 active-region soup census, the source in `../tc.zip`, unpacked unmodified in
`orig/`) from B3/S23 to Life-like, isotropic non-totalistic (Hensel/INT) and isotropic MAP rules.

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

## Speed (M4 Pro, measured 2026-10-22, machine otherwise idle)

| run | time |
|---|---|
| original, B3/S23, 200k settled, 1 thread | 13.5 s |
| generic build, B3/S23, 1 thread | 44.9 s |
| generic build, B3/S23, `-j 10` | 6.6 s |
| `make fixed`, B3/S23, 1 thread | 14.8 s |
| generic / fixed, B3-cnqy4z/S23-a4i, 20k settled, 1 thread | 3.2 s / 0.9 s |

Single runs, so expect a few percent of noise. The full default census (10M settled) takes about 11 minutes
single-threaded for Life with the original or the fixed build.

**Watch the settle rate** in the progress lines (`Testing soup N, M used`). In rules whose soups rarely drop
to ≤ 12 cells within 256 generations, almost nothing settles. B2-a3/S23 settled 1,325 of 18.6M soups in
22 minutes on 10 threads. For such rules, change `-a/-b/-g` or pick another rule; the defaults are tuned for Life.

## Rules

- Golly syntax (`B3/S23`, `b3s23`, `23/3`, `B2-a3/S23`, `B2-a3_S23`, `MAP...`). The parser is a C port of
  casim's, which follows Golly's `liferules.cpp`. Checked equal to casim on 225 rule strings (`tests/test_rules.py`).
- **B0 rules are rejected.** The grid has a fixed dead background and no strobing.
- **Non-isotropic MAP rules are rejected**, because the census merges all 8 orientations. `--evolve` and
  `--emit-c` still accept them (for tests).
- Bounded-grid suffixes (`:T...`), other neighbourhoods and multistate rules are rejected.

## How the port works

- **Changed:** the only change to the library is the evolve step. `lib/golgrid.c` is `orig/lib/golgrid.c`
  except for the four `GoLGrid_int_evolve_*` loops (`diff orig/lib/golgrid.c lib/golgrid.c`). Strip layout,
  edge masking and bounding-box tracking are untouched. The other `lib/` files are copied unchanged.
- **Rule evaluation:** `lib/rule.c` turns the 512-entry rule table into a reduced ordered BDD. At start-up it
  tries 2000 variable orders and keeps the smallest (Life 26 nodes, typical INT rules 40–80). The BDD is
  evaluated bitsliced, one node at a time over 32 rows (vectorised by clang).
  `make fixed` writes the BDD as straight-line C instead (`--emit-c`), which is close to the original's
  hand-written Life adder.
- **Bit order:** in a grid word bit 63 is the leftmost cell, so `w >> 1` gives each cell's west neighbour.
  Rule table index is Golly MAP order (NW 256 … SE 1).
- **Threads:** the main thread draws soups from the generator in the original order, chunks of 4096 run in
  parallel, and results are merged in soup order. The output does not depend on `-j` and equals the
  single-threaded original.
- **Other changes:**
  - A key pattern store overflow now stops the run and still reports. The original exited with no output.
    For chaotic rules raise `-k` (each slot costs about 1.2 KB, touched only when used).
  - The output grid grows past 16384 rows only if more than 4080 patterns are reported. The original
    silently clipped them.
  - The x86 `verify_cpu_type` check is dropped.

## Tests (`tests/`)

- `test_rules.py`: rule parsing against casim.
- `test_engine.py [--generated]`: `--evolve` against casim on random soups. Covers Life, totalistic, INT,
  a random non-isotropic MAP rule and a pure NE-shift rule (orientation), on and across word boundaries.
  `--generated` repeats it for the `make fixed` path. A planted NW/NE swap fails 74 of 108 cases.
- `kat.sh`: known answers. B3/S23 must equal the unmodified original byte for byte (200k soups, built from
  `orig/` with smaller constants in a temp dir). The INT runs must give the recorded output digests with
  -j 1, -j 10 and the fixed-rule build.
