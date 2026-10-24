# transcensus-int: development notes

How the port differs from the original, and what the tests check. Usage is in [README.md](README.md).

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

## Notes

- `--evolve` and `--emit-c` accept non-isotropic MAP rules (the census itself rejects them); they are used by the tests.
- Speed table in README.md measured 2026-10-22 on an otherwise idle M4 Pro.
- **Dependency to replace later:** `tests/test_rules.py` and `tests/test_engine.py` use casim (a separate Python simulator
  that follows Golly, not included in this repository) as the reference parser and simulator. Without it only
  `tests/kat.sh` runs. Replacing it needs another independent Golly-compatible INT/MAP implementation.
