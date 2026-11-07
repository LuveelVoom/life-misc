#!/bin/sh
# Known-answer test, a few seconds. test/period.py checks the patterns found.
cd "$(dirname "$0")/.." || exit 1
ROOT=$PWD
OUT=$(mktemp -d)
python3 test/wrap_agreement.py "$ROOT"
check() { # file rule extra-args...   prints the period; compare with the expected line below
    f=$1; r=$2; shift 2
    (cd lls && python3 lls "$ROOT/test/$f" --rule "$r" -S kissat --population '>0' "$@" -o "$OUT/$f.rle" -v 0)
    echo "$f $r: $(python3 test/period.py "$OUT/$f.rle" "$r")"
}
check glider_nw.jdf B3/S23            # expected: spaceship, period 4, dx -1, dy -1
check flip_vertical.jdf B3/S23        # expected: spaceship, period 4, dx 2, dy 0 (p2 glide reflection)
check p3.jdf B3/S23-a4 --force_change 0 1   # expected: oscillator, period 3
(python3 backtrack.py test/glider_nw.jdf --threads 2 --rule B3/S23 --population '>0' -o "$OUT/bt.rle" >/dev/null)
echo "backtracker: $(python3 test/period.py "$OUT/bt.rle" B3/S23)"   # expected: spaceship, period 4, dx -1, dy -1
rm -rf "$OUT"
