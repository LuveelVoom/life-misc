#!/bin/sh
# Known-answer test, a few seconds. Needs casim, for checking the patterns found.
cd "$(dirname "$0")/.." || exit 1
ROOT=$PWD
CASIM="${CASIM:-$HOME/Desktop/Life tools/Vibecoded life tools/casim/casim.py}"
OUT=$(mktemp -d)
python3 test/wrap_agreement.py "$ROOT"
check() { # file rule extra-args...   prints casim's verdict; compare with the expected line below
    f=$1; r=$2; shift 2
    (cd lls && python3 lls "$ROOT/test/$f" --rule "$r" -S kissat --population '>0' "$@" -o "$OUT/$f.rle" -v 0)
    echo "$f $r: $(python3 "$CASIM" period "$OUT/$f.rle")"
}
check glider_nw.jdf B3/S23            # expected: period 4, displacement (-1, -1), up-left
check flip_vertical.jdf B3/S23        # expected: period 4, displacement (2, 0); p2 glide reflection
check p3.jdf B3/S23-a4 --force_change 0 1   # expected: period 3
(python3 backtrack.py test/glider_nw.jdf --threads 2 --rule B3/S23 --population '>0' -o "$OUT/bt.rle" >/dev/null)
echo "backtracker: $(python3 "$CASIM" period "$OUT/bt.rle")"   # expected: period 4, dx -1, dy -1
rm -rf "$OUT"
