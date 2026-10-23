#!/bin/sh
# Known-answer tests for transcensus-int (about 30 s on an M4 Pro). Run from anywhere: sh tests/kat.sh
# 1. B3/S23: byte-identical stdout and stderr to the unmodified original, built from orig/ with smaller constants
#    (200000 settled soups, report threshold 20) in a temp dir. The original's AVX2 warning lines are dropped.
# 2. INT rules: stdout+stderr digests recorded 2026-10-22, at -j 1 and -j 10, and from the fixed-rule build.
set -e
cd "$(dirname "$0")/.."
make -s transcensus-int
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
fail=0

sed -e 's/#define USED_SOUP_CNT ((s64) 10000000)/#define USED_SOUP_CNT ((s64) 200000)/' \
    -e 's/#define REPORT_THRESHOLD 200/#define REPORT_THRESHOLD 20/' orig/transcensus.c > "$T/tc_small.c"
cp -R orig/lib "$T/lib"
clang -std=c99 -O3 -Wno-unused-function "$T/tc_small.c" -o "$T/tc-orig-small" -lm
"$T/tc-orig-small" > "$T/orig.out" 2> "$T/orig.err.raw"
grep -v -e '^Warning: ' -e '^         ' "$T/orig.err.raw" > "$T/orig.err" || true
./transcensus-int -j 10 -n 200000 -t 20 > "$T/int.out" 2> "$T/int.err"
if cmp -s "$T/orig.out" "$T/int.out" && cmp -s "$T/orig.err" "$T/int.err"; then
    echo "ok   B3/S23 identical to the original"
else
    echo "FAIL B3/S23 differs from the original"; fail=1
fi

check() {   # rule, expected digest, executable, threads
    got=$("$3" -j "$4" -r "$1" -n 20000 -t 20 2>&1 | shasum | cut -d' ' -f1)
    if [ "$got" = "$2" ]; then echo "ok   $1 ($3 -j $4)"; else echo "FAIL $1 ($3 -j $4): $got"; fail=1; fi
}
check B3-cnqy4z/S23-a4i 213f7a33e51ff189ee0438bcb16f35ffbc1c3709 ./transcensus-int 1
check B3-cnqy4z/S23-a4i 213f7a33e51ff189ee0438bcb16f35ffbc1c3709 ./transcensus-int 10
check B3/S23-a4i 7ec3625d6c4d1d24ceabb7704b8a9cae99a4a37a ./transcensus-int 1
check B3/S23-a4i 7ec3625d6c4d1d24ceabb7704b8a9cae99a4a37a ./transcensus-int 10
make -s fixed RULE=B3-cnqy4z/S23-a4i > /dev/null
check B3-cnqy4z/S23-a4i 213f7a33e51ff189ee0438bcb16f35ffbc1c3709 ./transcensus-int-B3-cnqy4z_S23-a4i 10

exit $fail
