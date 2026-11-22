#!/bin/sh
# CEGAR sweep for B3/S23-a5: 12 background sfcegar.py jobs, one per configuration.
# usage: sh sweep.sh OUTDIR [extra sfcegar options, e.g. --rigid 2]
# Logs: OUTDIR/c_<box>_<inset>.log; a hit is written to OUTDIR/found_<box>_<inset>.txt.
cd "$(dirname "$0")"
O=${1:?usage: sh sweep.sh OUTDIR [sfcegar options]}; shift
mkdir -p "$O"
R=B3/S23-a5
A=.oooo../.oo..o./o.o..oo/...oooo/o.ooo../oo..o../oo...oo/
run() {  # name, then sfcegar arguments
  n=$1; shift
  nohup python3 sfcegar.py "$R" "$@" --out "$O/found_$n.txt" > "$O/c_$n.log" 2>&1 &
}
run 12_0 12 12 --inset 0 --minpop 48  "$@"
run 16_0 16 16 --inset 0 --minpop 85  "$@"
run 12_1 12 12 --inset 1 --minpop 48  "$@"
run 14_1 14 14 --inset 1 --minpop 65  "$@"
run 16_1 16 16 --inset 1 --minpop 85  "$@"
run 20_1 20 20 --inset 1 --minpop 133 "$@"
run 14_2 14 14 --inset 2 --minpop 65  "$@"
run 16_2 16 16 --inset 2 --minpop 85  "$@"
run 20_2 20 20 --inset 2 --minpop 133 "$@"
run 24_2 24 24 --inset 2 --minpop 192 "$@"
# 7x7 agar fixed on the central 14x14, walls free
run agar24_1 24 24 --inset 1 --agar $A --core 5 5 19 19 "$@"
run agar24_2 24 24 --inset 2 --agar $A --core 5 5 19 19 "$@"
