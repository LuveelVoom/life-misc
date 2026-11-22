#!/bin/sh
# 1) watch qfind (PID 91891) for its first ship; save it and kill qfind
# 2) run cf_batch.py (30 min per search)  3) Barrister lone-object searches, then casim checks
P="$(cd "$(dirname "$0")" && pwd)"
LT="$(cd "$P/../../../../.." && pwd)"          # .../Life tools
QPID=91891
SHIP="$LT/qfind_c7_w9_ship.rle"
log(){ echo "$(date '+%F %T') $*"; }
log "pipeline start; waiting for qfind $QPID"
while kill -0 $QPID 2>/dev/null; do
  osascript "$P/readhist.applescript" > "$P/qfind_history.txt" 2>"$P/osascript.err" || log "osascript failed: $(cat "$P/osascript.err")"
  if python3 "$P/extract_ship.py" "$P/qfind_history.txt" "$SHIP"; then
    log "ship found, saved to $SHIP; killing qfind"; kill $QPID; break
  fi
  sleep 60
done
while kill -0 $QPID 2>/dev/null; do sleep 5; done
osascript "$P/readhist.applescript" > "$P/qfind_history_final.txt" 2>/dev/null
if [ ! -f "$SHIP" ]; then
  if python3 "$P/extract_ship.py" "$P/qfind_history_final.txt" "$SHIP"; then log "ship found in final output, saved"; else log "qfind ended without a ship"; fi
fi
log "starting cf_batch"
cd "$LT" && python3 cf_batch.py --catforce catforce-int/CatForce/CatForce \
  --transcensus "Vibecoded life tools/life-misc/transcensus/transcensus-int" --out runs3 --max-time 30m > "$LT/cf_batch_runs3.log" 2>&1
log "cf_batch exit $?"
cd "$P"
log "starting Barrister searches"
for t in piq4i_N piq4i_E a5_N a5_E a5_W; do
  ( perl -e 'alarm 1800; exec @ARGV' ../../barrister-int/Barrister $t.toml > $t.out 2>&1; echo "done $? $(date +%T)" >> $t.out ) &
done
wait
log "Barrister searches done"
python3 lonecheck.py B3/S23-q4i piq4i_N.base.rle pi_q4i_summary.rle piq4i_N.out piq4i_E.out > pi_q4i_summary.txt 2>&1
python3 lonecheck.py B3/S23-a5 a5_N.base.rle a5_summary.rle a5_N.out a5_E.out a5_W.out > a5_summary.txt 2>&1
cat pi_q4i_summary.txt a5_summary.txt
log "pipeline finished"
