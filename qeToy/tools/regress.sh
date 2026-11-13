#!/bin/bash
# Regression harness for qe_davidson.x (correctness only, timing belongs to the target GPU).
# usage: regress.sh <bin_dir> <work_dir> <al_nc_dump> <al_us_dump> <al_nc_xml> <al_us_xml> [save|check]
#   save : store the eigenvalue files as the baseline (bit-for-bit reference)
#   check: run again, compare with the baseline (cmp) and with pw.x (compare_xml.py), repeat each case REPEAT times (default 3)
BIN=$1; W=$2; NC=$3; US=$4; XNC=$5; XUS=$6; MODE=${7:-check}; REPEAT=${REPEAT:-3}
HERE=$(cd "$(dirname "$0")" && pwd); mkdir -p "$W/base" "$W/run"; fail=0
cases="nc:1 nc:4 nc:16 us:1 us:4 us:8"
for c in $cases; do
  k=${c%%:*}; nb=${c##*:}; dump=$NC; xml=$XNC; ethr=1e-8; nks=0
  [ $k = us ] && { dump=$US; xml=$XUS; ethr=1e-9; }
  [ $k = nc ] && nks=0
  for r in $(seq 1 $REPEAT); do
    out=$W/run/${k}_${nb}_$r.txt
    printf " &qe\n dumpdir='%s', nk_batches=%s, david=4, ethr=%s, out='%s'\n /\n" "$dump" "$nb" "$ethr" "$out" > "$W/in.txt"
    $BIN/qe_davidson.x < "$W/in.txt" > "$W/run/${k}_${nb}_$r.log" 2>&1 || { echo "FAIL (run error) $k nk_batches=$nb run $r"; fail=1; continue; }
    if [ "$MODE" = save ]; then cp "$out" "$W/base/${k}_${nb}.txt"; break; fi
    if cmp -s "$out" "$W/base/${k}_${nb}.txt"; then s=identical; else s=DIFFERENT; fail=1; fi
    d=$(python3 "$HERE/compare_xml.py" "$out" "$xml" | sed 's/.*max |diff| = //; s/ eV.*//')
    echo "$k nk_batches=$nb run $r: vs baseline $s, vs pw.x $d eV"
  done
done
[ "$MODE" = save ] && echo "baseline saved in $W/base" || { [ $fail = 0 ] && echo "REGRESSION: PASSED" || echo "REGRESSION: FAILED"; }
exit $fail
