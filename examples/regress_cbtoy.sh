#!/bin/bash
# Regression harness for cb_davidson.x: the eigenvalues printed for several inputs, batch sizes and the overlap mode
# must stay identical (byte for byte) to the baseline. usage: regress_cbtoy.sh <bin_dir> <work_dir> [save|check]
BIN=$1; W=$2; MODE=${3:-check}; REPEAT=${REPEAT:-2}
HERE=$(cd "$(dirname "$0")" && pwd); mkdir -p "$W/base" "$W/run"; fail=0
extract() { grep -E "^ *k = |^ +-?[0-9]+\.[0-9]+( +-?[0-9]+\.[0-9]+)* *$" "$1"; }
mk() { # input, nk_batches, overlap(0/1) -> modified input on stdout
  python3 - "$1" "$2" "$3" <<'P'
import sys,re
f,nb,ov=sys.argv[1:4]
l=open(f).read().split('\n')
s=l[0]
s=re.sub(r'nk_batches *= *\d+','',s); s=re.sub(r'use_overlap *= *\S+','',s)
s=re.sub(r'david *= *(\d+) *,?',r'david=\1, ',s)
s=s.rstrip().rstrip('/').rstrip().rstrip(',')
s+=', nk_batches=%s'%nb+(', use_overlap=.true.' if ov=='1' else '')+' /'
print(s); print('\n'.join(l[1:]),end='')
P
}
cases="si2_11_points.in:1:0 si2_11_points.in:2:0 si2_11_points.in:4:0 si2_11_points.in:11:0 si2_11_points.in:4:1 si2_11_points.in:3:1 si2.in:2:0 si2.in.use_overlap:2:1 gaas2.in:1:0"
for c in $cases; do
  IFS=: read f nb ov <<< "$c"; tag=$(echo "$f" | tr . _)_${nb}_$ov
  for r in $(seq 1 $REPEAT); do
    mk "$HERE/$f" $nb $ov > "$W/in_$tag.txt"
    $BIN/cb_davidson.x < "$W/in_$tag.txt" > "$W/run/${tag}_$r.log" 2>&1 || { echo "FAIL (run error) $tag run $r"; fail=1; continue; }
    extract "$W/run/${tag}_$r.log" > "$W/run/${tag}_$r.txt"
    [ -s "$W/run/${tag}_$r.txt" ] || { echo "FAIL (no eigenvalues) $tag"; fail=1; continue; }
    if [ "$MODE" = save ]; then cp "$W/run/${tag}_$r.txt" "$W/base/$tag.txt"; break; fi
    if cmp -s "$W/run/${tag}_$r.txt" "$W/base/$tag.txt"; then s=identical; else s=DIFFERENT; fail=1; fi
    echo "$tag run $r: $s"
  done
done
[ "$MODE" = save ] && echo "baseline saved in $W/base ($(ls $W/base | wc -l) cases)" || { [ $fail = 0 ] && echo "REGRESSION cbToy: PASSED" || echo "REGRESSION cbToy: FAILED"; }
exit $fail
