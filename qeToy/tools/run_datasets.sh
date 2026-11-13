#!/bin/bash
# Runs qe_davidson.x on the dataset package (build_marco/qeToy_datasets) and compares with the stored pw.x results (no QE needed):
# eigenvalue differences to the CPU and GPU pw.x xml files, Davidson time, steady-state ms per k-point, and the pw.x cegterg times.
# usage: run_datasets.sh <bin_dir> <datasets_dir> [slots="16 32"] [david="2 4"] [cases="al_nc_24 al_us_8 AlNC40 CuUS CuPAW FeLSDA"]
BIN=$1; D=$2; SLOTS=${3:-"16 32"}; DAV=${4:-"2 4"}; CASES=${5:-"al_nc_24 al_us_8 AlNC40 CuUS CuPAW FeLSDA"}
HERE=$(cd "$(dirname "$0")" && pwd); W=${WORK:-$(mktemp -d)}; mkdir -p "$W"
for c in $CASES; do
  . "$D/$c/case.conf"
  echo "== $c: $DESC (ethr $ETHR)"
  for f in cpu gpu; do [ -f "$D/$c/pw_${f}_nscf.out" ] && echo "   pw.x $f build: cegterg $(grep -E '^ +cegterg +:' $D/$c/pw_${f}_nscf.out | head -1 | awk '{for(i=1;i<=NF;i++) if($i=="WALL") print $(i-1)}') , c_bands $(grep -E '^ +c_bands +:' $D/$c/pw_${f}_nscf.out | head -1 | awk '{for(i=1;i<=NF;i++) if($i=="WALL") print $(i-1)}') (WALL, recorded on the GB10 machine)"; done
  for dv in $DAV; do for nb in $SLOTS; do
    out=$W/${c}_${dv}_$nb.txt
    printf " &qe\n dumpdir='%s', nk_batches=%s, david=%s, ethr=%s, out='%s'\n /\n" "$D/$c/dump" "$nb" "$dv" "$ETHR" "$out" > $W/in.txt
    log=$($BIN/qe_davidson.x < $W/in.txt 2>&1)
    tt=$(echo "$log" | grep "Davidson time" | sed 's/Davidson time *\([0-9.]*\) s.*/\1/')
    med=$(echo "$log" | grep -E "^group +[0-9]" | python3 -c "
import sys,re,statistics as st
t=[float(re.search(r'time\s+([\d.]+)',l).group(1)) for l in sys.stdin]
n=$nb
print('%.2f'%(1e3*st.median(t[1:])/n) if len(t)>1 else 'n/a')")
    dc=$(python3 $HERE/compare_xml.py $out $D/$c/ref_cpu.xml 2>&1 | tail -1 | sed 's/.*max |diff| = //; s/ eV.*//')
    dg="-"; [ -f $D/$c/ref_gpu.xml ] && dg=$(python3 $HERE/compare_xml.py $out $D/$c/ref_gpu.xml 2>&1 | tail -1 | sed 's/.*max |diff| = //; s/ eV.*//')
    printf "   david=%s slots=%3s: Davidson %7s s, steady %6s ms/k, max|diff| vs CPU pw.x %s eV, vs GPU pw.x %s eV\n" $dv $nb "$tt" "$med" "$dc" "$dg"
  done; done
done
