#!/bin/bash
# Comparison of qeToy with pw.x (CPU and GPU builds) on a metal: scf, nscf with both pw.x, dump, qeToy at several group sizes.
# usage: compare_pw.sh <work_dir> <name> <element> <pseudo> <lattice: fcc|bcc> <alat_bohr> <ecutwfc> <ecutrho> <nbnd> <k_scf> <k_nscf> [nspin=1] [starting_magnetization]
# needs: the environment of the QE build (compilers, mpirun), QE_CPU, QE_GPU (pw.x), DUMPER (dump_ham.x), QETOY (qe_davidson.x), PSEUDO_DIR
# prints a table with the cegterg and c_bands times of pw.x and the Davidson time of qeToy, and the eigenvalue differences to both xml files
set -e
W=$1; NAME=$2; EL=$3; PS=$4; LAT=$5; A=$6; EW=$7; ER=$8; NB=$9; KS=${10}; KN=${11}; NSPIN=${12:-1}; MAG=${13:-0.5}
HERE=$(cd "$(dirname "$0")" && pwd)
: "${QE_CPU:?}" "${QE_GPU:?}" "${DUMPER:?}" "${QETOY:?}" "${PSEUDO_DIR:?}"
IBRAV=2; [ "$LAT" = bcc ] && IBRAV=3
MP=1; [ "$NSPIN" = 2 ] && MP=2
mkdir -p "$W/$NAME"; cd "$W/$NAME"
mkin() { # calculation k outfile
cat > $3 <<EOF
 &control
    calculation='$1', prefix='$NAME', outdir='./tmp', pseudo_dir='$PSEUDO_DIR'
 /
 &system
    ibrav=$IBRAV, celldm(1)=$A, nat=1, ntyp=1, ecutwfc=$EW, ecutrho=$ER, nbnd=$NB,
    occupations='smearing', smearing='marzari-vanderbilt', degauss=0.05$( [ "$NSPIN" = 2 ] && echo ", nspin=2, starting_magnetization(1)=$MAG" )
 /
 &electrons
$( [ "$1" = scf ] && echo "    conv_thr=1.d-10" )
$( [ "$1" = nscf ] && [ -n "$EXTRA_NSCF" ] && echo "    $EXTRA_NSCF" )
 /
ATOMIC_SPECIES
 $EL 1.0 $PS
ATOMIC_POSITIONS (alat)
 $EL 0.0 0.0 0.0
K_POINTS automatic
 $2 $2 $2 0 0 0
EOF
}
mkin scf $KS scf.in; mkin nscf $KN nscf.in
sed -i "s/nbnd=$NB,/nbnd=$NB,/" scf.in
t() { /usr/bin/time -f "%e" "$@" 2>&1 >/dev/null | tail -1; }
done_ok() { [ -s "$1" ] && grep -q "JOB DONE" "$1"; }
if done_ok scf.out && [ -z "$FORCE" ]; then echo "[$NAME] scf: reusing"; else echo "[$NAME] scf (CPU pw.x)"; rm -rf tmp cpu gpu; mpirun --allow-run-as-root -np 1 $QE_CPU -i scf.in > scf.out 2>&1; fi
for v in cpu gpu; do [ -d $v/tmp ] || { mkdir -p $v; cp -r tmp $v/tmp; }; cp nscf.in $v/; done
if done_ok cpu/nscf.out && [ -z "$FORCE" ]; then echo "[$NAME] nscf CPU: reusing"; else echo "[$NAME] nscf CPU"; (cd cpu && mpirun --allow-run-as-root -np 1 $QE_CPU -i nscf.in > nscf.out 2>&1); fi
if done_ok gpu/nscf.out && [ -z "$FORCE" ]; then echo "[$NAME] nscf GPU: reusing"; else echo "[$NAME] nscf GPU"; (cd gpu && mpirun --allow-run-as-root -np 1 $QE_GPU -i nscf.in > nscf.out 2>&1); fi
ETHR=$(grep -m1 "ethr =" cpu/nscf.out | sed 's/.*ethr = *\([0-9.E+-]*\),.*/\1/')
NK=$(grep "number of k points" cpu/nscf.out | head -1 | sed 's/.*points= *\([0-9]*\).*/\1/')
echo "[$NAME] dump, ethr=$ETHR, k-points=$NK"
rm -rf dump; printf " &inputpp\n  outdir='./cpu/tmp', prefix='$NAME', dumpdir='./dump', dump_format='hdf5', dump_tests=.false.\n /\n" > dump.in
mpirun --allow-run-as-root -np 1 $DUMPER -i dump.in > dump.out 2>&1
grep -E "wrote" dump.out | tr -s ' '
pwt() { grep -E "^ +$2 +:" $1 | head -1 | sed 's/.*CPU *\([0-9.]*\)s WALL.*/\1/; s/.*: *[0-9.]*s CPU *\([0-9.]*\)s WALL.*/\1/'; }
tm() { grep -E "^ +$2 +:" $1 | head -1 | awk '{for(i=1;i<=NF;i++) if($i=="WALL") print $(i-1)}' | tr -d s; }
echo; echo "pw.x timings (WALL, s)   CPU build   GPU build"
for k in cegterg c_bands PWSCF; do printf "  %-10s %12s %12s\n" $k "$(tm cpu/nscf.out $k)" "$(tm gpu/nscf.out $k)"; done
AVG=$(grep -m1 "avg # of iterations" cpu/nscf.out | sed 's/.*iterations = *//')
echo "  Davidson avg iterations (CPU pw.x): $AVG"
echo; echo "qeToy (Davidson time only, ethr=$ETHR, all bands tight)"
for dv in 2 4; do for nb in 16 32; do
  printf " &qe\n dumpdir='./dump', nk_batches=$nb, david=$dv, ethr=$ETHR, out='eig_${dv}_$nb.txt'\n /\n" > in.txt
  out=$($QETOY < in.txt 2>&1 | grep -E "^group +[0-9]|Davidson time|more than one call")
  tt=$(echo "$out" | grep "Davidson time" | sed 's/Davidson time *\([0-9.]*\) s.*/\1/')
  med=$(echo "$out" | grep -E "^group +[0-9]" | python3 -c "
import sys,re,statistics as st
t=[float(re.search(r'time\s+([\d.]+)',l).group(1)) for l in sys.stdin]
n=$nb
print('%.2f'%(1e3*st.median(t[1:])/n) if len(t)>1 else 'n/a')")
  dc=$(python3 $HERE/compare_xml.py eig_${dv}_$nb.txt cpu/tmp/$NAME.save/data-file-schema.xml | sed 's/.*max |diff| = //; s/ eV.*//')
  dg=$(python3 $HERE/compare_xml.py eig_${dv}_$nb.txt gpu/tmp/$NAME.save/data-file-schema.xml | sed 's/.*max |diff| = //; s/ eV.*//')
  printf "  david=%s slots=%2s: Davidson %6s s (steady %5s ms/k), max|diff| vs CPU pw.x %s eV, vs GPU pw.x %s eV\n" $dv $nb "$tt" "$med" "$dc" "$dg"
done; done
echo "CPU vs GPU pw.x eigenvalues: $(python3 /home/qespresso/KS_SOLVER_CB/build_marco/qe_miniapp_startup/scripts/compare_pp_xml.py cpu/tmp/$NAME.save/data-file-schema.xml gpu/tmp/$NAME.save/data-file-schema.xml 2>&1 | tail -1)"
