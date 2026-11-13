# Dump format: data of a QE run for the batched miniapp

`dump_ham.x` (source: `PP/src/dump_ham.f90` in the Quantum ESPRESSO tree, derived from `pp_example.f90`) reads a finished QE run and writes, for every k-point, what is needed to
build H|psi> and S|psi> outside QE, plus QE's own results on test vectors. `read_dump_check.f90` rebuilds H|psi>, S|psi>, `h_diag`
and `s_diag` from the dumped data alone (QE-independent, gfortran) and compares them with QE's.

## Build and run
The dumper is part of Quantum ESPRESSO (`PP/src/dump_ham.f90`; target `dump_ham.x` in the PP Makefile and `qe_pp_dump_ham_exe` in PP/CMakeLists.txt).
```
cd <qe>/PP && make dump_ham.x                    # after building pw and the PP libraries
printf " &inputpp\n  outdir='./tmp', prefix='al', dumpdir='./dump_al'\n /\n" > dump.in
mpirun -np 1 <qe>/bin/dump_ham.x -i dump.in      # after pw.x scf + nscf with the same outdir/prefix
gfortran -O2 -ffree-line-length-none -o read_dump_check tools/read_dump_check.f90
./read_dump_check dump_al [k ...]                # default: first, middle and last k-point
```
(`mpirun` prints an error message after the final STOP; the exit code 1 is harmless.) Not supported: non-collinear and gamma-only runs.
Files are Fortran unformatted, `access='stream'`, native endianness, 32-bit integers, 64-bit reals.

## HDF5 variant (`dump_format='hdf5'`)
The dumper can write the same data as one file `<dumpdir>/ham.h5` (QE built with HDF5, e.g. `-DQE_ENABLE_HDF5=ON`; serial HDF5 is enough).
`qe_ham_module` uses `ham.h5` automatically when it is present in the dump directory (this repo built with `-DQE_ENABLE_HDF5=ON`).
Layout (source: `qe_lib/qe_dump_h5.f90`, the same file is in `PP/src` of QE): root attributes `magic`, `version`, `nspin`, `nkstot`, `nbnd`, `npwx`, `nkb`, `nat`, `nsp`, `nhm`,
`nr1..nr3`, `nr1x..nr3x`, `nnr`, `ngms`, `okvan`, `doublegrid`, `has_tests`, `v_of_0`, `omega`, `alat`, `at(9)`, `bg(9)`; root datasets `nl`, `vrs`, `ityp`, `nh`, `ofsbeta`,
`deeq` (stored as `(nhm,nhm,nat*nspin)`), `qq_at`, `xk`, `wk`, `isk`, `ngk`; one group `k00001`, `k00002`, ... per k-point with `igk`, `g2kin`, `h_diag`, `s_diag`, `vkb`, `et` and, if
`has_tests`, `evc`, `hpsi`, `spsi`, `hnl`. Complex arrays have a leading dimension 2 (real, imaginary). Contents are identical to the stream files (checked: same eigenvalues from the Davidson).
`tools/dump_to_h5.f90` converts an existing stream dump (`dump_to_h5.x indir outdir [notests]`).
Al (413 k-points): 38.7 MB with test vectors, 5.2 MB without.

## ham_global.bin (once)
| item | type and shape |
|---|---|
| `magic` (20251113), `nspin`, `nkstot`, `nbnd`, `npwx`, `nkb`, `nat`, `nsp`, `nhm`, `nr1,nr2,nr3`, `nr1x,nr2x,nr3x`, `nnr`, `ngms`, `okvan`, `doublegrid` | 17 int32 (smooth grid) |
| `v_of_0`, `omega`, `alat`, `at(3,3)`, `bg(3,3)` | real64 |
| `nl(ngms)` | int32: 1-based index in the smooth FFT array of G vector number ig (ig = `igk`) |
| `vrs(nnr,nspin)` | real64: local potential on the smooth grid (what QE's `h_psi` multiplies by) |
| `ityp(nat)`, `nh(nsp)`, `ofsbeta(nat)` | int32: atom types, projectors per type, offset of the projectors of each atom in `vkb` |
| `deeq(nhm,nhm,nat,nspin)` | real64: D_ij per atom and spin (= `dvan` for norm-conserving) |
| `qq_at(nhm,nhm,nat)` | real64: overlap integrals (zero for norm-conserving) |
| `xk(3,nkstot)`, `wk(nkstot)`, `isk(nkstot)`, `ngk(nkstot)` | real64, real64, int32, int32 |

## ham_k.bin (one block per k-point, in the order of the run)
`ik, npw` (int32); `igk(npw)` (int32); `g2kin(npw)`, `h_diag(npw)`, `s_diag(npw)` (real64); `vkb(npw,nkb)` (complex128);
`et(nbnd)` eigenvalues of the run in Ry; `evc(npw,nbnd)` converged wavefunctions (test vectors);
`hpsi(npw,nbnd)`, `spsi(npw,nbnd)` = QE's `h_psi`, `s_psi` of `evc`; `hnl(npw,nbnd)` = the non local part alone (`add_vuspsi`).

## The Hamiltonian in these terms
`H psi = g2kin*psi + IFFT-multiply-FFT(vrs, psi placed at nl(igk)) + vkb * D * (vkb^H psi)`,
`S psi = psi + vkb * qq * (vkb^H psi)`, with D and qq block diagonal per atom (`ofsbeta(na)+ih`).
Preconditioner (QE's `usnldiag`): `h_diag = g2kin + v_of_0 + sum D_ij Re(vkb_i vkb_j^*)`, `s_diag = 1 + sum qq_ij Re(vkb_i vkb_j^*)`.

## Verified
| run | k-points | smooth grid | `H psi`, `S psi`, `V_NL psi` vs QE | `h_diag`, `s_diag` vs QE |
|---|---|---|---|---|
| Al, `Al.pz-vbc` (norm-conserving), ecut 15, nscf 24x24x24 | 413 | 15^3 | 2.7e-15, 0, 1.7e-16 (max) | exact |
| Al, `Al.pbe-nl-rrkjus_psl` (ultrasoft), ecut 25/200, nscf 8x8x8 | 29 | 18^3 (dense 24^3) | 2.8e-15, 2.2e-16, 5.6e-16 | exact |
(checked on k-points 1, middle and last; the NC dump is 37 MB.)
