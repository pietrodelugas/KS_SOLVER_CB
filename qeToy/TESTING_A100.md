# Testing and debugging qeToy and the single-thread cegterg on an A100

Hand-over written on 13 Nov 2026 after the work on the GB10 development machine (20 CPU cores, small GPU).
Everything below that is stated as measured was measured there; nothing has been run on an A100 yet.

## 1. What exists and where

| item | where |
|---|---|
| repository | `KS_SOLVER_CB`, local tags `v1.0` (cbToy before the QE miniapp, 9725d7c), `v1.1` (+ qeToy, one OpenMP thread per k-point, 0c3fd73), `v2.0` (single-thread `cegterg`, 49db2b0) |
| working branch | `single-thread-cegterg` (HEAD 882e3dc). `v1.1` is also on `fix-batched-cegterg-multibatch`. Nothing is pushed. |
| QE side | `~/q-e`, branch `dump-ham-pp` (2c393186e, b97536597): `PP/src/dump_ham.f90`, `PP/src/qe_dump_h5.f90`, Makefile and CMake entries. Not pushed. |
| miniapp | `qeToy/` (qe_davidson.x), toy miniapp `cbToy/` (cb_davidson.x), solver `KS_Solvers/Davidson/cegterg.f90` |
| datasets for testing without QE | `build_marco/qeToy_datasets/` (57 MB, git-ignored) and `build_marco/qeToy_datasets.tar.gz`; copy it to the A100 machine |
| documentation | `qeToy/DESIGN.md`, `qeToy/DUMP_FORMAT.md`, the artifact https://claude.ai/artifact/15cqc6BgRRGuGDpNvs7PSv |

What each tag means for testing: `v1.1` is the reference behaviour (thread per k-point, group size up to about the core count). `v2.0` has the same numerics with one host thread and one queue per group and no limit from the core count. The next step (stage B below) is not done.

## 2. Build (as done on the GB10; adapt the compiler path and architecture)

HPC SDK 25.9 here (`nvfortran`, `nvc`), serial HDF5 built with the same Fortran compiler (its `.mod` files are compiler specific):

```
qeToy/tools/build_hdf5.sh 25.9            # installs ~/opt/hdf5-1.14.6-nvhpc25.9 (needs network once for the tarball, gcc)
mkdir build_a100 && cd build_a100
cmake .. -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_C_COMPILER=nvc -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DQE_ENABLE_CUDA=ON -DQE_ENABLE_OPENACC=ON -DQE_ENABLE_OPENMP=ON -DQE_ENABLE_MPI=OFF -DQE_ENABLE_PROFILE_NVTX=ON \
  -DQE_ENABLE_HDF5=ON -DHDF5_ROOT=$HOME/opt/hdf5-1.14.6-nvhpc25.9 -DHDF5_USE_STATIC_LIBRARIES=TRUE \
  -DNVFORTRAN_CUDA_CC=80 -DQE_FFTW_VENDOR=<NVPL on aarch64, FFTW3 or MKL on x86_64>
make -j8 cb_davidson qe_davidson qe_hpsi_test qe_read_test qe_kpoint_test dump_to_h5
```
`NVFORTRAN_CUDA_CC=80` is the A100; the GB10 build did not set it (compute capability 12.1 picked by default). The HDF5 library needs zlib at link time (handled in the root CMake). Without `-DQE_ENABLE_HDF5=ON` only the stream-file dumps are readable; the datasets are HDF5, so keep it on. Notes that bit me: CUDA Fortran is needed for the cuBLAS/cuSOLVER calls; `ndims = 64` in `FFTXlib/src/fft_param.f90` limits the number of cuFFT plans (one per power of two up to `fft_cols_eff`).

To compare v1.1 and v2.0 on the same machine, check out each tag in its own worktree and build each in its own directory:
`git worktree add ../KS_v11 v1.1` and `git worktree add ../KS_v20 v2.0` (the later commits on the branch, 366f7a7 and 882e3dc, add the retry loop, the LSDA comparison and tools; use the branch head for the current qeToy).

## 3. Correctness gates (do these first)

Both harnesses compare each case with a saved baseline byte for byte and run every case twice. Byte-identical results held on the GB10 for v2.0 against v1.1. **On the A100 the cuBLAS and cuSOLVER kernels differ, so baselines saved on the GB10 will not match bit for bit: save the baseline again on the A100 with the v1.1 build, then compare later builds with it.** (`build_marco/qeToy_datasets/baselines_gb10/` has the GB10 baselines, for reading only.)

```
D=<path>/qeToy_datasets
# cbToy: 9 cases (several batch sizes, with and without the overlap mode)
examples/regress_cbtoy.sh <v1.1 build>/bin /tmp/rb save
examples/regress_cbtoy.sh <v2.0 build>/bin /tmp/rb check            # expect: REGRESSION cbToy: PASSED
# qeToy: Al norm-conserving at 1, 4, 16 slots and ultrasoft at 1, 4, 8, each also against pw.x
qeToy/tools/regress.sh <v1.1 build>/bin /tmp/rq $D/al_nc_24/dump $D/al_us_8/dump $D/al_nc_24/ref_cpu.xml $D/al_us_8/ref_cpu.xml save
REPEAT=2 qeToy/tools/regress.sh <v2.0 build>/bin /tmp/rq $D/al_nc_24/dump $D/al_us_8/dump $D/al_nc_24/ref_cpu.xml $D/al_us_8/ref_cpu.xml check
```
Expected on the GB10: differences to pw.x 7.5e-6 eV (Al nc), 5.1e-6 eV (Al ultrasoft); results independent of the number of slots and identical from run to run. If a baseline comparison fails on the A100 but the pw.x differences are still a few 1e-6 eV and two runs agree, the cause is a different summation order, not a bug. A real failure shows as a difference to pw.x above about 1e-4 eV, a change between two identical runs, or different results for different slot counts.

Unit-type tests (stream or HDF5 dumps, GPU): `qe_read_test.x <dump>` (reader, QE's S|psi> from the dumped data), `qe_hpsi_test.x <dump> [k0 [nb]]` (H|psi> and S|psi> on the device against QE's own, 2e-16 to 5e-16 expected; needs a dump with test vectors: only `al_us_8` in the package has them), `qe_kpoint_test.x <dump>` (k-point data computed from the setup against QE's: same G vectors and order, 2e-16 to 4e-15).

## 4. Comparison with pw.x, without QE (dataset runner)

```
qeToy/tools/run_datasets.sh <build>/bin $D "16 32 64 128" "2 4"
```
Six cases (details in each `case.conf`): `al_nc_24` (413 k-points, 12 bands), `al_us_8` (ultrasoft, 29 k-points), `AlNC40` (413 k-points, 16 bands, 40 Ry), `CuUS` and `CuPAW` (145 k-points, 20 bands), `FeLSDA` (29 k-points x 2 spins). It prints the recorded pw.x `cegterg` times (CPU and GPU build, recorded on the GB10 machine, so they are **not** A100 numbers), the qeToy Davidson time, the steady-state ms per k-point (median over groups without the first one) and the maximum eigenvalue difference to the CPU and GPU pw.x results stored in the package. To get pw.x times on the A100, rerun `qeToy/tools/compare_pw.sh` there (needs QE, a GPU build of pw.x and the dumper; see its header) or run pw.x on the `scf.in`/`nscf.in` stored in each case.

GB10 reference numbers (v2.0, 32 slots, Davidson seconds; pw.x `cegterg` seconds, CPU, GPU):

| case | pw.x CPU | pw.x GPU | qeToy david=2 | qeToy david=4 | max diff to pw.x |
|---|---|---|---|---|---|
| AlNC40 | 21.9 | 25.4 | 3.47 | 2.67 | 5.8e-6 eV |
| CuUS | 8.7 | 14.2 | 2.07 | 1.85 | 2.2e-6 eV |
| CuPAW | 6.8 | 11.1 | 2.08 | 1.17 | 1.2e-6 eV |
| FeLSDA | 5.0 | 8.0 | 1.35 | 0.99 | 3.3e-6 eV, one top band 6e-3 eV |

Known exception: FeLSDA k-point 12, highest spin-down band (52 eV). CPU pw.x and qeToy david=4 give 52.0356 eV, GPU pw.x and qeToy david=2 give 52.0420 eV. Probably a competing state just outside the computed set (not verified; try one more band). Every other of the 928 values agrees to better than 1e-5 eV. The first group of every qeToy run includes about 0.5 s of CUDA module loading on the GB10, a large share of the short runs.

## 5. Performance measurements worth taking on the A100

The goal there is to judge stage B and the synchronisation points, which were deliberately not tuned on the GB10.
1. Steady-state ms per k-point against the number of slots (16, 32, 64, 128) for `AlNC40` and `CuUS`, `david=2` and `4`. GB10 values, Al 12 bands, 413 k-points: 3.8 ms (v1.1, 16 threads), 3.3 to 3.4 (v2.0, 16 slots), 2.59 (32), 2.33 (64).
2. Launch counts and GPU busy fraction:
   ```
   nsys profile -t cuda,nvtx -o prof --force-overwrite true <build>/bin/qe_davidson.x < input.txt
   python3 qeToy/tools/nsys_summary.py prof.nsys-rep
   ```
   GB10 (64 k-points of `al_nc_24`, 16 slots): 16325 launches for v2.0 (16075 for v1.1), GPU busy about 59% in steady state, mean kernel 8.6 to 12 microseconds. The pw.x GPU numbers show it is dominated by launch and synchronisation costs for these small problems.
3. Limits to probe: the batched cuSOLVER eigensolver (`laxlib_cdiaghg_gpu_batched`, requires leading dimension equal to the matrix size) for `david*nbnd` above 64 and up to 128; memory of the FFT work array (`fft_cols_eff`, a 1 GiB budget is built in) and of `psi`, `hpsi`, `spsi` for many slots; slots above 64.
4. The input namelist `&qe` of `qe_davidson.x` (stdin): `dumpdir`, `nk_batches` (slots per group), `david`, `ethr`, `fft_cols_max`, `nks_max`, `nonlocal_batched`, `kmode` (`dump`, `list`, `path`, `mp`) with `nk1..nk3`, `k1..k3`. The k-points of `list` and `path` follow the namelist.

## 6. Debugging playbook

Tools: `compute-sanitizer --tool memcheck|racecheck|synccheck <build>/bin/qe_davidson.x < input` on a small case (`nks_max=32`, 4 slots), `nsys` as above. Run any change twice and compare the output files byte for byte; compare 1, 4 and 16 slots. Switches that isolate a part: `nonlocal_batched=.false.` (per-slot nonlocal GEMMs and per-thread S|psi> instead of the batched ones), `start_from_dump=.true.` (start from QE's converged vectors: if the Davidson does not stay on them, the operator is wrong), `kmode='dump'` against a computed k-point set (isolates the k-point setup), `qe_hpsi_test.x` for the operators alone.

Rules that cost time when broken (all from the batching work):
- Neutral zero padding of every contracted dimension when slots have different sizes: rows beyond `npw` of `vkb` and of the vectors are zero; the padded part of the reduced matrices has `sc = 1` and `hc = i*1e3` on the diagonal so that the padded eigenpairs decouple and sort last.
- Arrays are slot-indexed; the slot is the position among the still-active k-points and `slot2batch` maps it to the batch (the index of `vkb_batched`, `igk_batched`, `g2kin_batched`, ...). A slot count that changes between iterations needs matching arrays.
- `cegterg` writes `npw_comp`, `nb1_comp` and `slot2batch` only on the device: the batched nonlocal routine refreshes them on the host with `acc update host` (`qe_h_psi_batched_iter`), because it needs them as GEMM sizes. Forgetting this gave garbage sizes, illegal addresses and ghost eigenvalues at the start of qeToy.
- The batched cuSOLVER call needs leading dimension equal to the matrix size. The cuFFT batched path needs `nr1x = nr1` (checked at start-up).
- The `S|psi>` callback relies on the batched H|psi> having been called just before on the same columns (`sp_valid`); otherwise it falls back to the per-batch path.
- Library calls must run on the stream of the queue in use: with one host thread and one queue (v2.0) the handles of index 1 are used (`clock_thread = 1`).
- The table of QE's radial beta functions ends just beyond the largest |k+G| of the run and `interp_beta` returns zero for the last three nodes of any table: the dumper extends the table first (`init_tab_beta`).
- Run-to-run reproducibility and agreement between different slot counts are the cheapest correctness test; a race shows up there first.

## 7. Open work

Stage B (not started): replace the per-slot kernel loops in `KS_Solvers/Davidson/cegterg.f90` by kernels over all slots or batched calls. The list from the Stage A report (line numbers of `49db2b0`): copy of the start vectors (222); first `hc` symmetrisation (237); the `lrot` branch (252, 263); reordering of `vc` through `vc_tmp` (294); zeroing of `vc` rows and the pointer arrays (314, 369); scaling `psi` by `-ew` (403); zeroing the garbage rows of `psi` (419); copies of `psi` and `hpsi` into the shared slot arrays and back (429, 457); the per-slot `s_psi_ptr` call (478); second GEMM pointer updates and the `hc`/`sc` symmetrisation (558, 568, 579, 595); filling the padded diagonalisation arrays (620); retrieval of `vc` and `ew` (662); convergence bookkeeping and the `e` copy (716, 730); the final `evc` GEMM per slot and the restart refresh (747, 768). The shared copies `psi_fft_comp`/`hpsi_fft_comp` still exist and could go once the H|psi> callbacks accept slot-indexed `psi` with masked slots.

Synchronisation: the points left are the `acc wait` before host reads of device data (convergence flags, `notcnv`) and at the end of the routine. On the GB10 about 27 000 `cuStreamSynchronize` and 8 000 `cudaDeviceSynchronize` calls were seen per 128 k-points in v1.1. A cheap, safe first step on the A100: pack the small host-to-device transfers of `qe_nonlocal_hpsi_batched` (the `copyin(c0, nb1, m, bt)` clauses and seven pointer-array updates per call) into one transfer; a device-side dependency `!$acc wait(q1) async(q2)` can replace a host `wait` where two queues exchange data (not needed with one queue).

Other open points: QE-side dumper not built through QE's own Makefile or CMake and tested only against libraries of July (old `interp_beta`/`init_tab_beta` calls without `omega`; the committed source uses the current API); only one-atom cells tested (no several species, no non-collinear spin, no gamma-only); FFT without stick pruning (deliberately dropped: the batched 3-D cuFFT keeps the launch count low); no symmetry (every k-point is computed directly, by decision).

## 8. Decisions taken by the owner (do not revisit without asking)
- No symmetry anywhere: every k-point of the set is computed directly.
- Synchronisation removal and anything about the thread structure is judged on the A100, not on the GB10.
- The GB10 is small: keep groups and runs modest there.
