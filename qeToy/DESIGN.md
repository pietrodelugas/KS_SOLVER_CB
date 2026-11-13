# qeToy: threading and batching design (draft, nothing implemented yet)

`qeToy` is the miniapp that runs the batched Davidson of this repo on the Hamiltonian of a real QE run (data from the QE dumper, see `DUMP_FORMAT.md`).
It keeps the **threading and batching skeleton of `cbToy`** (`cbToy/cb_davidson_main.f90`, `KS_Solvers/Davidson/cegterg.f90`);
what changes is where the data comes from and the H|psi>, S|psi> and preconditioner stages.

## 0. Dependencies and autonomy
- **qeToy needs nothing from QE.** It is built with this repo's CMake only (`FFTXlib`, `UtilXlib`, `LAXlib`, `KS_Solvers`, which are vendored here, plus OpenACC/OpenMP/CUDA)
  and reads plain stream files. No QE source, module file or library is used at compile time or at run time.
- **Only the dumper needs QE**: it lives in the QE tree (`PP/src/dump_ham.f90`, target `dump_ham.x`), not in this repo; format in `DUMP_FORMAT.md`, checker in `tools/read_dump_check.f90`.
  It is the data-generation step, run once per system.
- **To run without QE present, the dump files are shipped with the miniapp.** Production dumps do not need the test vectors (`evc`, `hpsi`, `spsi`, `hnl` are four fifths of the size
  of the Al file: 37 MB for 413 k-points), so the dumper gets an option to omit them (to be added) and keeps them only for a few k-points of the validation dumps.
  Expected sizes: about 7 MB for the 413 k-point Al run, well under 1 MB for the validation dumps; the reference eigenvalues travel with them.
- **Not autonomous, by design:** a new system (another metal, another cutoff) needs QE once to produce its dump. Making the miniapp generate its own data would need
  a UPF reader, the radial Bessel tables, spherical harmonics, the local potential and a self-consistent density (or a non-self-consistent atomic-superposition
  potential with its own XC functional), which is a project of its own. Not planned.

## 1. What is taken over from cbToy unchanged
- One OpenMP thread per batch (`!$omp parallel num_threads(n_k)` for every group of `nk_batches` k-points, `n_k` smaller for the last
  group). Thread `i_batch`: `clock_thread = i_batch`, own OpenACC queue and CUDA stream (`acc_get_cuda_stream`), own cuBLAS and
  cuSOLVER handles (`initialize_cublas_handles`, `initialize_cusolver_handles`, `initialize_laxlib_cuda_stream`).
- The same `cegterg` (batched), called once per thread with the same shared arrays: `evc_batched`, `eig_batched`, `hc_c/sc_c/vc_c/ew_c`,
  `done_comp`, `hc_comp_itr`..., `psi_fft_comp/hpsi_fft_comp`, `npw_comp`, `nb1_comp`, `slot2batch`, `eig_new_comp`, `conv_comp`,
  `ew_old_comp`; sizes from `nbnd` and `nbndx = david*nbnd`.
- The batched stages inside `cegterg` (diagonalisations, ZGEMMs, convergence test, `g_psi`, packed-column FFT H|psi>) and the rules of
  the Batching Map (zero padding, slot indexing, waits, run-to-run reproducibility).
- cuFFT plans created at start-up for 1, 2, 4, ... `fft_cols_eff` columns on the grid of the run; `fft_cols_max` input option.

## 2. What is new (cbToy element -> qeToy element)
| cbToy | qeToy |
|---|---|
| `input`, `ggen`, `set_cb_potential` (model built analytically) | `qe_read_global`: `ham_global.bin` (grid, `nl`, `vrs`, `deeq`, `qq_at`, `ityp`, `nh`, `ofsbeta`, `xk`, `wk`, `isk`, `ngk`) |
| `init_k(ik, i_batch)`: `igk`, `ekin` of one k-point | `qe_read_group`: reads the k-points of a group from `ham_k.bin` by offset (record size follows from `ngk`) into batch arrays |
| `ekin_batched(npwx,nb)`, `igk_batched(npwx,nb)`, `vloc(nnr)` | `g2kin_batched`, `igk_batched`, `vrs(nnr,nspin)` with the spin of the k-point per batch (`spin_batched`) |
| none | `vkb_batched(npwx,nkb,nb)`, `hdiag_batched`, `sdiag_batched` (QE's preconditioner data) |
| `init_random_wfcs` | the same at first (QE starts from atomic + random; eigenvalues must agree, iteration counts need not) |
| `cb_g_psi_batched_iter` (toy preconditioner) | same interface with `x = h_diag - e*s_diag`, `denm = 0.5*(1+x+sqrt(1+(x-1)^2))` read through `slot2batch` |
| toy overlap `D` factors | removed: S is `1 + vkb qq vkb^H` (zero for norm-conserving, `uspp = okvan`) |

All batch arrays are device resident (`acc enter data`), refreshed for every group with `acc update device`, as in cbToy.
Sizes are tiny for Al (`vkb`: 7 KB per k-point norm-conserving, 29 KB ultrasoft); the largest object is `fft_array_batched`.

## 3. H|psi> with the QE ingredients: how it is split
`H psi = g2kin*psi + V_loc psi + vkb D vkb^H psi`.
1. **Kinetic + local potential: batched, packed (slot, band) columns** (the routine `cbToy/cb_operators/my_h_psi_batched_iter.f90`
   adapted): `ekin_batched -> g2kin_batched(ig, slot2batch(ik))`; `dfft%nl(igk) -> nl(igk_batched(ig, b))`; the multiplication by `vloc(ir)`
   becomes `vrs(ir, spin of the column's batch)`. Used for the first H|psi> too (every k-point active, `nb1 = 1`).
2. **Non local part: one GEMM pair per slot**, in the same single-threaded phase (no change of `cegterg`): `becp = vkb^H psi_new`
   (`MYZGEMM_H`, thread stream), `ps = D becp` (OpenACC kernel over atom blocks: `ps(ofs+ih,col) = sum_jh deeq(ih,jh,na,spin) becp(ofs+jh,col)`),
   `hpsi += vkb ps` (`MYZGEMM_H` with `beta = 1`). Cost: 3 launches per active slot per H|psi> call.
3. **`S psi`** (ultrasoft/PAW only): per thread, in the existing `s_psi_ptr` callback: `spsi = psi + vkb (qq becp)`; the same two GEMMs and a kernel on the
   thread's stream. Norm-conserving runs use `uspp = .false.`.
4. **`g_psi`**: batched callback as above.

Alternatives, to try only if the profile asks:
- (a) nonlocal as one `cublasZgemmBatched` over the slots with the neutral-padding rule (existing machinery; pointer arrays);
- (b) nonlocal as a new per-thread callback after the copy-back in `cegterg` (threads then work concurrently, but `cegterg` gains an argument).

## 4. FFT grid
`invfft('Wave', f, dfft, howmany)` needs a `fft_type_descriptor`. Build it with `fft_type_init` from the cell and the same cutoffs, then **check** `nr1,nr2,nr3,nr1x,nr2x,nr3x,nnr` against the dump and stop on any difference. `cfft3d_gpu_batched` requires `nr1x = nr1` etc. (leading
dimensions equal to the grid size): true for Al (15^3, 18^3), to be tested for the other metals. QE's `nl` is used as is: it is the index map of the dump.

## 5. Convergence control and outputs
- `ethr`, `empty_ethr = max(5*ethr, 1e-5)` as QE; `btype(ibnd)`: bands QE treats as occupied. **Open:** dump `wg` or the Fermi energy to set it, or use `btype = 1` for the bands below E_F + a margin.
- Eigenvalues written in the layout of `reference/al_nscf_24x24x24_eigenvalues_eV.txt` (k, weight, bands) so the xml and the QE runs can be compared directly; bands and DOS from the same file later.
- LSDA: the spin of a k-point selects `vrs(:,spin)` and `deeq(:,:,:,spin)`; nothing else changes.

## 6. Threads and group size (the main risk)
cbToy needs `n_k` threads for a group of `n_k` k-points; Al has 413 k-points of 12 bands and about 110 plane waves, so one would like groups of 32 to 128 k-points, more than the 20 cores. First version: the cbToy model with `nk_batches` up to about 20 to 32 (threads mostly wait on barriers and the GPU), then measure.
If throughput is limited by thread count or per-thread launches, the next step is the structural change already identified on the Batching Map: all state in shared
3-D arrays indexed by slot, single kernels, no thread per batch. The batch-indexed layout above is chosen so this change does not alter the data files or the H|psi> stages.

## 7. Build
New directory `qeToy/` with its own `CMakeLists.txt` like `cbToy` (target `qe_davidson`, output `qe_davidson.x`), linking `qe_kssolver_davidson`, `qe_lax`, FFTXlib, UtilXlib,
OpenACC and OpenMP; the infrastructure files of `cb_lib` (`mp_*`, `io_global`, clocks, random numbers) are reused through `cb_lib`, the model-specific ones (`ggen`, `set_cb_potential`, `init_k`) are not. Common code (packed FFT routine) is copied and adapted first, and factored out once both miniapps stabilise.

## 8. Steps and checks
| step | content | check |
|---|---|---|
| M1 | module, global and group readers, grid check | sizes and `nl` against the dump; `tools/read_dump_check` numbers reproduced on the device |
| M2 | batched kinetic + local stage | against `hpsi - hnl` of the dump to round-off (all k-points of a group at once) |
| M3 | nonlocal stage, `S psi`, `g_psi` | `hpsi`, `spsi` against the dump to round-off; Al norm-conserving, then ultrasoft |
| M4 | main program with the thread/batch skeleton and `cegterg` | eigenvalues against the xml (to the threshold), several `nk_batches`, identical from run to run, sanitizer clean |
| M5 | Al 413 k-points: bands, Fermi energy, DOS | against `pw.x` CPU and GPU; timing against both |
| M6 | profile, kernel counts, decide on section 6 | numbers on the same pages as for cbToy |
