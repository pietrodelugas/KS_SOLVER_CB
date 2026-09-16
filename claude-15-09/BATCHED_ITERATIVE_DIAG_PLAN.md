# Plan: batch the iterative-loop diagonalization in `cegterg`

Branch: `first_diag_single_batched_call` (on top of tag `nishas_run`).
Status of this file: planning only, written at the end of the session, no
further code changes made after this. Resume from here.

## Where things stand

The first reduced-Hamiltonian diagonalization in `cegterg` (the one before
the `iterate: DO kter = 1, maxter` loop, where `nbase == nvec` is identical
for every k-point in the batch) is already batched into a single
`laxlib_cdiaghg_gpu_batched` call, executed by one thread, and verified
bit-for-bit identical to an unmodified `nishas_run` build across all
shipped `cbToy` examples plus an 11 k-point / `nk_batches=4` case (commit
`86a76c6`). Key lesson learned there, which applies to everything below:

> **The batched cuSOLVER eigensolver (`cusolverDnZheevjBatched`) requires
> `ldh == n` for every slot.** If the shared array's declared leading
> dimension differs from the size actually passed as `n`, only slot 1 comes
> back correct; every other slot is silently corrupted (not NaN — either
> stale zeros or garbage, depending on when the bug was hit). This was the
> `nbndx`(576) vs `nbnd`(64) bug. Any array we build for a batched call
> **must** be allocated with leading dimension exactly equal to whatever
> `n` we pass that call, no bigger.

The second diagonalization, inside `iterate: DO kter = 1, maxter`, is
**untouched** — still solved independently per thread, each holding its own
lock (`omp_set_lock(cegterg_locker)`). This plan is about batching that one
too.

## The new problem this call has that the first one didn't

`nbase` grows independently per k-point/thread every iteration
(`nbase = nbase + notcnv`, and `notcnv` — how many roots are still
unconverged — evolves independently per k-point). Worse: on `notcnv == 0`
a thread currently does `EXIT iterate` immediately, i.e. **different
k-points in the same batch leave the loop on different iterations.** A
team-wide `!$omp barrier`/`!$omp single` inside the loop cannot tolerate
that — a thread that has exited never reaches a later iteration's barrier,
and everyone else hangs.

One thing *is* already synchronized for free: `dav_iter == maxter` fires on
the exact same `kter` for every thread, since `maxter` is the same
compile-time constant everywhere. Only the individual early-convergence
exit (`notcnv == 0`) needs to change.

## Design decided (confirmed with the user)

- **Lockstep, not team-shrinking.** A k-point that converges does not
  leave the OpenMP team. It sets a per-thread `my_done` flag, freezes its
  own `dav_iter`, and skips all its own per-iteration compute (`h_psi`,
  `s_psi`, basis expansion, overlap build) from then on — but it still
  reaches every barrier/single every remaining iteration, so the team
  never has a missing member.
- **Converged members are excluded from the batched call itself**
  ("their matrix is not sent as input to the batched diagonalizer") —
  not merely padded and included. This needs a *compacted* slot index:
  active threads get a contiguous `1..n_active` slot number (computed as
  `1 + count of active threads with a smaller i_batch`), and the batched
  call is made with `n_k_arg = n_active`, not the full `n_k`.
- **`nbase_max` sizing, not fixed `nvecx`.** Discussed and explicitly
  chosen over "always pad to `nvecx`": the batched Cholesky/triangular
  solve/Jacobi-eigensolver steps all pay for the *full* matrix size handed
  to them regardless of how much of it is real data, and `nvecx` (=
  `david * nbnd`) can be many times larger than the actual `nbase` early
  in the iteration (9x in the `si2_11_points.in` test case). Accepted
  trade-off: this requires reallocating the shared work arrays to
  `(nbase_max, nbase_max, n_active)` whenever `nbase_max` changes
  (possibly every iteration), which `nvecx`-fixed sizing would have
  avoided. Chosen anyway — "the philosophy is to go in batches", i.e.
  actually save the compute, not just avoid corruption.
- **Collective loop exit.** `iterate: DO kter = 1, maxter` now exits when
  `ALL(done_comp(1:n_k)) .OR. kter == maxter`, computed once by the
  executing thread and read by everyone, so the team leaves the loop
  together.

## Staged implementation plan

Each stage should build, run, and be checked against the unmodified
`nishas_run` baseline (same procedure as before: `git worktree add` a
clean `nishas_run` checkout, configure with the same NVPL BLAS/LAPACK/FFTW
+ `nvc` C compiler flags recorded below, run the same example inputs,
diff the printed bands and `dav_iter` counts) before moving to the next
stage. Use inputs where k-points in the same batch visibly converge at
different iteration counts — `examples/si2_11_points.in` with
`nk_batches=4` (11 k-points, partial last batch of 3, already used before)
is a reasonable starting point; also worth constructing an input where
`ecutwfc`/`ncell` differ enough between k-points in one batch that
`dav_iter` clearly differs across the batch (check the per-batch
`dav_iter, nhpsi, notcnv` print at the end of each round in
`cb_davidson_main.f90` — that is the ground truth to diff against).

### Stage 1 — lockstep control flow only, no batching yet

Goal: get the loop-control rewrite right in isolation, before touching the
diagonalization call at all. This is the highest-risk part (touches every
thread's control flow, not just the diagonalization), so validate it on
its own first, with the second diagonalization *left exactly as it is
today* (per-thread, locked, unbatched).

Changes, all in `KS_Solvers/Davidson/cegterg.f90`:

1. Add `LOGICAL :: my_done` (init `.FALSE.` before `iterate: DO`).
2. Guard `dav_iter = kter` (currently unconditional, top of loop) with
   `IF (.NOT. my_done) dav_iter = kter`.
3. Wrap the existing per-iteration body (basis expansion / `h_psi` /
   `s_psi` / overlap build — everything between `CALL start_clock('cegterg:update')`
   and the existing `CALL diaghg(...)` call) in `IF (.NOT. my_done) THEN ... END IF`.
   Note `notcnv` naturally goes to 0 for a done thread and most of the
   existing loops/`ZGEMM`s are already bounded by `notcnv`, so this guard
   is mostly a safety/clarity measure, not strictly required for those —
   but it *is* required for the "refresh evc" `ZGEMM` block, which
   currently runs unconditionally inside the `IF (notcnv==0 .OR. ...)`
   check every time that check is true.
4. In that `IF (notcnv == 0 .OR. nbase+notcnv > nvecx .OR. dav_iter == maxter) THEN`
   block: replace `EXIT iterate` under `IF (notcnv == 0)` with
   `my_done = .TRUE.` (no exit). Leave the `dav_iter == maxter` branch's
   `EXIT iterate` as-is (already synchronized, see above) — actually,
   with lockstep this branch should probably also just fall through to
   the new collective-exit check at the bottom rather than exiting
   directly, for symmetry; decide during implementation whether keeping
   both exit points causes any issue (it shouldn't, since all threads hit
   `dav_iter==maxter` together) but simplifying to one exit point (the
   new collective check) is cleaner.
5. At the very bottom of the loop body (after the still-untouched,
   still-per-thread-locked second diagonalization call), add the
   collective exit check. Since this stage does not yet have the shared
   `done_comp` array (that is Stage 2's territory), a simple first cut:
   reuse `mp_bcast`/`mp_sum` over... no — this needs an OpenMP-visible
   shared array, not an MPI reduction (MPI reduces across *processes*,
   not across the OpenMP threads within one process that are each
   handling a different k-point). So Stage 1 actually needs *some*
   shared array to test "is everyone done" even before real batching
   exists. Simplest: add a small shared `LOGICAL, INTENT(INOUT) :: done_comp(n_k)`
   dummy argument now (allocated once in `cb_davidson_main.f90`, fixed
   size `nk_batches`, just like `hc_comp` etc. were for the first call),
   write `done_comp(i_batch) = my_done` every iteration (cheap, no
   `!$acc` involved, plain host array), `!$omp barrier`, then
   `IF (ALL(done_comp(1:n_k)) .OR. kter == maxter) EXIT iterate` — every
   thread computes this identically and independently, no `!$omp single`
   needed for the check itself (each thread does the same read-only
   reduction over a small host array).
6. Thread through `cb_davidson_main.f90`: allocate `done_comp` (`LOGICAL`,
   size `nk_batches`), reset it to `.FALSE.` at the top of each `cegterg`
   call's relevant scope (needs to be reset per-outer-`ik`-iteration,
   *not* per Davidson `kter` iteration — a thread must not enter a new
   k-point's `cegterg` call with a stale `.TRUE.` from a previous batch;
   since `done_comp` is indexed by `i_batch` and shared across the whole
   `nk_batches`-sized array, and it's declared fresh, this is naturally
   handled as long as it's set `.FALSE.` at loop entry inside `cegterg`,
   not in the driver — simplest: initialize `done_comp(i_batch) = .FALSE.`
   as the very first thing each `cegterg` call does, before the iterate
   loop, since each thread only ever writes its own index).

Validation for Stage 1: results must be **identical** to the current
committed state (86a76c6) — this stage changes *when* a thread exits the
loop internally (later, via idle spinning) but must not change the
*numerical outcome*, since the diagonalization itself is still done
exactly as before, just possibly on later iterations for k-points that
would have exited early (they now idle instead, doing nothing different
numerically). `dav_iter` for a "done early" thread must still report its
*true* convergence iteration (frozen), not the batch's overall iteration
count — this is the main thing to check in the "batch N, dav_iter,
nhpsi, notcnv" print at the end of each round; it should be unchanged
from the current committed baseline.

### Stage 2 — batch the call, compacted, fixed-size (`nvecx`) arrays

Goal: get the compaction/batching mechanics right, deferring the
allocation-lifecycle complexity of `nbase_max` by using a fixed
`(nvecx, nvecx, n_k)` array for now (same shape discipline as `hc`/`sc`
themselves, so no reallocation needed — but note this reintroduces the
`ldh == n` risk if not careful, see below).

Changes:

1. Add shared work arrays sized `(nvecx, nvecx, n_k)` (yes, `nvecx`, not
   `nvec` — deliberately different shape from the first call's arrays)
   plus `INTEGER, INTENT(INOUT) :: nbase_comp(n_k)`, allocated once in
   `cb_davidson_main.f90`, `!$acc enter data create`'d once up front
   (same as the first-diag arrays — recall the earlier bug where
   per-call dynamic `enter/exit data` from multiple concurrent threads
   corrupted all but one slot; do **not** repeat that mistake here).
2. Each active thread reports `nbase_comp(i_batch) = nbase` (or leaves
   it alone / marks it excluded some other way for done threads — a
   `done_comp` array already exists from Stage 1, so "active" is simply
   `.NOT. done_comp(k)`).
3. `!$omp barrier`.
4. Every thread independently computes:
   - `n_active = COUNT(.NOT. done_comp(1:n_k))`
   - `nbase_max = MAXVAL(nbase_comp(1:n_k), MASK=.NOT. done_comp(1:n_k))`
     (guard the case `n_active == 0`, though that should already have
     triggered the collective exit before reaching here)
   - `my_slot = 1 + COUNT(.NOT. done_comp(1:i_batch-1))` if
     `.NOT. done_comp(i_batch)`, otherwise irrelevant/unused for a done
     thread.
5. Active threads stage into `hc_comp(1:nbase,1:nbase,my_slot)` /
   `sc_comp(...)`, then zero-pad the rest up to `nbase_max` with the
   decoupled block-diagonal scheme from `demo_diaghg_gpu_batched.f90`
   (large diagonal on `hc`, `1` on `sc`, zero cross-terms) for
   `nbase+1 .. nbase_max` — **but only up to `nbase_max`, not `nvecx`**,
   even though the array's declared shape is `(nvecx,nvecx,n_k)`. This
   means the actual `diaghg` call passes `n = nbase_max`, `ldh = nvecx`.
   **This is exactly the mismatch that corrupted the first call.** So
   Stage 2, as scoped ("fixed `nvecx`-shaped array, no realloc"), is
   only safe if we additionally pass `n = nvecx` (i.e. pad *all the way*
   to `nvecx`, not just to `nbase_max`) — which defeats the purpose of
   `nbase_max` and is exactly the "always use `nvecx`" alternative that
   was explicitly rejected. **Conclusion reached while writing this
   plan: Stage 2 as "fixed-size array + `nbase_max`" is not actually a
   safe intermediate step — the `ldh==n` constraint means `nbase_max`
   sizing and dynamic array (re)allocation are not separable.** Revise:
   Stage 2 should instead validate *compaction only*, padding all the
   way to `nvecx` (fixed size, `ldh = n = nvecx`, safe per the known
   constraint), and Stage 3 replaces the fixed `nvecx` pad target with a
   reallocated-each-time `nbase_max` target. This isolates "does
   compaction/indexing work" (Stage 2) from "does dynamic reallocation
   work" (Stage 3), which was the actual point of staging.
6. Single thread (`!$omp single`) calls
   `diaghg(nvecx, nvec, hc_comp, sc_comp, nvecx, ew_comp, vc_comp, n_active, ...)`.
7. Active threads retrieve from their `my_slot`.
8. `!$omp barrier`.

Validation for Stage 2: identical numerical results to Stage 1 (and thus
to 86a76c6) on every example; additionally construct/verify a case where
`n_active < n_k` actually occurs mid-run (i.e. confirm via a temporary
print that compaction is really being exercised, not just falling back
to `n_active == n_k` every time because nothing happens to converge
early in the test inputs used).

### Stage 3 — replace fixed `nvecx` padding with real `nbase_max`

Goal: the actual performance win. Only attempt this once Stage 2's
compaction is proven correct.

1. Change `hc_comp`/`sc_comp`/`vc_comp`/`ew_comp` from fixed
   `(nvecx,nvecx,n_k)` dummy arguments to `ALLOCATABLE` dummy arguments
   (`COMPLEX(DP), INTENT(INOUT), ALLOCATABLE :: hc_comp(:,:,:)`, etc.),
   passed from `cb_davidson_main.f90` as `ALLOCATABLE` actuals (initially
   unallocated, or allocated at some placeholder small size).
2. In the single-thread block, before staging: if `nbase_max` differs
   from the array's current allocated size (`SIZE(hc_comp,1)` or track a
   separate saved-size variable to avoid the `SIZE` call racing with
   in-flight device work), do, in order: `!$acc exit data delete(hc_comp,...)`
   (only if already allocated), `DEALLOCATE`, `ALLOCATE(hc_comp(nbase_max,nbase_max,n_active), ...)`,
   `!$acc enter data create(hc_comp,...)`. All of this strictly inside
   the `!$omp single`/`!$omp end single`, with barriers before (so no
   other thread is mid-stage into the old array) and after (so no other
   thread touches the array before the new one exists).
3. Note `n_active` itself can also change iteration to iteration (as
   more k-points converge), so the array's *third* dimension also needs
   to potentially shrink — reallocate on either `nbase_max` or
   `n_active` changing.
4. Pass `n = nbase_max`, `ldh = nbase_max` (now genuinely equal, since
   the array's actual leading dimension *is* `nbase_max`) to `diaghg`.

Validation for Stage 3: identical numerical results to Stage 2 (and thus
to 86a76c6); additionally, add temporary instrumentation (like the
`DEBUG stage`/`DEBUG retrieve` prints used earlier this session) to
confirm `nbase_max` actually varies run-to-run as expected and the
reallocation path is exercised (not just hit once and then never
change).

## Known risks / things to re-check while implementing

- **`cegterg_locker`**: the second diagonalization currently uses
  `omp_set_lock(cegterg_locker)`/`omp_unset_lock` around the per-thread
  call. Once batched, that lock is no longer needed for this call (only
  one thread calls `diaghg` at all, via `!$omp single`) — remove it here,
  but do not touch its other use (the *first* diagonalization never used
  this lock either, by contrast; double check there is no other caller
  depending on this lock's side effects).
- **`nbgrp > 1` (band-group) broadcasts**: the existing code does
  `mp_bcast(vc, ...)`/`mp_bcast(ew, ...)` over `inter_bgrp_comm` after
  the diagonalization. This must still happen per-thread, for whichever
  threads are active, after retrieval — make sure this isn't
  accidentally done for `my_done` threads (they have nothing new to
  broadcast) or skipped for active ones.
- **MPI band-group root (`my_bgrp_id == root_bgrp_id`)**: same pattern
  as the first call — this condition should be uniform across threads
  (it was for the first call); re-confirm that still holds combined with
  the new per-thread `my_done`/compaction logic (i.e. `my_bgrp_id` must
  not itself become thread-divergent — it shouldn't, it's an MPI
  band-group property, orthogonal to per-k-point convergence).
- **`nk_batches` team-size fix from the first call** (OMP team sized to
  `n_k = min(nk_batches, nks-ik+1)`, not fixed `nk_batches`) must remain;
  this new lockstep logic is an *additional* synchronization requirement
  on top of that, not a replacement for it.
- **`i_batch` vs compacted slot**: be careful never to mix the two up —
  `i_batch` is the stable per-thread identity (1..n_k, fixed for the
  whole `cegterg` call), `my_slot` is only meaningful for active threads
  and changes every iteration as other threads finish. Any per-thread
  local array indexed by `i_batch` (none currently planned) vs the
  shared batched array indexed by `my_slot` must not be confused.
- **`hc`/`sc`/`vc`/`ew` are `!$acc declare device_resident`, fixed
  `(nvecx,nvecx)`/`(nvecx)` shape** — staging always reads
  `hc(1:nbase,1:nbase)` (a sub-block of the full local array), same as
  the first call's `hc(1:nvec,1:nvec)` pattern already committed.
- **Testing must include a case that actually exercises different
  `dav_iter` per k-point within one batch** — none of the currently used
  examples were checked for this explicitly; before trusting Stage 1's
  result, confirm (e.g. via the per-batch printed `dav_iter, nhpsi,
  notcnv` line) that at least one example has k-points converging on
  different iteration counts within the same batch, otherwise the
  lockstep path is never actually tested.

## Environment notes (for resuming without re-deriving them)

- Working build: `/home/qespresso/KS_SOLVER_CB/build_gpu` (already
  configured for this repo checkout; `cmake --build . --target
  cb_davidson -j4` rebuilds).
- Clean-baseline comparison build recipe (`nishas_run`, unmodified),
  needed again for Stage-by-stage diffing:
  ```
  git worktree add /tmp/nishas_run_check nishas_run
  mkdir -p /tmp/nishas_run_check/build_gpu && cd /tmp/nishas_run_check/build_gpu
  cmake .. \
    -DCMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_aarch64/25.9/compilers/bin/nvfortran \
    -DCMAKE_C_COMPILER=/opt/nvidia/hpc_sdk/Linux_aarch64/25.9/compilers/bin/nvc \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DQE_ENABLE_CUDA=ON -DQE_ENABLE_OPENACC=ON -DQE_ENABLE_OPENMP=ON -DQE_ENABLE_MPI=OFF \
    -DQE_ENABLE_PROFILE_NVTX=ON -DQE_ENABLE_TEST=ON \
    -DBLA_VENDOR=NVPL -Dnvpl_DIR=/opt/nvidia/hpc_sdk/Linux_aarch64/25.9/math_libs/nvpl/lib/cmake/nvpl \
    -DQE_FFTW_VENDOR=NVPL -Dnvpl_fft_DIR=/opt/nvidia/hpc_sdk/Linux_aarch64/25.9/math_libs/nvpl/lib/cmake/nvpl_fft
  cmake --build . --target cb_davidson -j4
  ```
  (The `nvc` C compiler and NVPL BLAS/LAPACK/FFTW pins are all required —
  without them, configure either fails outright or silently pulls in
  system FFTW3 + GNU OpenMP, which crashes at runtime with
  `libgomp: TODO` due to an OpenMP runtime conflict with nvfortran's own.)
  Remove the worktree when done: `git worktree remove /tmp/nishas_run_check --force`.
- `examples/si2_11_points.in` now uses `david=4` (was `9`) per this
  session's discussion — keep `nvecx = david*nbnd` modest for testing
  the padding-cost tradeoffs.
- Test loop used this session for a given input:
  ```
  cd /home/qespresso/KS_SOLVER_CB/build_gpu/bin
  ./cb_davidson.x < /home/qespresso/KS_SOLVER_CB/examples/<name>.in > /tmp/run.log 2>&1
  grep -c NaN /tmp/run.log
  grep -A3 "bands (ev)" /tmp/run.log
  ```
  plus `diff` against the equivalent baseline run for exact-match checks.

## Current commit state

- Branch `first_diag_single_batched_call`, commit `86a76c6`: first
  diagonalization batched (single call, single thread), verified against
  `nishas_run`. This plan starts from there; no code for the iterative
  call has been written yet.
