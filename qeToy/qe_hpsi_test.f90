!
! qe_hpsi_test: M2/M3 check. H|psi> (kinetic + local + non local) for a group of k-points at once (one slot per k-point, all bands),
! on the GPU, against hpsi of the dump.   usage: qe_hpsi_test.x dumpdir [k0 [nb]]
!
PROGRAM qe_hpsi_test
  USE qe_ham_module
  USE qe_batch_module
  USE fft_scalar_cuFFT
  USE mytime, ONLY : clock_thread, clock_cuda_stream
  USE openacc, ONLY : acc_get_cuda_stream
  USE laxlib_cusolver_handles, ONLY : initialize_cublas_handles
  IMPLICIT NONE
  CHARACTER(LEN=256) :: dumpdir, arg
  INTEGER :: k0, nb, ib, ik, npw, i, cols, i_b, ind_p, nbndx
  INTEGER, ALLOCATABLE :: npw_comp(:), nbnd_comp(:), nb1_comp(:), slot2batch(:)
  COMPLEX(dp), ALLOCATABLE :: psi(:,:,:), hpsi(:,:,:), hp(:,:), sp(:,:), hn(:,:), ev(:,:)
  INTEGER, ALLOCATABLE :: igk(:)
  REAL(dp), ALLOCATABLE :: g(:), h(:), s(:), et(:)
  COMPLEX(dp), ALLOCATABLE :: vk(:,:)
  REAL(dp) :: err, emax
  IF (COMMAND_ARGUMENT_COUNT() < 1) STOP 'usage: qe_hpsi_test.x dumpdir [k0 [nb]]'
  CALL GET_COMMAND_ARGUMENT(1, dumpdir)
  k0 = 0; nb = 4
  IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN; CALL GET_COMMAND_ARGUMENT(2, arg); READ(arg,*) k0; END IF
  IF (COMMAND_ARGUMENT_COUNT() >= 3) THEN; CALL GET_COMMAND_ARGUMENT(3, arg); READ(arg,*) nb; END IF
  CALL qe_read_global(TRIM(dumpdir))
  IF (.NOT. dump_has_tests) STOP 'the test needs a dump with test vectors'
  IF (k0 + nb > nkstot) STOP 'group beyond the number of k-points'
  CALL qe_setup_grid()
  CALL init_clocks(.TRUE.)
  CALL initialize_cublas_handles(1)
  clock_thread = 1; clock_cuda_stream = acc_get_cuda_stream(clock_thread)
  cols = 1
  DO WHILE (2*cols <= nb*nbnd); cols = 2*cols; END DO
  dim_fft = nnr
  i_b = 1
  DO WHILE (i_b <= cols)
     CALL init_plan_cfft3d_gpu(nr1, nr2, nr3, nr1x, nr2x, nr3x, i_b, ind_p)
     i_b = 2*i_b
  END DO
  CALL qe_alloc_batch(nb, cols, nbnd)
  CALL qe_open_k(TRIM(dumpdir))
  CALL qe_load_group(k0, nb)
  nbndx = nbnd
  ALLOCATE(psi(npwx,nbndx,nb), hpsi(npwx,nbndx,nb), npw_comp(nb), nbnd_comp(nb), nb1_comp(nb), slot2batch(nb))
  ALLOCATE(igk(npwx), g(npwx), h(npwx), s(npwx), et(nbnd), vk(npwx,MAX(nkb,1)))
  ALLOCATE(hp(npwx,nbnd), sp(npwx,nbnd), hn(npwx,nbnd), ev(npwx,nbnd))
  psi = 0; hpsi = 0
  DO ib = 1, nb
     ik = k0 + ib; npw = ngk(ik)
     CALL qe_read_k(ik, igk, g, h, s, vk, et, ev, hp, sp, hn)
     psi(1:npw,:,ib) = ev(1:npw,:)
     npw_comp(ib) = npw; nbnd_comp(ib) = nbnd; nb1_comp(ib) = 1; slot2batch(ib) = ib
  END DO
  !$acc enter data copyin(psi, npw_comp, nb1_comp, slot2batch) create(hpsi)
  CALL qe_h_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_comp, nb1_comp, psi, hpsi, 1, nb, slot2batch)
  !$acc update host(hpsi)
  emax = 0
  DO ib = 1, nb
     ik = k0 + ib; npw = ngk(ik)
     CALL qe_read_k(ik, igk, g, h, s, vk, et, ev, hp, sp, hn)
     err = MAXVAL(ABS(hpsi(1:npw,:,ib) - hp(1:npw,:)))
     WRITE(*,'("k-point ",i5," npw ",i5,"  max|H psi - QE| = ",es10.2)') ik, npw, err
     emax = MAX(emax, err)
  END DO
  IF (emax > 1.d-12) STOP 'FAILED'
  WRITE(*,'("M3 H psi check (kinetic + local + nonlocal): PASSED")')
END PROGRAM qe_hpsi_test
