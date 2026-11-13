!
! dump_to_h5: converts a dump in stream files (ham_global.bin, ham_k.bin) into ham.h5 (needs a build with QE_ENABLE_HDF5=ON).
! usage: dump_to_h5.x indir outdir [notests]      (notests: leave the test vectors out)
!
PROGRAM dump_to_h5
  USE qe_ham_module
  USE qe_dump_h5
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  CHARACTER(LEN=256) :: indir, outdir, arg
  LOGICAL :: tests
  INTEGER :: ik, npw, ih(nihead)
  INTEGER, ALLOCATABLE :: igk(:)
  REAL(dp), ALLOCATABLE :: g2kin(:), h_diag(:), s_diag(:), et(:)
  COMPLEX(dp), ALLOCATABLE :: vkb(:,:), evc(:,:), hpsi(:,:), spsi(:,:), hnl(:,:)
#if !defined(__HDF5)
  STOP 'dump_to_h5: build with QE_ENABLE_HDF5=ON'
#else
  IF (COMMAND_ARGUMENT_COUNT() < 2) STOP 'usage: dump_to_h5.x indir outdir [notests]'
  CALL GET_COMMAND_ARGUMENT(1, indir); CALL GET_COMMAND_ARGUMENT(2, outdir)
  CALL qe_read_global(TRIM(indir))
  IF (use_h5) STOP 'the input is already HDF5'
  tests = dump_has_tests
  IF (COMMAND_ARGUMENT_COUNT() >= 3) THEN
     CALL GET_COMMAND_ARGUMENT(3, arg); IF (TRIM(arg) == 'notests') tests = .FALSE.
  END IF
  ih = [20251113, nspin, nkstot, nbnd, npwx, nkb, nat, nsp, nhm, nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr, ngms, &
        MERGE(1,0,okvan), MERGE(1,0,doublegrid)]
  CALL EXECUTE_COMMAND_LINE('mkdir -p '//TRIM(outdir))
  CALL dump_h5_create(TRIM(outdir)//'/ham.h5')
  CALL dump_h5_write_global(ih, MERGE(1,0,tests), v_of_0, omega, alat, at, bg, nl, vrs, ityp, nh, ofsbeta, deeq, qq_at, xk, wk, isk, ngk)
  ALLOCATE(igk(npwx), g2kin(npwx), h_diag(npwx), s_diag(npwx), et(nbnd), vkb(npwx,MAX(nkb,1)))
  ALLOCATE(evc(npwx,nbnd), hpsi(npwx,nbnd), spsi(npwx,nbnd), hnl(npwx,nbnd))
  CALL qe_open_k(TRIM(indir))
  DO ik = 1, nkstot
     npw = ngk(ik)
     IF (dump_has_tests) THEN
        CALL qe_read_k(ik, igk, g2kin, h_diag, s_diag, vkb, et, evc, hpsi, spsi, hnl)
     ELSE
        CALL qe_read_k(ik, igk, g2kin, h_diag, s_diag, vkb, et)
     END IF
     IF (tests) THEN
        CALL dump_h5_write_k(ik, npw, nbnd, nkb, igk, g2kin, h_diag, s_diag, vkb, et, evc, hpsi, spsi, hnl)
     ELSE
        CALL dump_h5_write_k(ik, npw, nbnd, nkb, igk, g2kin, h_diag, s_diag, vkb, et)
     END IF
  END DO
  CALL qe_close_k()
  CALL dump_h5_close()
  WRITE(*,'("wrote ",a,"/ham.h5: ",i5," k-points, tests=",l1)') TRIM(outdir), nkstot, tests
#endif
END PROGRAM dump_to_h5
