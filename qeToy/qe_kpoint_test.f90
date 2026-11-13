!
! qe_kpoint_test: the k dependent data computed from the setup of the dump (qe_kpoint_module) against the dumped ones
! (computed by QE), for the k-points of the run.   usage: qe_kpoint_test.x dumpdir [step]   (every step-th k-point)
!
PROGRAM qe_kpoint_test
  USE qe_ham_module
  USE qe_kpoint_module
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  CHARACTER(LEN=256) :: dumpdir, arg
  INTEGER :: nprint = 0, j, ik, npw, step, nchk, ndiffn, ndiffg
  INTEGER, ALLOCATABLE :: igk(:), igk2(:)
  REAL(dp), ALLOCATABLE :: g2kin(:), h_diag(:), s_diag(:), et(:), g2(:), h2(:), s2(:)
  COMPLEX(dp), ALLOCATABLE :: vkb(:,:), vkb2(:,:)
  REAL(dp) :: eg, eh, es, ev, mg, mh, ms, mv
  IF (COMMAND_ARGUMENT_COUNT() < 1) STOP 'usage: qe_kpoint_test.x dumpdir [step]'
  CALL GET_COMMAND_ARGUMENT(1, dumpdir); step = 1
  IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN; CALL GET_COMMAND_ARGUMENT(2, arg); READ(arg,*) step; END IF
  CALL qe_read_global(TRIM(dumpdir))
  IF (.NOT. has_setup) STOP 'the dump has no setup data (dumper v2)'
  ALLOCATE(igk(npwx), g2kin(npwx), h_diag(npwx), s_diag(npwx), et(nbnd), vkb(npwx,MAX(nkb,1)))
  ALLOCATE(igk2(npwx), g2(npwx), h2(npwx), s2(npwx), vkb2(npwx,MAX(nkb,1)))
  CALL qe_open_k(TRIM(dumpdir))
  mg = 0; mh = 0; ms = 0; mv = 0; nchk = 0; ndiffn = 0; ndiffg = 0
  DO ik = 1, nkstot, step
     CALL qe_read_k(ik, igk, g2kin, h_diag, s_diag, vkb, et)
     CALL qe_make_kpoint(xk(:,ik), isk(ik), npwx, npw, igk2, g2, h2, s2, vkb2)
     nchk = nchk + 1
     IF (npw /= ngk(ik)) THEN; ndiffn = ndiffn + 1; CYCLE; END IF
     IF (ANY(igk2(1:npw) /= igk(1:npw))) ndiffg = ndiffg + 1      ! same G vectors, different order
     IF (ANY(igk2(1:npw) /= igk(1:npw))) CYCLE
     eg = MAXVAL(ABS(g2(1:npw) - g2kin(1:npw))); eh = MAXVAL(ABS(h2(1:npw) - h_diag(1:npw)))
     es = MAXVAL(ABS(s2(1:npw) - s_diag(1:npw))); ev = MAXVAL(ABS(vkb2(1:npw,1:nkb) - vkb(1:npw,1:nkb)))
     IF (ev > 1.d-4 .AND. nprint < 4) THEN
        nprint = nprint + 1
        WRITE(*,'("k-point",i5," xk=",3f9.5," npw",i4)') ik, xk(:,ik), npw
        DO j = 1, nkb
           WRITE(*,'("  k",i4," proj",i2," max diff",es10.2,"  max |vkb|",es10.2," at ig",i5)') ik, j, MAXVAL(ABS(vkb2(1:npw,j)-vkb(1:npw,j))), MAXVAL(ABS(vkb(1:npw,j))), MAXLOC(ABS(vkb2(1:npw,j)-vkb(1:npw,j)),1)
        END DO
     END IF
     mg = MAX(mg, eg); mh = MAX(mh, eh); ms = MAX(ms, es); mv = MAX(mv, ev)
  END DO
  CALL qe_close_k()
  WRITE(*,'("k-points checked: ",i5,"; different npw: ",i4,"; different G order: ",i4)') nchk, ndiffn, ndiffg
  WRITE(*,'("max |g2kin|",es10.2,"  |h_diag|",es10.2,"  |s_diag|",es10.2,"  |vkb|",es10.2)') mg, mh, ms, mv
  IF (ndiffn == 0 .AND. MAX(mg, mh, ms, mv) < 1.d-10) WRITE(*,'("k-point setup check: PASSED")')
END PROGRAM qe_kpoint_test
