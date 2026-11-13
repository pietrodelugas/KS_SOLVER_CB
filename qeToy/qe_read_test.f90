!
! qe_read_test: M1 check of the readers. Reads a dump with qe_ham_module and, for all k-points, verifies from the dumped data
! alone: h_diag, s_diag (as usnldiag), and the non local part and S|psi> against QE's hnl and spsi.
! usage: qe_read_test.x dumpdir
!
PROGRAM qe_read_test
  USE qe_ham_module
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  CHARACTER(LEN=256) :: dumpdir
  INTEGER :: ik, npw, na, nt, ih, jh, ib, ig, spin, o
  INTEGER, ALLOCATABLE :: igk(:)
  REAL(dp), ALLOCATABLE :: g2kin(:), h_diag(:), s_diag(:), et(:), hd(:), sd(:)
  COMPLEX(dp), ALLOCATABLE :: vkb(:,:), evc(:,:), hpsi(:,:), spsi(:,:), hnl(:,:), becp(:,:), ps(:,:), qs(:,:), hn(:,:), sc(:,:)
  REAL(dp) :: ar, e_hd, e_sd, e_nl, e_s, m_hd, m_sd, m_nl, m_s
  IF (COMMAND_ARGUMENT_COUNT() < 1) STOP 'usage: qe_read_test.x dumpdir'
  CALL GET_COMMAND_ARGUMENT(1, dumpdir)
  CALL qe_read_global(TRIM(dumpdir))
  WRITE(*,'("nk=",i5," nbnd=",i3," npwx=",i5," nkb=",i3," grid ",3i4," nnr=",i7," ngms=",i6," okvan=",l1," tests in dump=",l1)') &
       nkstot, nbnd, npwx, nkb, nr1, nr2, nr3, nnr, ngms, okvan, dump_has_tests
  IF (MINVAL(nl) < 1 .OR. MAXVAL(nl) > nnr) STOP 'nl out of range'
  ALLOCATE(igk(npwx), g2kin(npwx), h_diag(npwx), s_diag(npwx), et(nbnd), vkb(npwx,MAX(nkb,1)), hd(npwx), sd(npwx))
  ALLOCATE(evc(npwx,nbnd), hpsi(npwx,nbnd), spsi(npwx,nbnd), hnl(npwx,nbnd), becp(MAX(nkb,1),nbnd), ps(MAX(nkb,1),nbnd), &
       qs(MAX(nkb,1),nbnd), hn(npwx,nbnd), sc(npwx,nbnd))
  CALL qe_open_k(TRIM(dumpdir))
  m_hd = 0; m_sd = 0; m_nl = 0; m_s = 0
  DO ik = nkstot, 1, -1           ! backwards: exercises the random access
     npw = ngk(ik); spin = isk(ik)
     IF (dump_has_tests) THEN
        CALL qe_read_k(ik, igk, g2kin, h_diag, s_diag, vkb, et, evc, hpsi, spsi, hnl)
     ELSE
        CALL qe_read_k(ik, igk, g2kin, h_diag, s_diag, vkb, et)
     END IF
     hd = 0; sd = 0
     hd(1:npw) = g2kin(1:npw) + v_of_0; sd(1:npw) = 1.d0
     DO na = 1, nat; nt = ityp(na); o = ofsbeta(na)
        DO ig = 1, npw; DO ih = 1, nh(nt); DO jh = 1, nh(nt)
           ar = DBLE(vkb(ig,o+ih) * CONJG(vkb(ig,o+jh)))
           hd(ig) = hd(ig) + deeq(ih,jh,na,spin) * ar
           sd(ig) = sd(ig) + qq_at(ih,jh,na) * ar
        END DO; END DO; END DO
     END DO
     e_hd = MAXVAL(ABS(hd(1:npw) - h_diag(1:npw))); e_sd = MAXVAL(ABS(sd(1:npw) - s_diag(1:npw)))
     m_hd = MAX(m_hd, e_hd); m_sd = MAX(m_sd, e_sd)
     IF (dump_has_tests) THEN
        DO ib = 1, nbnd; DO o = 1, nkb
           becp(o,ib) = SUM(CONJG(vkb(1:npw,o)) * evc(1:npw,ib))
        END DO; END DO
        ps = 0; qs = 0
        DO na = 1, nat; nt = ityp(na)
           DO ib = 1, nbnd; DO ih = 1, nh(nt); DO jh = 1, nh(nt)
              ps(ofsbeta(na)+ih,ib) = ps(ofsbeta(na)+ih,ib) + deeq(ih,jh,na,spin) * becp(ofsbeta(na)+jh,ib)
              qs(ofsbeta(na)+ih,ib) = qs(ofsbeta(na)+ih,ib) + qq_at(ih,jh,na) * becp(ofsbeta(na)+jh,ib)
           END DO; END DO; END DO
        END DO
        DO ib = 1, nbnd
           DO ig = 1, npw
              hn(ig,ib) = SUM(vkb(ig,1:nkb) * ps(1:nkb,ib))
              sc(ig,ib) = evc(ig,ib) + SUM(vkb(ig,1:nkb) * qs(1:nkb,ib))
           END DO
        END DO
        e_nl = MAXVAL(ABS(hn(1:npw,:) - hnl(1:npw,:))); e_s = MAXVAL(ABS(sc(1:npw,:) - spsi(1:npw,:)))
        m_nl = MAX(m_nl, e_nl); m_s = MAX(m_s, e_s)
     END IF
  END DO
  CALL qe_close_k()
  WRITE(*,'("max over all k-points: |h_diag|",es10.2," |s_diag|",es10.2," |V_NL psi|",es10.2," |S psi|",es10.2)') m_hd, m_sd, m_nl, m_s
  IF (MAX(m_hd, m_sd, m_nl, m_s) > 1.d-12) STOP 'FAILED'
  WRITE(*,'("M1 reader check: PASSED")')
END PROGRAM qe_read_test
