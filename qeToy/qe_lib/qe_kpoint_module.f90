!
! qe_kpoint_module: the k dependent part of the Hamiltonian for ARBITRARY k-points, from the setup data of the dump
! (no QE code, no symmetry: every k-point is computed directly). For a k-point (cartesian, units 2pi/alat) it gives
!   igk     the G vectors (indices in the list of the dump) with |k+G|^2 tpiba^2 <= ecutwfc, ordered by |k+G|
!   g2kin   tpiba^2 |k+G|^2 (Ry)
!   vkb     beta projectors, as init_us_2 of QE:  vkb = beta(|k+G|) Y_lm(k+G) exp(-i (k+G).tau) (-i)^l
!   h_diag, s_diag   preconditioner data, as usnldiag of QE
! beta(q) is interpolated (cubic Lagrange, 4 nodes) from the table of the dump, as interp_beta of QE.
!
MODULE qe_kpoint_module
  USE qe_ham_module
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  REAL(dp), PARAMETER :: pi = 3.14159265358979323846d0, tpi = 2.d0*pi, fpi = 4.d0*pi
  PUBLIC :: qe_count_pw, qe_make_kpoint
CONTAINS
  !
  INTEGER FUNCTION qe_count_pw(k)
    REAL(dp), INTENT(IN) :: k(3)
    INTEGER :: ig
    REAL(dp) :: gk(3), gcutw
    gcutw = ecutwfc / tpiba**2
    qe_count_pw = 0
    DO ig = 1, ngms
       gk = k + MATMUL(bg, REAL(mill(:,ig), dp))
       IF (SUM(gk*gk) <= gcutw) qe_count_pw = qe_count_pw + 1
    END DO
  END FUNCTION qe_count_pw
  !
  SUBROUTINE qe_make_kpoint(k, spin, npwmax, npw, igk, g2kin, h_diag, s_diag, vkb)
    REAL(dp), INTENT(IN) :: k(3)
    INTEGER, INTENT(IN) :: spin, npwmax
    INTEGER, INTENT(OUT) :: npw
    INTEGER, INTENT(OUT) :: igk(npwmax)
    REAL(dp), INTENT(OUT) :: g2kin(npwmax), h_diag(npwmax), s_diag(npwmax)
    COMPLEX(dp), INTENT(OUT) :: vkb(npwmax,*)
    INTEGER :: ig, n, i, j, na, nt, ih, jh, ofs, nb, lm, iq
    INTEGER, ALLOCATABLE :: list(:), perm(:)
    REAL(dp), ALLOCATABLE :: qq2(:), gk(:,:), qg(:), ylm(:,:), vq(:,:,:), vkb1(:,:)
    COMPLEX(dp), ALLOCATABLE :: sk(:)
    REAL(dp) :: gcutw, g3(3), arg, ar, px, ux, vx, wx, qgr
    INTEGER :: i0
    gcutw = ecutwfc / tpiba**2
    ! G vectors inside the cutoff, ordered by |k+G|^2 (ties: by index, as a stable sort)
    ALLOCATE(list(ngms), qq2(ngms)); n = 0
    DO ig = 1, ngms
       g3 = k + MATMUL(bg, REAL(mill(:,ig), dp))
       IF (SUM(g3*g3) <= gcutw) THEN; n = n + 1; list(n) = ig; qq2(n) = SUM(g3*g3); END IF
    END DO
    IF (n > npwmax) STOP 'qe_make_kpoint: npwmax too small'
    ALLOCATE(perm(n))
    perm = [(i, i=1,n)]
    DO i = 2, n                      ! insertion sort, stable
       j = i - 1; ig = perm(i)
       DO WHILE (j >= 1)
          IF (qq2(perm(j)) > qq2(ig) + 1.d-9) THEN; perm(j+1) = perm(j); j = j - 1; ELSE; EXIT; END IF
       END DO
       perm(j+1) = ig
    END DO
    npw = n
    igk = 1; g2kin = 0; h_diag = 1; s_diag = 1
    IF (nkb > 0) vkb(1:npwmax,1:nkb) = (0.d0,0.d0)
    ALLOCATE(gk(3,n), qg(n))
    DO i = 1, n
       igk(i) = list(perm(i))
       gk(:,i) = k + MATMUL(bg, REAL(mill(:,igk(i)), dp))
       qg(i) = SUM(gk(:,i)**2)
       g2kin(i) = tpiba**2 * qg(i)
    END DO
    IF (nkb > 0 .AND. lmaxkb >= 0) THEN
       ALLOCATE(ylm(n,(lmaxkb+1)**2), vkb1(n,nhm), sk(n), vq(n,nbetam,nsp))
       CALL qe_ylmr2((lmaxkb+1)**2, n, gk, qg, ylm)
       ! radial functions beta(|k+G|)/sqrt(omega) from the table of the dump (a.u.)
       DO nt = 1, nsp
          DO i = 1, n
             qgr = SQRT(qg(i)) * tpiba
             px = qgr/dq_tab - DBLE(INT(qgr/dq_tab)); ux = 1.d0 - px; vx = 2.d0 - px; wx = 3.d0 - px
             i0 = INT(qgr/dq_tab) + 1
             DO nb = 1, nbetam
                IF (i0 + 3 <= nqx) THEN
                   vq(i,nb,nt) = tab(i0,nb,nt)*ux*vx*wx/6.d0 + tab(i0+1,nb,nt)*px*vx*wx/2.d0 &
                               - tab(i0+2,nb,nt)*px*ux*wx/2.d0 + tab(i0+3,nb,nt)*px*ux*vx/6.d0
                ELSE
                   vq(i,nb,nt) = 0.d0
                END IF
             END DO
          END DO
       END DO
       DO na = 1, nat
          nt = ityp(na); ofs = ofsbeta(na)
          DO ih = 1, nh(nt)
             nb = indv(ih,nt); lm = nhtolm(ih,nt)
             vkb1(:,ih) = ylm(:,lm) * vq(:,nb,nt)
          END DO
          DO i = 1, n
             arg = tpi * SUM((k + MATMUL(bg, REAL(mill(:,igk(i)), dp))) * tau(:,na))
             sk(i) = CMPLX(COS(arg), -SIN(arg), dp)
          END DO
          DO ih = 1, nh(nt)
             vkb(1:n,ofs+ih) = vkb1(1:n,ih) * sk(1:n) * (CMPLX(0.d0,-1.d0,dp)**nhtol(ih,nt))
          END DO
       END DO
       ! preconditioner data, as usnldiag
       h_diag(1:n) = g2kin(1:n) + v_of_0; s_diag(1:n) = 1.d0
       DO na = 1, nat
          nt = ityp(na); ofs = ofsbeta(na)
          DO ih = 1, nh(nt); DO jh = 1, nh(nt)
             DO i = 1, n
                ar = DBLE(vkb(i,ofs+ih) * CONJG(vkb(i,ofs+jh)))
                h_diag(i) = h_diag(i) + deeq(ih,jh,na,spin) * ar
                s_diag(i) = s_diag(i) + qq_at(ih,jh,na) * ar
             END DO
          END DO; END DO
       END DO
    ELSE
       h_diag(1:n) = g2kin(1:n) + v_of_0; s_diag(1:n) = 1.d0
    END IF
  END SUBROUTINE qe_make_kpoint
  !
  ! Real spherical harmonics, same numbering and algorithm as ylmr2 of Quantum ESPRESSO (upflib/ylmr2.f90, CPU branch).
  SUBROUTINE qe_ylmr2(lmax2, ng, g, gg, ylm)
    INTEGER, INTENT(IN) :: lmax2, ng
    REAL(dp), INTENT(IN) :: g(3,ng), gg(ng)
    REAL(dp), INTENT(OUT) :: ylm(ng,lmax2)
    REAL(dp), PARAMETER :: eps = 1.0d-9
    REAL(dp) :: cost, sent, phi, c, gmod
    INTEGER :: lmax, ig, l, m, lm, lm1, lm2
    IF (ng < 1 .OR. lmax2 < 1) RETURN
    DO lmax = 0, 20
       IF ((lmax+1)**2 == lmax2) EXIT
    END DO
    IF (lmax == 0) THEN
       ylm(:,1) = SQRT(1.d0/fpi); RETURN
    END IF
    DO ig = 1, ng
       gmod = SQRT(gg(ig))
       IF (gmod < eps) THEN; cost = 0.d0; ELSE; cost = g(3,ig)/gmod; END IF
       sent = SQRT(MAX(0.d0, 1.d0-cost*cost))
       ylm(ig,1) = 1.d0
       ylm(ig,2) = cost
       ylm(ig,4) = -sent/SQRT(2.d0)
       DO l = 2, lmax
          DO m = 0, l-2
             lm = l**2 + 1 + 2*m; lm1 = (l-1)**2 + 1 + 2*m; lm2 = (l-2)**2 + 1 + 2*m
             ylm(ig,lm) = cost*(2*l-1)/SQRT(DBLE(l*l-m*m))*ylm(ig,lm1) - SQRT(DBLE((l-1)*(l-1)-m*m))/SQRT(DBLE(l*l-m*m))*ylm(ig,lm2)
          END DO
          lm = l**2 + 1 + 2*l; lm1 = l**2 + 1 + 2*(l-1); lm2 = (l-1)**2 + 1 + 2*(l-1)
          ylm(ig,lm1) = cost*SQRT(DBLE(2*l-1))*ylm(ig,lm2)
          ylm(ig,lm) = -SQRT(DBLE(2*l-1))/SQRT(DBLE(2*l))*sent*ylm(ig,lm2)
       END DO
       IF (g(1,ig) > eps) THEN
          phi = ATAN(g(2,ig)/g(1,ig))
       ELSE IF (g(1,ig) < -eps) THEN
          phi = ATAN(g(2,ig)/g(1,ig)) + pi
       ELSE
          phi = SIGN(pi/2.d0, g(2,ig))
       END IF
       lm = 1
       ylm(ig,1) = ylm(ig,1)/SQRT(fpi)
       DO l = 1, lmax
          c = SQRT(DBLE(2*l+1)/fpi)
          lm = lm + 1
          ylm(ig,lm) = c*ylm(ig,lm)
          DO m = 1, l
             lm = lm + 2
             ylm(ig,lm-1) = c*SQRT(2.d0)*ylm(ig,lm)*COS(m*phi)
             ylm(ig,lm) = c*SQRT(2.d0)*ylm(ig,lm)*SIN(m*phi)
          END DO
       END DO
    END DO
  END SUBROUTINE qe_ylmr2
END MODULE qe_kpoint_module
