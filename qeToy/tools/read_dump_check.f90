!
! read_dump_check: reads the files written by dump_ham (no QE needed) and rebuilds H|psi> and S|psi> for a few k-points
! from the dumped ingredients only, with a direct (slow) Fourier transform, then compares with QE's own result:
!    H psi = g2kin * psi + V_loc psi + vkb D vkb^H psi        S psi = psi + vkb qq vkb^H psi
! usage: read_dump_check dumpdir [k1 k2 ...]   (default: first, middle and last k-point)
!
PROGRAM read_dump_check
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  INTEGER :: magic, nspin, nkstot, nbnd, npwx, nkb, nat, nsp, nhm, nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr, ngms, okvan, dbl
  REAL(dp) :: v_of_0, omega, alat, at(3,3), bg(3,3)
  INTEGER, ALLOCATABLE :: nl(:), ityp(:), nh(:), ofsbeta(:), isk(:), ngk(:)
  REAL(dp), ALLOCATABLE :: vrs(:,:), deeq(:,:,:,:), qq_at(:,:,:), xk(:,:), wk(:)
  CHARACTER(LEN=256) :: dumpdir, arg
  INTEGER :: iun1 = 31, iun2 = 32, ik, narg, i, nsel, sel(64), ikr, npw, nbytes_per_k
  INTEGER, ALLOCATABLE :: igk(:)
  REAL(dp), ALLOCATABLE :: g2kin(:), h_diag(:), s_diag(:), et(:)
  COMPLEX(dp), ALLOCATABLE :: vkb(:,:), evc(:,:), hpsi(:,:), spsi(:,:), hnl(:,:), hc(:,:), sc(:,:), hn(:,:)
  REAL(dp) :: e_h, e_s, e_nl, e_diag, e_sd
  LOGICAL :: check
  !
  narg = COMMAND_ARGUMENT_COUNT()
  IF (narg < 1) STOP 'usage: read_dump_check dumpdir [k ...]'
  CALL GET_COMMAND_ARGUMENT(1, dumpdir)
  OPEN(iun1, FILE=TRIM(dumpdir)//'/ham_global.bin', FORM='unformatted', ACCESS='stream', STATUS='old')
  READ(iun1) magic, nspin, nkstot, nbnd, npwx, nkb, nat, nsp, nhm, nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr, ngms, okvan, dbl
  IF (magic /= 20251113) STOP 'wrong magic number'
  READ(iun1) v_of_0, omega, alat, at, bg
  ALLOCATE(nl(ngms), vrs(nnr,nspin), ityp(nat), nh(nsp), ofsbeta(nat), deeq(nhm,nhm,nat,nspin), qq_at(nhm,nhm,nat))
  ALLOCATE(xk(3,nkstot), wk(nkstot), isk(nkstot), ngk(nkstot))
  READ(iun1) nl
  READ(iun1) vrs
  READ(iun1) ityp, nh, ofsbeta
  READ(iun1) deeq
  READ(iun1) qq_at
  READ(iun1) xk, wk, isk, ngk
  CLOSE(iun1)
  WRITE(*,'("global: nspin=",i1," nk=",i5," nbnd=",i3," npwx=",i5," nkb=",i3," nat=",i3," nsp=",i2," nhm=",i2,/,"        grid ",3i4," (padded ",3i4,") nnr=",i7," ngms=",i6," okvan=",i1," doublegrid=",i1)') &
       nspin, nkstot, nbnd, npwx, nkb, nat, nsp, nhm, nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr, ngms, okvan, dbl
  ! k-points to check
  nsel = 0
  DO i = 2, narg
     CALL GET_COMMAND_ARGUMENT(i, arg); nsel = nsel + 1; READ(arg,*) sel(nsel)
  END DO
  IF (nsel == 0) THEN
     sel(1:3) = [1, (nkstot+1)/2, nkstot]; nsel = 3
  END IF
  ALLOCATE(igk(npwx), g2kin(npwx), h_diag(npwx), s_diag(npwx), et(nbnd), vkb(npwx,MAX(nkb,1)))
  ALLOCATE(evc(npwx,nbnd), hpsi(npwx,nbnd), spsi(npwx,nbnd), hnl(npwx,nbnd), hc(npwx,nbnd), sc(npwx,nbnd), hn(npwx,nbnd))
  OPEN(iun2, FILE=TRIM(dumpdir)//'/ham_k.bin', FORM='unformatted', ACCESS='stream', STATUS='old')
  WRITE(*,'(/,"  k-point    npw  max|H psi-QE|  max|S psi-QE|  max|VNL psi-QE|  max|h_diag-QE|  max|s_diag-QE|")')
  DO ik = 1, nkstot
     READ(iun2) ikr, npw
     IF (ikr /= ik) STOP 'k-point order'
     READ(iun2) igk(1:npw)
     READ(iun2) g2kin(1:npw), h_diag(1:npw), s_diag(1:npw)
     READ(iun2) vkb(1:npw,1:nkb)
     READ(iun2) et(1:nbnd)
     READ(iun2) evc(1:npw,1:nbnd)
     READ(iun2) hpsi(1:npw,1:nbnd), spsi(1:npw,1:nbnd), hnl(1:npw,1:nbnd)
     check = ANY(sel(1:nsel) == ik)
     IF (.NOT. check) CYCLE
     CALL rebuild(ik, npw)
     e_h = MAXVAL(ABS(hc(1:npw,:) - hpsi(1:npw,:))); e_s = MAXVAL(ABS(sc(1:npw,:) - spsi(1:npw,:)))
     e_nl = MAXVAL(ABS(hn(1:npw,:) - hnl(1:npw,:)))
     CALL diag_check(ik, npw, e_diag, e_sd)
     WRITE(*,'(i8,i8,5es17.3)') ik, npw, e_h, e_s, e_nl, e_diag, e_sd
  END DO
  CLOSE(iun2)
CONTAINS
  SUBROUTINE diag_check(ik, npw, e_hd, e_sd)
    ! h_diag = g2kin + v_of_0 + sum D_ij Re(vkb_i vkb_j^*),  s_diag = 1 + sum qq_ij Re(vkb_i vkb_j^*)   (as QE's usnldiag)
    INTEGER, INTENT(IN) :: ik, npw
    REAL(dp), INTENT(OUT) :: e_hd, e_sd
    REAL(dp), ALLOCATABLE :: hd(:), sd(:)
    INTEGER :: ig, na, nt, ih, jh, spin
    REAL(dp) :: ar
    spin = isk(ik)
    ALLOCATE(hd(npw), sd(npw))
    hd = g2kin(1:npw) + v_of_0; sd = 1.d0
    DO na = 1, nat
       nt = ityp(na)
       DO ig = 1, npw; DO ih = 1, nh(nt); DO jh = 1, nh(nt)
          ar = DBLE(vkb(ig,ofsbeta(na)+ih) * CONJG(vkb(ig,ofsbeta(na)+jh)))
          hd(ig) = hd(ig) + deeq(ih,jh,na,spin) * ar
          sd(ig) = sd(ig) + qq_at(ih,jh,na) * ar
       END DO; END DO; END DO
    END DO
    e_hd = MAXVAL(ABS(hd - h_diag(1:npw))); e_sd = MAXVAL(ABS(sd - s_diag(1:npw)))
  END SUBROUTINE diag_check
  SUBROUTINE rebuild(ik, npw)
    INTEGER, INTENT(IN) :: ik, npw
    COMPLEX(dp), ALLOCATABLE :: e1(:), e2(:), e3(:), psir(:), becp(:,:), ps(:,:)
    REAL(dp), PARAMETER :: tpi = 6.283185307179586d0
    INTEGER :: ig, ib, ix, iy, iz, ii, jj, kk, n, na, nt, ih, jh, ikb, jkb, spin
    INTEGER, ALLOCATABLE :: gi(:), gj(:), gk(:)
    COMPLEX(dp) :: ph, acc
    spin = isk(ik)
    ALLOCATE(e1(0:nr1-1), e2(0:nr2-1), e3(0:nr3-1), psir(nnr), gi(npw), gj(npw), gk(npw))
    DO ix = 0, nr1-1; e1(ix) = CMPLX(COS(tpi*ix/nr1), SIN(tpi*ix/nr1), dp); END DO
    DO iy = 0, nr2-1; e2(iy) = CMPLX(COS(tpi*iy/nr2), SIN(tpi*iy/nr2), dp); END DO
    DO iz = 0, nr3-1; e3(iz) = CMPLX(COS(tpi*iz/nr3), SIN(tpi*iz/nr3), dp); END DO
    DO ig = 1, npw       ! Miller indices (mod grid size) of the plane waves from the smooth-grid index
       n = nl(igk(ig)) - 1
       gi(ig) = MOD(n, nr1x); gj(ig) = MOD(n/nr1x, nr2x); gk(ig) = n/(nr1x*nr2x)
    END DO
    ! local potential: psi(r) = sum_G c_G exp(+i G r); (V psi)(r); back: (1/N) sum_r f(r) exp(-i G r)
    hc = (0.d0,0.d0)
    DO ib = 1, nbnd
       psir = (0.d0,0.d0)
       DO iz = 0, nr3-1; DO iy = 0, nr2-1; DO ix = 0, nr1-1
          acc = (0.d0,0.d0)
          DO ig = 1, npw
             acc = acc + evc(ig,ib) * e1(MOD(gi(ig)*ix,nr1)) * e2(MOD(gj(ig)*iy,nr2)) * e3(MOD(gk(ig)*iz,nr3))
          END DO
          n = 1 + ix + nr1x*(iy + nr2x*iz)
          psir(n) = acc * vrs(n,spin)
       END DO; END DO; END DO
       DO ig = 1, npw
          acc = (0.d0,0.d0)
          DO iz = 0, nr3-1; DO iy = 0, nr2-1; DO ix = 0, nr1-1
             n = 1 + ix + nr1x*(iy + nr2x*iz)
             acc = acc + psir(n) * CONJG(e1(MOD(gi(ig)*ix,nr1)) * e2(MOD(gj(ig)*iy,nr2)) * e3(MOD(gk(ig)*iz,nr3)))
          END DO; END DO; END DO
          hc(ig,ib) = acc / REAL(nr1*nr2*nr3, dp) + g2kin(ig) * evc(ig,ib)
       END DO
    END DO
    ! non local part: vkb D vkb^H psi  (D per atom, projectors of atom na start at ofsbeta(na)+1)
    ALLOCATE(becp(nkb,nbnd), ps(nkb,nbnd))
    becp = (0.d0,0.d0); ps = (0.d0,0.d0)
    DO ib = 1, nbnd; DO ikb = 1, nkb
       becp(ikb,ib) = SUM(CONJG(vkb(1:npw,ikb)) * evc(1:npw,ib))
    END DO; END DO
    DO na = 1, nat
       nt = ityp(na)
       DO ib = 1, nbnd; DO ih = 1, nh(nt); DO jh = 1, nh(nt)
          ps(ofsbeta(na)+ih,ib) = ps(ofsbeta(na)+ih,ib) + deeq(ih,jh,na,spin) * becp(ofsbeta(na)+jh,ib)
       END DO; END DO; END DO
    END DO
    hn = (0.d0,0.d0)
    DO ib = 1, nbnd; DO ikb = 1, nkb; hn(1:npw,ib) = hn(1:npw,ib) + vkb(1:npw,ikb) * ps(ikb,ib); END DO; END DO
    hc(1:npw,:) = hc(1:npw,:) + hn(1:npw,1:nbnd)
    ! overlap: psi + vkb qq vkb^H psi
    ps = (0.d0,0.d0)
    DO na = 1, nat
       nt = ityp(na)
       DO ib = 1, nbnd; DO ih = 1, nh(nt); DO jh = 1, nh(nt)
          ps(ofsbeta(na)+ih,ib) = ps(ofsbeta(na)+ih,ib) + qq_at(ih,jh,na) * becp(ofsbeta(na)+jh,ib)
       END DO; END DO; END DO
    END DO
    sc = evc
    DO ib = 1, nbnd; DO ikb = 1, nkb; sc(1:npw,ib) = sc(1:npw,ib) + vkb(1:npw,ikb) * ps(ikb,ib); END DO; END DO
    DEALLOCATE(e1, e2, e3, psir, gi, gj, gk, becp, ps)
  END SUBROUTINE rebuild
END PROGRAM read_dump_check
