!
! qe_ham_module: Hamiltonian data of a Quantum ESPRESSO run, read from the files of the dumper (see DUMP_FORMAT.md).
! No QE code is used. ham_k.bin is read by random access: the byte offset of every k-point block follows from ngk.
!
MODULE qe_ham_module
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  INTEGER, PARAMETER :: magic_number = 20251113
  INTEGER, PUBLIC :: nspin, nkstot, nbnd, npwx, nkb, nat, nsp, nhm, nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr, ngms
  LOGICAL, PUBLIC :: okvan, doublegrid
  REAL(dp), PUBLIC :: v_of_0, omega, alat, at(3,3), bg(3,3)
  INTEGER, ALLOCATABLE, PUBLIC :: nl(:), ityp(:), nh(:), ofsbeta(:), isk(:), ngk(:)
  REAL(dp), ALLOCATABLE, PUBLIC :: vrs(:,:), deeq(:,:,:,:), qq_at(:,:,:), xk(:,:), wk(:)
  LOGICAL, PUBLIC :: dump_has_tests = .TRUE.   ! test vectors (evc, hpsi, spsi, hnl) present in ham_k.bin
  INTEGER(8), ALLOCATABLE :: koffset(:)         ! first byte of each k-point block
  INTEGER :: iun_k = -1
  PUBLIC :: qe_read_global, qe_open_k, qe_read_k, qe_close_k
CONTAINS
  SUBROUTINE qe_read_global(dumpdir)
    CHARACTER(LEN=*), INTENT(IN) :: dumpdir
    INTEGER :: iun, magic, ok, dbl, ik
    INTEGER(8) :: pos, sz
    OPEN(NEWUNIT=iun, FILE=TRIM(dumpdir)//'/ham_global.bin', FORM='unformatted', ACCESS='stream', STATUS='old')
    READ(iun) magic, nspin, nkstot, nbnd, npwx, nkb, nat, nsp, nhm, nr1, nr2, nr3, nr1x, nr2x, nr3x, nnr, ngms, ok, dbl
    IF (magic /= magic_number) STOP 'qe_read_global: wrong magic number'
    okvan = ok /= 0; doublegrid = dbl /= 0
    READ(iun) v_of_0, omega, alat, at, bg
    ALLOCATE(nl(ngms), vrs(nnr,nspin), ityp(nat), nh(nsp), ofsbeta(nat), deeq(nhm,nhm,nat,nspin), qq_at(nhm,nhm,nat))
    ALLOCATE(xk(3,nkstot), wk(nkstot), isk(nkstot), ngk(nkstot))
    READ(iun) nl
    READ(iun) vrs
    READ(iun) ityp, nh, ofsbeta
    READ(iun) deeq
    READ(iun) qq_at
    READ(iun) xk, wk, isk, ngk
    CLOSE(iun)
    npwx = MAXVAL(ngk)   ! consistent with the file even if the header was written with a larger value
    ! byte offsets of the k-point blocks, with the test vectors if the file size says so
    OPEN(NEWUNIT=iun, FILE=TRIM(dumpdir)//'/ham_k.bin', FORM='unformatted', ACCESS='stream', STATUS='old')
    INQUIRE(UNIT=iun, SIZE=sz); CLOSE(iun)
    ALLOCATE(koffset(nkstot+1))
    dump_has_tests = (sz == total_size(.TRUE.))
    IF (.NOT. dump_has_tests .AND. sz /= total_size(.FALSE.)) STOP 'qe_read_global: size of ham_k.bin does not match'
    pos = 1
    DO ik = 1, nkstot
       koffset(ik) = pos; pos = pos + block_size(ngk(ik), dump_has_tests)
    END DO
    koffset(nkstot+1) = pos
  END SUBROUTINE qe_read_global
  !
  INTEGER(8) FUNCTION block_size(npw, tests)
    INTEGER, INTENT(IN) :: npw
    LOGICAL, INTENT(IN) :: tests
    INTEGER(8) :: n
    n = npw
    block_size = 8 + 4*n + 3*8*n + 16*n*nkb + 8*nbnd
    IF (tests) block_size = block_size + 4*16*n*nbnd
  END FUNCTION block_size
  INTEGER(8) FUNCTION total_size(tests)
    LOGICAL, INTENT(IN) :: tests
    INTEGER :: ik
    total_size = 0
    DO ik = 1, nkstot; total_size = total_size + block_size(ngk(ik), tests); END DO
  END FUNCTION total_size
  !
  SUBROUTINE qe_open_k(dumpdir)
    CHARACTER(LEN=*), INTENT(IN) :: dumpdir
    OPEN(NEWUNIT=iun_k, FILE=TRIM(dumpdir)//'/ham_k.bin', FORM='unformatted', ACCESS='stream', STATUS='old', ACTION='read')
  END SUBROUTINE qe_open_k
  SUBROUTINE qe_close_k()
    IF (iun_k >= 0) CLOSE(iun_k); iun_k = -1
  END SUBROUTINE qe_close_k
  !
  ! Reads the block of k-point ik. Arrays must have first dimension >= ngk(ik); optional ones are the test vectors.
  SUBROUTINE qe_read_k(ik, igk, g2kin, h_diag, s_diag, vkb, et, evc, hpsi, spsi, hnl)
    INTEGER, INTENT(IN) :: ik
    INTEGER, INTENT(OUT) :: igk(:)
    REAL(dp), INTENT(OUT) :: g2kin(:), h_diag(:), s_diag(:), et(:)
    COMPLEX(dp), INTENT(OUT) :: vkb(:,:)
    COMPLEX(dp), INTENT(OUT), OPTIONAL :: evc(:,:), hpsi(:,:), spsi(:,:), hnl(:,:)
    INTEGER :: ikr, npw, npw_file
    npw = ngk(ik)
    READ(iun_k, POS=koffset(ik)) ikr, npw_file
    IF (ikr /= ik .OR. npw_file /= npw) STOP 'qe_read_k: k-point block mismatch'
    READ(iun_k) igk(1:npw)
    READ(iun_k) g2kin(1:npw), h_diag(1:npw), s_diag(1:npw)
    READ(iun_k) vkb(1:npw,1:nkb)
    READ(iun_k) et(1:nbnd)
    IF (PRESENT(evc)) THEN
       IF (.NOT. dump_has_tests) STOP 'qe_read_k: test vectors requested but not in the dump'
       READ(iun_k) evc(1:npw,1:nbnd)
       READ(iun_k) hpsi(1:npw,1:nbnd), spsi(1:npw,1:nbnd), hnl(1:npw,1:nbnd)
    END IF
  END SUBROUTINE qe_read_k
END MODULE qe_ham_module
