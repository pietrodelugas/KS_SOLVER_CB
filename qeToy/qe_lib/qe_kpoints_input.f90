!
! qe_kpoints_input: k-point sets for qeToy. No symmetry: every k-point of the set is computed directly.
!   kmode = 'dump' : the k-points of the dumped run (nothing to do)
!   kmode = 'list' : read from the input unit:  n  then n lines   kx ky kz weight        (cartesian, units 2pi/alat)
!   kmode = 'path' : read from the input unit:  nv then nv lines  kx ky kz npts           (cartesian, units 2pi/alat)
!                    npts points from each vertex to the next one (the last vertex is added once), weight 1
!   kmode = 'mp'   : full Monkhorst-Pack grid nk(1) x nk(2) x nk(3), shifted by shift(i)/2 grid steps, folded in (-1/2,1/2]
! For LSDA every k-point is repeated for the two spins (all spin up first, then all spin down, as QE).
! It then sets nkstot, xk, wk, isk, ngk and npwx for the new set (the data of the dump setup are needed).
!
MODULE qe_kpoints_input
  USE qe_ham_module
  USE qe_kpoint_module
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  PUBLIC :: qe_set_kpoints
CONTAINS
  SUBROUTINE qe_set_kpoints(kmode, nk, shift, iun)
    CHARACTER(LEN=*), INTENT(IN) :: kmode
    INTEGER, INTENT(IN) :: nk(3), shift(3), iun
    REAL(dp), ALLOCATABLE :: xl(:,:), wl(:), v(:,:)
    INTEGER, ALLOCATABLE :: np(:)
    INTEGER :: n, nv, i, j, k, iv, ip, is, nl
    REAL(dp) :: x(3)
    IF (TRIM(kmode) == 'dump') RETURN
    IF (.NOT. has_setup) STOP 'qe_set_kpoints: the dump has no setup data (dumper v2): only kmode=dump is possible'
    SELECT CASE (TRIM(kmode))
    CASE ('list')
       READ(iun,*) n
       ALLOCATE(xl(3,n), wl(n))
       DO i = 1, n; READ(iun,*) xl(:,i), wl(i); END DO
    CASE ('path')
       READ(iun,*) nv
       ALLOCATE(v(3,nv), np(nv))
       DO i = 1, nv; READ(iun,*) v(:,i), np(i); END DO
       n = 1 + SUM(np(1:nv-1))
       ALLOCATE(xl(3,n), wl(n)); wl = 1.d0
       n = 0
       DO iv = 1, nv-1
          DO ip = 0, np(iv)-1
             n = n + 1; xl(:,n) = v(:,iv) + (v(:,iv+1)-v(:,iv)) * REAL(ip,dp)/REAL(np(iv),dp)
          END DO
       END DO
       n = n + 1; xl(:,n) = v(:,nv)
    CASE ('mp')
       IF (ANY(nk < 1)) STOP 'qe_set_kpoints: nk1, nk2, nk3 must be positive'
       n = nk(1)*nk(2)*nk(3)
       ALLOCATE(xl(3,n), wl(n)); wl = 1.d0/REAL(n,dp)
       nl = 0
       DO i = 1, nk(1); DO j = 1, nk(2); DO k = 1, nk(3)
          nl = nl + 1
          x = [ REAL(i-1,dp)/nk(1) + REAL(shift(1),dp)/2/nk(1), REAL(j-1,dp)/nk(2) + REAL(shift(2),dp)/2/nk(2), &
                REAL(k-1,dp)/nk(3) + REAL(shift(3),dp)/2/nk(3) ]
          x = x - NINT(x)
          xl(:,nl) = MATMUL(bg, x)
       END DO; END DO; END DO
    CASE DEFAULT
       STOP "qe_set_kpoints: kmode must be 'dump', 'list', 'path' or 'mp'"
    END SELECT
    n = SIZE(xl,2)
    DEALLOCATE(xk, wk, isk, ngk)
    nkstot = n * nspin
    ALLOCATE(xk(3,nkstot), wk(nkstot), isk(nkstot), ngk(nkstot))
    DO is = 1, nspin
       xk(:,(is-1)*n+1:is*n) = xl; wk((is-1)*n+1:is*n) = wl; isk((is-1)*n+1:is*n) = is
    END DO
    DO i = 1, nkstot
       ngk(i) = qe_count_pw(xk(:,i))
    END DO
    npwx = MAXVAL(ngk)
    kcompute = .TRUE.
  END SUBROUTINE qe_set_kpoints
END MODULE qe_kpoints_input
