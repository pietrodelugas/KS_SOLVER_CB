!
! dump_ham: reads a finished QE run (as pp_example does) and writes, for every k-point, everything needed to rebuild
! H|psi> and S|psi> outside QE, together with QE's own results for a set of test vectors.
!
! Input (stdin): &inputpp  outdir='./tmp', prefix='pwscf', dumpdir='./dump' /
!
! Files (Fortran unformatted, access='stream', native endianness, 32 bit integers, 64 bit reals):
!   <dumpdir>/ham_global.bin   written once
!   <dumpdir>/ham_k.bin        one block per k-point, in the order of the k-points of the run (spin up then down for LSDA)
! The layout is documented in README.md and read by qeToy/ (to be written).
!
PROGRAM dump_ham
  USE io_global,  ONLY : ionode, ionode_id
  USE io_files,   ONLY : tmp_dir, prefix
  USE mp_global,  ONLY : mp_startup
  USE mp_images,  ONLY : intra_image_comm
  USE mp,         ONLY : mp_bcast
  USE environment,ONLY : environment_start, environment_end
  IMPLICIT NONE
  LOGICAL :: needwf = .true.
  INTEGER :: ios
  CHARACTER(LEN=256) :: outdir, dumpdir
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  NAMELIST / inputpp / outdir, prefix, dumpdir
  CALL mp_startup ( )
  CALL environment_start ( 'DUMPHAM' )
  IF ( ionode )  CALL input_from_file ( )
  prefix = 'pwscf'
  dumpdir = './dump'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  IF ( ionode )  THEN
     READ (5, inputpp, iostat = ios)
     tmp_dir = trimcheck ( outdir )
  ENDIF
  CALL mp_bcast (ios, ionode_id, intra_image_comm)
  IF ( ios /= 0) CALL errore ('dump_ham', 'reading inputpp namelist', abs(ios))
  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )
  CALL mp_bcast( dumpdir, ionode_id, intra_image_comm )
  !
  CALL read_file_new ( needwf )
  CALL dump_all ( TRIM(dumpdir) )
  CALL environment_end( )
  STOP
END PROGRAM dump_ham
!
!-----------------------------------------------------------------------
SUBROUTINE dump_all ( dumpdir )
  !-----------------------------------------------------------------------
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode
  USE constants,  ONLY : rytoev
  USE becmod,     ONLY : becp, calbec, allocate_bec_type
  USE cell_base,  ONLY : alat, at, bg, omega
  USE fft_base,   ONLY : dfftp, dffts
  USE klist,      ONLY : xk, wk, nks, nkstot, igk_k, ngk
  USE lsda_mod,   ONLY : nspin, isk, current_spin
  USE io_files,   ONLY : restart_dir
  USE scf,        ONLY : vrs, vltot, v, kedtau, v_of_0
  USE gvecs,      ONLY : doublegrid, ngms
  USE uspp,       ONLY : nkb, vkb, deeq, qq_at, okvan, ofsbeta
  USE uspp_param, ONLY : nh, nhm
  USE uspp_init,  ONLY : init_us_2
  USE ions_base,  ONLY : nat, ityp, nsp => nsp
  USE wvfct,      ONLY : npwx, nbnd, current_k, et, g2kin
  USE noncollin_module, ONLY : noncolin, npol
  USE control_flags,    ONLY : gamma_only
  USE wavefunctions,    ONLY : evc
  USE pw_restart_new,   ONLY : read_collected_wfc
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: dumpdir
  INTEGER, PARAMETER :: magic = 20251113
  COMPLEX(DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:), hnl(:,:)
  REAL(DP), ALLOCATABLE :: h_diag(:,:), s_diag(:,:)
  INTEGER :: ik, npw, iun1, iun2, ierr, i
  !
  IF ( noncolin )   CALL errore( 'dump_ham', 'non-collinear calculations are not supported', 1 )
  IF ( gamma_only ) CALL errore( 'dump_ham', 'gamma_only calculations are not supported', 1 )
  CALL execute_command_line( 'mkdir -p '//dumpdir )
  !
  ALLOCATE( hpsi(npwx,nbnd), spsi(npwx,nbnd), hnl(npwx,nbnd) )
  ALLOCATE( h_diag(npwx,1), s_diag(npwx,1) )
  CALL allocate_bec_type ( nkb, nbnd, becp )
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  !
  ! ---------------------------------------------------------------- global data
  iun1 = 31
  OPEN( UNIT=iun1, FILE=dumpdir//'/ham_global.bin', FORM='unformatted', ACCESS='stream', STATUS='replace' )
  WRITE(iun1) magic, nspin, nkstot, nbnd, npwx, nkb, nat, nsp, nhm, &
              dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x, dffts%nnr, ngms, &
              MERGE(1,0,okvan), MERGE(1,0,doublegrid)
  WRITE(iun1) v_of_0, omega, alat, at, bg
  WRITE(iun1) dffts%nl(1:ngms)                  ! smooth-grid index of the G vector number ig (ig = igk_k(:,ik))
  WRITE(iun1) vrs(1:dffts%nnr,1:nspin)          ! local potential on the smooth grid, per spin
  WRITE(iun1) ityp(1:nat), nh(1:nsp), ofsbeta(1:nat)
  WRITE(iun1) deeq(1:nhm,1:nhm,1:nat,1:nspin)   ! D_ij per atom and spin (= dvan for norm-conserving)
  WRITE(iun1) qq_at(1:nhm,1:nhm,1:nat)          ! overlap integrals (zero for norm-conserving)
  WRITE(iun1) xk(1:3,1:nkstot), wk(1:nkstot), isk(1:nkstot), ngk(1:nkstot)
  CLOSE(iun1)
  !
  ! ------------------------------------------------------------- per k-point data
  iun2 = 32
  OPEN( UNIT=iun2, FILE=dumpdir//'/ham_k.bin', FORM='unformatted', ACCESS='stream', STATUS='replace' )
  DO ik = 1, nkstot
     CALL read_collected_wfc ( restart_dir() , ik, evc )
     npw = ngk(ik)
     current_k = ik
     current_spin = isk(ik)
     CALL init_us_2 ( npw, igk_k(1,ik), xk(1,ik), vkb )
     CALL calbec ( npw, vkb, evc, becp )
     CALL g2_kin ( ik )
     CALL usnldiag ( npw, 1, h_diag, s_diag )
     !
     CALL h_psi( npwx, npw, nbnd, evc, hpsi )     ! H psi with QE's own routine
     CALL s_psi( npwx, npw, nbnd, evc, spsi )
     hnl = (0.d0,0.d0)
     CALL add_vuspsi( npwx, npw, nbnd, hnl )      ! only the non local part (becp is that of evc)
     !
     WRITE(iun2) ik, npw
     WRITE(iun2) igk_k(1:npw,ik)
     WRITE(iun2) g2kin(1:npw), h_diag(1:npw,1), s_diag(1:npw,1)
     WRITE(iun2) vkb(1:npw,1:nkb)
     WRITE(iun2) et(1:nbnd,ik)                    ! eigenvalues of the run (Ry)
     WRITE(iun2) evc(1:npw,1:nbnd)                ! converged wavefunctions of the run, used as test vectors
     WRITE(iun2) hpsi(1:npw,1:nbnd), spsi(1:npw,1:nbnd), hnl(1:npw,1:nbnd)
  END DO
  CLOSE(iun2)
  IF (ionode) THEN
     PRINT '(/,5x,"dump_ham: wrote ",i6," k-points, nbnd=",i4,", npwx=",i6,", nkb=",i4,", smooth grid ",3i5)', &
           nkstot, nbnd, npwx, nkb, dffts%nr1, dffts%nr2, dffts%nr3
  END IF
END SUBROUTINE dump_all
