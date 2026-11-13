!
! qe_dump_h5: HDF5 version of the dump files (layout in DUMP_FORMAT.md). One file, ham.h5:
!   root attributes : magic, version, nspin, nkstot, nbnd, npwx, nkb, nat, nsp, nhm, nr1..nr3, nr1x..nr3x, nnr, ngms,
!                     okvan, doublegrid, has_tests, v_of_0, omega, alat, at(9), bg(9)
!   root datasets   : nl, vrs, ityp, nh, ofsbeta, deeq, qq_at, xk, wk, isk, ngk
!                     (deeq and qq_at are stored with their last two indices merged: (nhm,nhm,nat*nspin) and (nhm,nhm,nat))
!   groups k00001.. : igk, g2kin, h_diag, s_diag, vkb, et and, if has_tests, evc, hpsi, spsi, hnl
! Complex arrays are stored by qeh5 with a leading dimension 2 (real, imaginary).
! The same source file is used by the dumper (dump_ham.f90, in the PP directory of Quantum ESPRESSO) and by the qeToy miniapp.
! It needs the module qeh5_base_module (Modules/qeh5_module.f90 of QE, or UtilXlib/qeh5_module.f90 of KS_SOLVER_CB),
! and is empty unless __HDF5 is defined.
!
MODULE qe_dump_h5
#if defined(__HDF5)
  USE qeh5_base_module
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  INTEGER, PARAMETER, PUBLIC :: nihead = 19
  TYPE(qeh5_file) :: hf
  INTERFACE put
     MODULE PROCEDURE put_i1, put_i2, put_i3, put_r1, put_r2, put_r3, put_c2
  END INTERFACE
  INTERFACE get
     MODULE PROCEDURE get_i1, get_i2, get_i3, get_r1, get_r2, get_r3, get_c2
  END INTERFACE
  PUBLIC :: dump_h5_create, dump_h5_write_global, dump_h5_write_k, dump_h5_close
  PUBLIC :: dump_h5_write_setup, dump_h5_has_setup, dump_h5_read_setup_head, dump_h5_read_setup
  PUBLIC :: dump_h5_open_read, dump_h5_read_head, dump_h5_read_global, dump_h5_read_k
  CHARACTER(LEN=*), PARAMETER :: ihead_names(nihead) = [CHARACTER(LEN=10) :: 'magic', 'nspin', 'nkstot', 'nbnd', 'npwx', 'nkb', &
       'nat', 'nsp', 'nhm', 'nr1', 'nr2', 'nr3', 'nr1x', 'nr2x', 'nr3x', 'nnr', 'ngms', 'okvan', 'doublegrid']
CONTAINS
  CHARACTER(LEN=6) FUNCTION kname(ik)
    INTEGER, INTENT(IN) :: ik
    WRITE(kname,'("k",i5.5)') ik
  END FUNCTION kname
  !
  SUBROUTINE dump_h5_create(fname)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    CALL initialize_hdf5()
    CALL qeh5_openfile(hf, TRIM(fname), 'write')
  END SUBROUTINE dump_h5_create
  SUBROUTINE dump_h5_close()
    CALL qeh5_close(hf)
    CALL finalize_hdf5()
  END SUBROUTINE dump_h5_close
  !
  ! ihead: the 19 integers in the order of ihead_names (as the header of ham_global.bin)
  SUBROUTINE dump_h5_write_global(ihead, has_tests, v_of_0, omega, alat, at, bg, nl, vrs, ityp, nh, ofsbeta, deeq, qq_at, &
                                  xk, wk, isk, ngk)
    INTEGER, INTENT(IN) :: ihead(nihead), has_tests
    REAL(dp), INTENT(IN) :: v_of_0, omega, alat, at(3,3), bg(3,3)
    INTEGER, INTENT(IN) :: nl(:), ityp(:), nh(:), ofsbeta(:), isk(:), ngk(:)
    REAL(dp), INTENT(IN) :: vrs(:,:), deeq(:,:,:,:), qq_at(:,:,:), xk(:,:), wk(:)
    INTEGER :: i
    DO i = 1, nihead
       CALL qeh5_add_attribute(hf%id, TRIM(ihead_names(i)), ihead(i))
    END DO
    CALL qeh5_add_attribute(hf%id, 'version', 1)
    CALL qeh5_add_attribute(hf%id, 'has_tests', has_tests)
    CALL qeh5_add_attribute(hf%id, 'v_of_0', v_of_0)
    CALL qeh5_add_attribute(hf%id, 'omega', omega)
    CALL qeh5_add_attribute(hf%id, 'alat', alat)
    CALL qeh5_add_attribute(hf%id, 'at', RESHAPE(at, [9]), 1, [9])
    CALL qeh5_add_attribute(hf%id, 'bg', RESHAPE(bg, [9]), 1, [9])
    CALL put(hf, 'nl', nl); CALL put(hf, 'vrs', vrs)
    CALL put(hf, 'ityp', ityp); CALL put(hf, 'nh', nh); CALL put(hf, 'ofsbeta', ofsbeta)
    CALL put(hf, 'deeq', RESHAPE(deeq, [SIZE(deeq,1), SIZE(deeq,2), SIZE(deeq,3)*SIZE(deeq,4)]))
    CALL put(hf, 'qq_at', qq_at)
    CALL put(hf, 'xk', xk); CALL put(hf, 'wk', wk); CALL put(hf, 'isk', isk); CALL put(hf, 'ngk', ngk)
  END SUBROUTINE dump_h5_write_global
  !
  ! Data to rebuild the k dependent part of the Hamiltonian for arbitrary k-points (group 'setup', see DUMP_FORMAT.md):
  ! ihs = [lmaxkb, nbetam, nqx], rhs = [ecutwfc, tpiba, dq]
  SUBROUTINE dump_h5_write_setup(ihs, rhs, tau, mill, tab, indv, nhtol, nhtolm)
    INTEGER, INTENT(IN) :: ihs(3), mill(:,:), indv(:,:), nhtol(:,:), nhtolm(:,:)
    REAL(dp), INTENT(IN) :: rhs(3), tau(:,:), tab(:,:,:)
    TYPE(qeh5_file) :: g
    CALL qeh5_open_group(hf, 'setup', g)
    CALL qeh5_add_attribute(g%id, 'lmaxkb', ihs(1)); CALL qeh5_add_attribute(g%id, 'nbetam', ihs(2))
    CALL qeh5_add_attribute(g%id, 'nqx', ihs(3))
    CALL qeh5_add_attribute(g%id, 'ecutwfc', rhs(1)); CALL qeh5_add_attribute(g%id, 'tpiba', rhs(2))
    CALL qeh5_add_attribute(g%id, 'dq', rhs(3))
    CALL put(g, 'tau', tau); CALL put(g, 'mill', mill); CALL put(g, 'tab', tab)
    CALL put(g, 'indv', indv); CALL put(g, 'nhtol', nhtol); CALL put(g, 'nhtolm', nhtolm)
    CALL qeh5_close(g)
  END SUBROUTINE dump_h5_write_setup
  !
  LOGICAL FUNCTION dump_h5_has_setup()
    USE hdf5, ONLY : h5lexists_f
    INTEGER :: ierr
    LOGICAL :: ex
    CALL h5lexists_f(hf%id, 'setup', ex, ierr)
    dump_h5_has_setup = ex .AND. ierr == 0
  END FUNCTION dump_h5_has_setup
  !
  SUBROUTINE dump_h5_read_setup_head(ihs, rhs)
    INTEGER, INTENT(OUT) :: ihs(3)
    REAL(dp), INTENT(OUT) :: rhs(3)
    TYPE(qeh5_file) :: g
    CALL qeh5_open_group(hf, 'setup', g)
    CALL qeh5_read_attribute(g%id, 'lmaxkb', ihs(1)); CALL qeh5_read_attribute(g%id, 'nbetam', ihs(2))
    CALL qeh5_read_attribute(g%id, 'nqx', ihs(3))
    CALL qeh5_read_attribute(g%id, 'ecutwfc', rhs(1)); CALL qeh5_read_attribute(g%id, 'tpiba', rhs(2))
    CALL qeh5_read_attribute(g%id, 'dq', rhs(3))
    CALL qeh5_close(g)
  END SUBROUTINE dump_h5_read_setup_head
  SUBROUTINE dump_h5_read_setup(tau, mill, tab, indv, nhtol, nhtolm)
    INTEGER, INTENT(OUT) :: mill(:,:), indv(:,:), nhtol(:,:), nhtolm(:,:)
    REAL(dp), INTENT(OUT) :: tau(:,:), tab(:,:,:)
    TYPE(qeh5_file) :: g
    CALL qeh5_open_group(hf, 'setup', g)
    CALL get(g, 'tau', tau); CALL get(g, 'mill', mill); CALL get(g, 'tab', tab)
    CALL get(g, 'indv', indv); CALL get(g, 'nhtol', nhtol); CALL get(g, 'nhtolm', nhtolm)
    CALL qeh5_close(g)
  END SUBROUTINE dump_h5_read_setup
  !
  SUBROUTINE dump_h5_write_k(ik, npw, nbnd, nkb, igk, g2kin, h_diag, s_diag, vkb, et, evc, hpsi, spsi, hnl)
    INTEGER, INTENT(IN) :: ik, npw, nbnd, nkb
    INTEGER, INTENT(IN) :: igk(:)
    REAL(dp), INTENT(IN) :: g2kin(:), h_diag(:), s_diag(:), et(:)
    COMPLEX(dp), INTENT(IN) :: vkb(:,:)
    COMPLEX(dp), INTENT(IN), OPTIONAL :: evc(:,:), hpsi(:,:), spsi(:,:), hnl(:,:)
    TYPE(qeh5_file) :: g
    CALL qeh5_open_group(hf, kname(ik), g)
    CALL put(g, 'igk', igk(1:npw)); CALL put(g, 'g2kin', g2kin(1:npw)); CALL put(g, 'h_diag', h_diag(1:npw))
    CALL put(g, 's_diag', s_diag(1:npw)); CALL put(g, 'vkb', vkb(1:npw,1:nkb)); CALL put(g, 'et', et(1:nbnd))
    IF (PRESENT(evc)) THEN
       CALL put(g, 'evc', evc(1:npw,1:nbnd)); CALL put(g, 'hpsi', hpsi(1:npw,1:nbnd))
       CALL put(g, 'spsi', spsi(1:npw,1:nbnd)); CALL put(g, 'hnl', hnl(1:npw,1:nbnd))
    END IF
    CALL qeh5_close(g)
  END SUBROUTINE dump_h5_write_k
  !
  SUBROUTINE dump_h5_open_read(fname)
    CHARACTER(LEN=*), INTENT(IN) :: fname
    CALL initialize_hdf5()
    CALL qeh5_openfile(hf, TRIM(fname), 'read')
  END SUBROUTINE dump_h5_open_read
  !
  SUBROUTINE dump_h5_read_head(ihead, has_tests, v_of_0, omega, alat, at, bg)
    INTEGER, INTENT(OUT) :: ihead(nihead), has_tests
    REAL(dp), INTENT(OUT) :: v_of_0, omega, alat, at(3,3), bg(3,3)
    REAL(dp) :: t(9)
    INTEGER :: i
    DO i = 1, nihead
       CALL qeh5_read_attribute(hf%id, TRIM(ihead_names(i)), ihead(i))
    END DO
    CALL qeh5_read_attribute(hf%id, 'has_tests', has_tests)
    CALL qeh5_read_attribute(hf%id, 'v_of_0', v_of_0)
    CALL qeh5_read_attribute(hf%id, 'omega', omega)
    CALL qeh5_read_attribute(hf%id, 'alat', alat)
    CALL qeh5_read_attribute(hf%id, 'at', t, 1, [9]); at = RESHAPE(t, [3,3])
    CALL qeh5_read_attribute(hf%id, 'bg', t, 1, [9]); bg = RESHAPE(t, [3,3])
  END SUBROUTINE dump_h5_read_head
  !
  SUBROUTINE dump_h5_read_global(nl, vrs, ityp, nh, ofsbeta, deeq, qq_at, xk, wk, isk, ngk)
    INTEGER, INTENT(OUT) :: nl(:), ityp(:), nh(:), ofsbeta(:), isk(:), ngk(:)
    REAL(dp), INTENT(OUT) :: vrs(:,:), deeq(:,:,:,:), qq_at(:,:,:), xk(:,:), wk(:)
    REAL(dp), ALLOCATABLE :: d3(:,:,:)
    ALLOCATE(d3(SIZE(deeq,1), SIZE(deeq,2), SIZE(deeq,3)*SIZE(deeq,4)))
    CALL get(hf, 'nl', nl); CALL get(hf, 'vrs', vrs)
    CALL get(hf, 'ityp', ityp); CALL get(hf, 'nh', nh); CALL get(hf, 'ofsbeta', ofsbeta)
    CALL get(hf, 'deeq', d3); deeq = RESHAPE(d3, SHAPE(deeq))
    CALL get(hf, 'qq_at', qq_at)
    CALL get(hf, 'xk', xk); CALL get(hf, 'wk', wk); CALL get(hf, 'isk', isk); CALL get(hf, 'ngk', ngk)
  END SUBROUTINE dump_h5_read_global
  !
  SUBROUTINE dump_h5_read_k(ik, npw, nbnd, nkb, igk, g2kin, h_diag, s_diag, vkb, et, evc, hpsi, spsi, hnl)
    INTEGER, INTENT(IN) :: ik, npw, nbnd, nkb
    INTEGER, INTENT(OUT) :: igk(:)
    REAL(dp), INTENT(OUT) :: g2kin(:), h_diag(:), s_diag(:), et(:)
    COMPLEX(dp), INTENT(OUT) :: vkb(:,:)
    COMPLEX(dp), INTENT(OUT), OPTIONAL :: evc(:,:), hpsi(:,:), spsi(:,:), hnl(:,:)
    TYPE(qeh5_file) :: g
    INTEGER, ALLOCATABLE :: ti(:)
    REAL(dp), ALLOCATABLE :: tr(:)
    COMPLEX(dp), ALLOCATABLE :: tc(:,:)
    CALL qeh5_open_group(hf, kname(ik), g)
    ALLOCATE(ti(npw), tr(npw), tc(npw,MAX(nkb,nbnd)))
    CALL get(g, 'igk', ti); igk(1:npw) = ti
    CALL get(g, 'g2kin', tr); g2kin(1:npw) = tr
    CALL get(g, 'h_diag', tr); h_diag(1:npw) = tr
    CALL get(g, 's_diag', tr); s_diag(1:npw) = tr
    DEALLOCATE(tc); ALLOCATE(tc(npw,nkb)); CALL get(g, 'vkb', tc); vkb(1:npw,1:nkb) = tc
    DEALLOCATE(tr); ALLOCATE(tr(nbnd)); CALL get(g, 'et', tr); et(1:nbnd) = tr
    IF (PRESENT(evc)) THEN
       DEALLOCATE(tc); ALLOCATE(tc(npw,nbnd))
       CALL get(g, 'evc', tc); evc(1:npw,1:nbnd) = tc
       CALL get(g, 'hpsi', tc); hpsi(1:npw,1:nbnd) = tc
       CALL get(g, 'spsi', tc); spsi(1:npw,1:nbnd) = tc
       CALL get(g, 'hnl', tc); hnl(1:npw,1:nbnd) = tc
    END IF
    CALL qeh5_close(g)
  END SUBROUTINE dump_h5_read_k
  !
  ! ---- one dataset in/out; the arrays are copied to contiguous temporaries because qeh5 takes explicit-shape dummies
  SUBROUTINE put_i1(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: a(:)
    INTEGER, ALLOCATABLE :: t(:)
    TYPE(qeh5_dataset) :: d
    t = a
    CALL qeh5_set_space(d, t(1), 1, [SIZE(t)])
    CALL qeh5_open_dataset(f, d, action='write', name=name)
    CALL qeh5_write_dataset(t, d)
    CALL qeh5_close(d)
  END SUBROUTINE put_i1
  SUBROUTINE put_i2(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: a(:,:)
    INTEGER, ALLOCATABLE :: t(:,:)
    TYPE(qeh5_dataset) :: d
    t = a
    CALL qeh5_set_space(d, t(1,1), 2, SHAPE(t))
    CALL qeh5_open_dataset(f, d, action='write', name=name)
    CALL qeh5_write_dataset(t, d)
    CALL qeh5_close(d)
  END SUBROUTINE put_i2
  SUBROUTINE get_i2(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: a(:,:)
    INTEGER, ALLOCATABLE :: t(:,:)
    TYPE(qeh5_dataset) :: d
    ALLOCATE(t(SIZE(a,1),SIZE(a,2)))
    CALL qeh5_open_dataset(f, d, action='read', name=name)
    CALL qeh5_read_dataset(t, d)
    CALL qeh5_close(d); a = t
  END SUBROUTINE get_i2
  SUBROUTINE put_i3(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(IN) :: a(:,:,:)
    INTEGER, ALLOCATABLE :: t(:,:,:)
    TYPE(qeh5_dataset) :: d
    t = a
    CALL qeh5_set_space(d, t(1,1,1), 3, SHAPE(t))
    CALL qeh5_open_dataset(f, d, action='write', name=name)
    CALL qeh5_write_dataset(t, d)
    CALL qeh5_close(d)
  END SUBROUTINE put_i3
  SUBROUTINE get_i3(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: a(:,:,:)
    INTEGER, ALLOCATABLE :: t(:,:,:)
    TYPE(qeh5_dataset) :: d
    ALLOCATE(t(SIZE(a,1),SIZE(a,2),SIZE(a,3)))
    CALL qeh5_open_dataset(f, d, action='read', name=name)
    CALL qeh5_read_dataset(t, d)
    CALL qeh5_close(d); a = t
  END SUBROUTINE get_i3
  SUBROUTINE put_r1(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(dp), INTENT(IN) :: a(:)
    REAL(dp), ALLOCATABLE :: t(:)
    TYPE(qeh5_dataset) :: d
    t = a
    CALL qeh5_set_space(d, t(1), 1, [SIZE(t)])
    CALL qeh5_open_dataset(f, d, action='write', name=name)
    CALL qeh5_write_dataset(t, d)
    CALL qeh5_close(d)
  END SUBROUTINE put_r1
  SUBROUTINE put_r2(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(dp), INTENT(IN) :: a(:,:)
    REAL(dp), ALLOCATABLE :: t(:,:)
    TYPE(qeh5_dataset) :: d
    t = a
    CALL qeh5_set_space(d, t(1,1), 2, SHAPE(t))
    CALL qeh5_open_dataset(f, d, action='write', name=name)
    CALL qeh5_write_dataset(t, d)
    CALL qeh5_close(d)
  END SUBROUTINE put_r2
  SUBROUTINE put_r3(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(dp), INTENT(IN) :: a(:,:,:)
    REAL(dp), ALLOCATABLE :: t(:,:,:)
    TYPE(qeh5_dataset) :: d
    t = a
    CALL qeh5_set_space(d, t(1,1,1), 3, SHAPE(t))
    CALL qeh5_open_dataset(f, d, action='write', name=name)
    CALL qeh5_write_dataset(t, d)
    CALL qeh5_close(d)
  END SUBROUTINE put_r3
  SUBROUTINE put_c2(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    COMPLEX(dp), INTENT(IN) :: a(:,:)
    COMPLEX(dp), ALLOCATABLE :: t(:,:)
    TYPE(qeh5_dataset) :: d
    t = a
    CALL qeh5_set_space(d, t(1,1), 2, SHAPE(t))
    CALL qeh5_open_dataset(f, d, action='write', name=name)
    CALL qeh5_write_dataset(t, d)
    CALL qeh5_close(d)
  END SUBROUTINE put_c2
  SUBROUTINE get_i1(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER, INTENT(OUT) :: a(:)
    INTEGER, ALLOCATABLE :: t(:)
    TYPE(qeh5_dataset) :: d
    ALLOCATE(t(SIZE(a)))
    CALL qeh5_open_dataset(f, d, action='read', name=name)
    CALL qeh5_read_dataset(t, d)
    CALL qeh5_close(d); a = t
  END SUBROUTINE get_i1
  SUBROUTINE get_r1(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(dp), INTENT(OUT) :: a(:)
    REAL(dp), ALLOCATABLE :: t(:)
    TYPE(qeh5_dataset) :: d
    ALLOCATE(t(SIZE(a)))
    CALL qeh5_open_dataset(f, d, action='read', name=name)
    CALL qeh5_read_dataset(t, d)
    CALL qeh5_close(d); a = t
  END SUBROUTINE get_r1
  SUBROUTINE get_r2(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(dp), INTENT(OUT) :: a(:,:)
    REAL(dp), ALLOCATABLE :: t(:,:)
    TYPE(qeh5_dataset) :: d
    ALLOCATE(t(SIZE(a,1),SIZE(a,2)))
    CALL qeh5_open_dataset(f, d, action='read', name=name)
    CALL qeh5_read_dataset(t, d)
    CALL qeh5_close(d); a = t
  END SUBROUTINE get_r2
  SUBROUTINE get_r3(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(dp), INTENT(OUT) :: a(:,:,:)
    REAL(dp), ALLOCATABLE :: t(:,:,:)
    TYPE(qeh5_dataset) :: d
    ALLOCATE(t(SIZE(a,1),SIZE(a,2),SIZE(a,3)))
    CALL qeh5_open_dataset(f, d, action='read', name=name)
    CALL qeh5_read_dataset(t, d)
    CALL qeh5_close(d); a = t
  END SUBROUTINE get_r3
  SUBROUTINE get_c2(f, name, a)
    TYPE(qeh5_file), INTENT(IN) :: f
    CHARACTER(LEN=*), INTENT(IN) :: name
    COMPLEX(dp), INTENT(OUT) :: a(:,:)
    COMPLEX(dp), ALLOCATABLE :: t(:,:)
    TYPE(qeh5_dataset) :: d
    ALLOCATE(t(SIZE(a,1),SIZE(a,2)))
    CALL qeh5_open_dataset(f, d, action='read', name=name)
    CALL qeh5_read_dataset(t, d)
    CALL qeh5_close(d); a = t
  END SUBROUTINE get_c2
#endif
END MODULE qe_dump_h5
