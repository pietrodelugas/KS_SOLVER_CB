!
! Copyright (C) 2002-2026 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Host scalars (alpha, beta) for the BLAS-like library calls made in "host pointer mode" on the device
! (cuBLAS, rocBLAS, oneMKL, ...). Plain Fortran, no device runtime needed.
!
! The scalar of a complex(DP) BLAS call is passed by address. The literal constants that the compiler creates
! for an expression like ( 1.D0, 0.D0 ) are only aligned to 8 bytes, but cuBLASLt (CUDA 12.x) reads them with
! aligned 16 byte SSE loads (movdqa) in some code paths (e.g. n=1 products) and crashes with a SIGSEGV on x86_64.
! The variables of this module are allocated on the heap (16 byte aligned), and the alignment is checked.
!
! usage:  USE device_blas_scalars_m, ONLY : blas_cone, blas_czero, ...
!         CALL device_blas_scalars_init()      ! once, from any thread (idempotent)
!         the arguments are used as the alpha and beta of the calls
!
MODULE device_blas_scalars_m
  !
  USE util_param, ONLY : DP
  USE iso_c_binding, ONLY : c_loc, c_intptr_t, c_ptr
  IMPLICIT NONE
  PRIVATE
  !
  COMPLEX(DP), ALLOCATABLE, TARGET, PUBLIC :: blas_cone, blas_czero, blas_cmone
  REAL(DP),    ALLOCATABLE, TARGET, PUBLIC :: blas_done, blas_dzero, blas_dmone
  LOGICAL, SAVE :: scalars_ready = .FALSE.
  INTEGER, PARAMETER :: required_alignment = 16
  !
  PUBLIC :: device_blas_scalars_init
  !
  CONTAINS
  !
  SUBROUTINE device_blas_scalars_init()
    IMPLICIT NONE
    IF ( scalars_ready ) RETURN
    !$omp critical (device_blas_scalars_init_cs)
    IF ( .NOT. scalars_ready ) THEN
       ALLOCATE( blas_cone, blas_czero, blas_cmone, blas_done, blas_dzero, blas_dmone )
       blas_cone  = ( 1.0_DP, 0.0_DP )
       blas_czero = ( 0.0_DP, 0.0_DP )
       blas_cmone = ( -1.0_DP, 0.0_DP )
       blas_done  = 1.0_DP
       blas_dzero = 0.0_DP
       blas_dmone = -1.0_DP
       CALL check_alignment( c_loc(blas_cone),  'blas_cone' )
       CALL check_alignment( c_loc(blas_czero), 'blas_czero' )
       CALL check_alignment( c_loc(blas_cmone), 'blas_cmone' )
       scalars_ready = .TRUE.
    END IF
    !$omp end critical (device_blas_scalars_init_cs)
  END SUBROUTINE device_blas_scalars_init
  !
  SUBROUTINE check_alignment( p, name )
    USE iso_c_binding, ONLY : c_ptr
    IMPLICIT NONE
    TYPE(c_ptr), INTENT(IN) :: p
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(c_intptr_t) :: addr
    addr = TRANSFER( p, addr )
    IF ( MOD( addr, INT(required_alignment, c_intptr_t) ) /= 0 ) &
       CALL errore( 'device_blas_scalars_init', TRIM(name)//' is not aligned to 16 bytes', 1 )
  END SUBROUTINE check_alignment
  !
END MODULE device_blas_scalars_m
