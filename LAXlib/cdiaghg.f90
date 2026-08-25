!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE laxlib_cdiaghg( n, m, h, s, ldh, e, v, me_bgrp, root_bgrp, intra_bgrp_comm )
  !----------------------------------------------------------------------------
  !
  !! Called by diaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem.
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! complex matrices version.
  !! On output both matrix are unchanged.
  !!
  !! LAPACK version - uses both ZHEGV and ZHEGVX
  !!
  !
  USE laxlib_parallel_include
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: m
  !! number of eigenstates to be calculated
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,n)
  !! matrix to be diagonalized
  COMPLEX(DP), INTENT(INOUT) :: s(ldh,n)
  !! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
  !! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
  !! eigenvectors (column-wise)
  INTEGER,  INTENT(IN)  :: me_bgrp
  !! index of the processor within a band group
  INTEGER,  INTENT(IN)  :: root_bgrp
  !! index of the root processor within a band group
  INTEGER,  INTENT(IN)  :: intra_bgrp_comm
  !! intra band group communicator
  !
  INTEGER                  :: lwork, nb, mm, info, i, j
    ! mm = number of calculated eigenvectors
  REAL(DP)                 :: abstol
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
    ! various work space
  LOGICAL                  :: all_eigenvalues
 ! REAL(DP), EXTERNAL       :: DLAMCH
  INTEGER,  EXTERNAL       :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  !
  CALL start_clock( 'cdiaghg' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
     ALLOCATE( sdiag( n ) )
     DO i = 1, n
        sdiag(i) = DBLE( s(i,i) )
     END DO
     !
     all_eigenvalues = ( m == n )
     !
     ! ... check for optimal block size
     !
     nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
     !
     IF ( nb < 1 .OR. nb >= n) THEN
        !
        lwork = 2*n
        !
     ELSE
        !
        lwork = ( nb + 1 )*n
        !
     END IF
     !
     ALLOCATE( work( lwork ) )
     !
     IF ( all_eigenvalues ) THEN
        !
        ALLOCATE( rwork( 3*n - 2 ) )
        !
        ! ... calculate all eigenvalues (overwritten to v)
        !
        v(:,:) = h(:,:)
        !
        CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                    s, ldh, e, work, lwork, rwork, info )
        !
     ELSE
        !
        ALLOCATE( rwork( 7*n ) )
        !
        ! ... save the diagonal of input H (it will be overwritten)
        !
        ALLOCATE( hdiag( n ) )
        DO i = 1, n
           hdiag(i) = DBLE( h(i,i) )
        END DO
        !
        ALLOCATE( iwork( 5*n ) )
        ALLOCATE( ifail( n ) )
        !
        ! ... calculate only m lowest eigenvalues
        !
        abstol = 0.D0
       ! abstol = 2.D0*DLAMCH( 'S' )
        !
        ! ... the following commented lines calculate optimal lwork
        !
        !lwork = -1
        !
        !CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
        !             0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
        !             work, lwork, rwork, iwork, ifail, info )
        !
        !lwork = INT( work(1) ) + 1
        !
        !IF( lwork > SIZE( work ) ) THEN
        !   DEALLOCATE( work )
        !   ALLOCATE( work( lwork ) )
        !END IF
        !
        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
        !
        DEALLOCATE( ifail )
        DEALLOCATE( iwork )
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
        DO i = 1, n
           h(i,i) = CMPLX( hdiag(i), 0.0_DP ,kind=DP)
           DO j = i + 1, n
              h(i,j) = CONJG( h(j,i) )
           END DO
           DO j = n + 1, ldh
              h(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
        END DO
        !
        DEALLOCATE( hdiag )
        !
     END IF
     !
     !
     DEALLOCATE( rwork )
     DEALLOCATE( work )
     !
     IF ( info > n ) THEN
        CALL lax_error__( 'cdiaghg', 'S matrix not positive definite', ABS( info ) )
     ELSE IF ( info > 0 ) THEN
        CALL lax_error__( 'cdiaghg', 'eigenvectors failed to converge', ABS( info ) )
     ELSE IF ( info < 0 ) THEN
        CALL lax_error__( 'cdiaghg', 'incorrect call to ZHEGV*', ABS( info ) )
     END IF
     !
     ! ... restore input S matrix from saved diagonal and lower triangle
     !
     DO i = 1, n
        s(i,i) = CMPLX( sdiag(i), 0.0_DP ,kind=DP)
        DO j = i + 1, n
           s(i,j) = CONJG( s(j,i) )
        END DO
        DO j = n + 1, ldh
           s(j,i) = ( 0.0_DP, 0.0_DP )
        END DO
     END DO
     !
     DEALLOCATE( sdiag )
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array e', ABS( info ) )
  CALL MPI_BCAST( v, SIZE(v), MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array v', ABS( info ) )
#endif
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE laxlib_cdiaghg
!
!----------------------------------------------------------------------------
SUBROUTINE laxlib_cdiaghg_gpu( n, m, h_d, s_d, ldh, e_d, v_d, me_bgrp, root_bgrp, intra_bgrp_comm)
  !----------------------------------------------------------------------------
  !!
  !! Called by diaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! complex matrices version.
  !! On output both matrix are unchanged.
  !!
  !! GPU VERSION.
  !
#if defined(_OPENMP)
  USE omp_lib
#endif
  !
#if defined(__CUDA)
  USE cudafor
  !
  USE cusolverdn
  USE laxlib_cusolver_handles, ONLY : cusolver_handle, cusolver_initialized, laxlib_cuda_stream, & 
                                      cusolver_thread 
#endif
  !
  USE laxlib_parallel_include
  !
  ! NB: the flag below can be used to decouple LAXlib from devXlib.
  !     This will make devXlib an optional dependency of LAXlib when
  !     the library will be decoupled from QuantumESPRESSO.
#if defined(__USE_GLOBAL_BUFFER) && defined(__CUDA)
  USE device_fbuff_m,        ONLY : dev=>dev_buf, pin=>pin_buf
#define VARTYPE POINTER
#else
#define VARTYPE ALLOCATABLE
#endif
  !
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized
  INTEGER, INTENT(IN) :: m
  !! number of eigenstates to be calculated
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h_d(ldh,n)
  !! matrix to be diagonalized, allocated on the GPU
  COMPLEX(DP), INTENT(INOUT) :: s_d(ldh,n)
  !! overlap matrix, allocated on the GPU
  REAL(DP), INTENT(OUT) :: e_d(n)
  !! eigenvalues, , allocated on the GPU
  COMPLEX(DP),  INTENT(OUT) :: v_d(ldh,n)
  !! eigenvectors (column-wise), , allocated on the GPU
  INTEGER,  INTENT(IN)  :: me_bgrp
  !! index of the processor within a band group
  INTEGER,  INTENT(IN)  :: root_bgrp
  !! index of the root processor within a band group
  INTEGER,  INTENT(IN)  :: intra_bgrp_comm
  !! intra band group communicator
  !
#if defined(__CUDA)
    ATTRIBUTES(DEVICE) :: h_d, s_d, e_d, v_d
#endif
  !
  INTEGER              :: lwork, info
  !
  REAL(DP)             :: abstol
  INTEGER, ALLOCATABLE :: ifail(:)
  INTEGER, VARTYPE     :: iwork(:)
  REAL(DP), VARTYPE    :: rwork(:)
  COMPLEX(DP), VARTYPE :: work(:)
  REAL(DP), ALLOCATABLE :: e_orig_check(:) !!!M.IOvine - debugging line
  !
  COMPLEX(DP), VARTYPE :: v_h(:,:)
  REAL(DP), VARTYPE    :: e_h(:)
#if (! defined(__USE_GLOBAL_BUFFER)) && defined(__CUDA)
  ATTRIBUTES( PINNED ) :: work, iwork, rwork, v_h, e_h
#endif
  !
  INTEGER              :: lwork_d, lrwork_d, liwork, lrwork
  REAL(DP), VARTYPE    :: rwork_d(:)
  COMPLEX(DP), VARTYPE :: work_d(:)
  ! various work space
  !
  ! Temp arrays to save H and S.
  REAL(DP), VARTYPE    :: h_diag_d(:), s_diag_d(:)
#if defined(__CUDA)
  ATTRIBUTES( DEVICE ) :: work_d, rwork_d, h_diag_d, s_diag_d
  INTEGER                      :: devInfo_d, h_meig
  ATTRIBUTES( DEVICE )         :: devInfo_d
  TYPE(cusolverDnHandle), SAVE :: cuSolverHandle
  LOGICAL, SAVE                :: cuSolverInitialized = .FALSE.
  !
  COMPLEX(DP), VARTYPE   :: h_bkp_d(:,:), s_bkp_d(:,:)
  ATTRIBUTES( DEVICE )   :: h_bkp_d, s_bkp_d
#endif
  INTEGER :: i, j
#undef VARTYPE
  !
  !
  !
  !
  CALL start_clock_gpu( 'cdiaghg' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
      !
      ! Keeping compatibility for both CUSolver and CustomEigensolver, CUSolver below
      !
#if defined(__CUDA)

#if ! defined(__USE_GLOBAL_BUFFER)
      ALLOCATE(h_bkp_d(n,n), s_bkp_d(n,n), STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate h_bkp_d or s_bkp_d ', ABS( info ) )
#else
      CALL dev%lock_buffer( h_bkp_d,  (/ n, n /), info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate h_bkp_d ', ABS( info ) )
      CALL dev%lock_buffer( s_bkp_d,  (/ n, n /), info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate s_bkp_d ', ABS( info ) )
#endif
      !
!$cuf kernel do(2) <<<*,*,0,laxlib_cuda_stream>>>
      DO j=1,n
         DO i=1,n
            h_bkp_d(i,j) = h_d(i,j)
            s_bkp_d(i,j) = s_d(i,j)
         ENDDO
      ENDDO
      !
#if defined(_OPENMP)
      IF (omp_get_num_threads() > 1) CALL lax_error__( ' cdiaghg_gpu ', 'cdiaghg_gpu is not thread-safe',  ABS( info ) )
#endif
      IF ( .NOT. cusolver_initialized(cusolver_thread) ) THEN
         info = cusolverDnCreate(cusolver_handle(cusolver_thread))
         IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnCreate',  ABS( info ) )
         cusolver_initialized(cusolver_thread) = .TRUE.
         info = cusolverDnSetStream(cusolver_handle(cusolver_thread), laxlib_cuda_stream )
         IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnSetStream',  ABS( info ) )   
      ENDIF
      !
      cuSolverHandle = cusolver_handle(cusolver_thread)
      info = cusolverDnZhegvdx_bufferSize(cuSolverHandle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                                               n, h_d, ldh, s_d, ldh, 0.D0, 0.D0, 1, m, h_meig, e_d, lwork_d)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', ' cusolverDnZhegvdx_bufferSize failed ', ABS( info ) )
      !
#if ! defined(__USE_GLOBAL_BUFFER)
      ALLOCATE(work_d(1*lwork_d), STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate work_d ', ABS( info ) )
#else
      CALL dev%lock_buffer( work_d,  lwork_d, info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate work_d ', ABS( info ) )
#endif
      !
      info = cusolverDnZhegvdx(cuSolverHandle, CUSOLVER_EIG_TYPE_1, CUSOLVER_EIG_MODE_VECTOR, CUSOLVER_EIG_RANGE_I, CUBLAS_FILL_MODE_UPPER, &
                                  n, h_d, ldh, s_d, ldh, 0.D0, 0.D0, 1, m, h_meig, e_d, work_d, lwork, devInfo_d)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', ' cusolverDnZhegvdx failed ', ABS( info ) )
!$cuf kernel do(2) <<<*,*,0,laxlib_cuda_stream>>>
      DO j=1,n
         DO i=1,n
            IF(j <= m) v_d(i,j) = h_d(i,j)
            h_d(i,j) = h_bkp_d(i,j)
            s_d(i,j) = s_bkp_d(i,j)
         ENDDO
      ENDDO
      !
      !
      ! Do not destroy the handle to save the (re)creation time on each call.
      !
      !info = cusolverDnDestroy(cuSolverHandle)
      !IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', ' cusolverDnDestroy failed ', ABS( info ) )
      !
#if ! defined(__USE_GLOBAL_BUFFER)
      DEALLOCATE(work_d)
      DEALLOCATE(h_bkp_d, s_bkp_d)
#else
      CALL dev%release_buffer( work_d,  info )
      CALL dev%release_buffer( h_bkp_d, info )
      CALL dev%release_buffer( s_bkp_d, info )
#endif
      !
      ! Keeping compatibility for both CUSolver and CustomEigensolver, CustomEigensolver below
      !
#else
     CALL lax_error__( 'cdiaghg', 'Called GPU eigensolver without GPU support', 1 )
#endif
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
#if defined __GPU_MPI
  info = cudaDeviceSynchronize()
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error synchronizing device (first)', ABS( info ) )
  CALL MPI_BCAST( e_d, n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array e_d', ABS( info ) )
  CALL MPI_BCAST( v_d, ldh*m, MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array v_d', ABS( info ) )
  info = cudaDeviceSynchronize() ! this is probably redundant...
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error synchronizing device (second)', ABS( info ) )
#else
  ALLOCATE(e_h(n), v_h(ldh,m))
  e_h(1:n) = e_d(1:n)
  v_h(1:ldh, 1:m) = v_d(1:ldh, 1:m)
  CALL MPI_BCAST( e_h, n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array e_d', ABS( info ) )
  CALL MPI_BCAST( v_h, ldh*m, MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error broadcasting array v_d', ABS( info ) )
  e_d(1:n) = e_h(1:n)
  v_d(1:ldh, 1:m) = v_h(1:ldh, 1:m)
  DEALLOCATE(e_h, v_h)
#endif
#endif
  !
  !!M.IOvine - debug line:
  IF (.NOT. ALLOCATED(e_orig_check)) ALLOCATE(e_orig_check(n))
   e_orig_check = e_d(1:n)    
  print *, 'ORIGINAL cdiaghg_gpu e_d =', e_orig_check(1:5)
  !
  CALL stop_clock_gpu( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE laxlib_cdiaghg_gpu
!----------------------------------------------------------------------------
!!! M.Iovine - created new subroutine for batched diagonalization:
SUBROUTINE laxlib_cdiaghg_gpu_batched( n, m, h_d, s_d, ldh, e_d, v_d, n_k, me_bgrp, root_bgrp, intra_bgrp_comm)
  !----------------------------------------------------------------------------
  !!
  !! Called by diaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! complex matrices version.
  !! On output both matrix are unchanged.
  !!
  !! GPU VERSION.
  !
#if defined(_OPENMP)
  USE omp_lib
#endif
  !
#if defined(__CUDA)
  USE cudafor
  USE ieee_arithmetic, ONLY: ieee_is_nan, ieee_is_finite !!! JUST FOR DEBUGGING!
  USE cublas !! M.Iovine - we need to explicitly include the use clublas
  !
  USE openacc, only: c_devptr !! M.Iovine - we include c_devptr
  USE cusolverdn
  USE laxlib_cusolver_handles, ONLY : cusolver_handle, cusolver_initialized, laxlib_cuda_stream, & 
                                      cusolver_thread 
#endif
  !
  USE laxlib_parallel_include
  !
  ! NB: the flag below can be used to decouple LAXlib from devXlib.
  !     This will make devXlib an optional dependency of LAXlib when
  !     the library will be decoupled from QuantumESPRESSO.
#if defined(__USE_GLOBAL_BUFFER) && defined(__CUDA)
  USE device_fbuff_m,        ONLY : dev=>dev_buf, pin=>pin_buf
#define VARTYPE POINTER
#else
#define VARTYPE ALLOCATABLE
#endif
  !
  IMPLICIT NONE
  include 'laxlib_kinds.fh'
  !
  !!!! M.Iovine - nbase is set equal to the maximum value among the 
  !!!! k-pints (threads)
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized
  !!!! number of desired root eigenstates for each k point is the same
  !!!! for each k-point!
  INTEGER, INTENT(IN) :: m
  !! number of eigenstates to be calculated
  !!!! M.Iovine - Batched change:
  INTEGER, INTENT(IN) :: n_k  
  !!!!
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h_d(ldh,n,n_k) !! M.Iovine - added a dimension -It is still a static array
  !! matrix to be diagonalized, allocated on the GPU
  COMPLEX(DP), INTENT(INOUT) :: s_d(ldh,n,n_k) !! M.Iovine - added a dimension -It is still a static array
  !! overlap matrix, allocated on the GPU
  !!!! M.Iovine - Batched change: the array was made 3D
  REAL(DP), INTENT(OUT) :: e_d(n,n_k)
  !!!!
  !! eigenvalues, , allocated on the GPU
  COMPLEX(DP),  INTENT(OUT) :: v_d(ldh,n,n_k) !!!! M.Iovine - v_d becomes a 3D array! It is still a static array
  !! eigenvectors (column-wise), , allocated on the GPU
  INTEGER,  INTENT(IN)  :: me_bgrp
  !! index of the processor within a band group
  INTEGER,  INTENT(IN)  :: root_bgrp
  !! index of the root processor within a band group
  INTEGER,  INTENT(IN)  :: intra_bgrp_comm
  !! intra band group communicator
  COMPLEX(DP), ALLOCATABLE :: h_sym_chk(:,:) !!M.Iovine - lines added for debugging
  REAL(DP) :: max_asym !!M.Iovine - lines added for debugging
  INTEGER :: ii_sym, jj_sym !!M.Iovine - lines added for debugging
  INTEGER, ALLOCATABLE :: dinfo_host(:)!! Added for debugging
  !
#if defined(__CUDA)
    ATTRIBUTES(DEVICE) :: h_d, s_d, e_d, v_d
#endif
  !
!!!! DEBUGGING LINES!1
#if defined(__CUDA)
    COMPLEX(DP), ALLOCATABLE :: nan_chk(:,:)
    LOGICAL, ALLOCATABLE :: nan_mask_r(:,:), nan_mask_i(:,:)
    LOGICAL :: has_nan_dbg
    REAL(DP), ALLOCATABLE :: nan_echk(:)
    LOGICAL, ALLOCATABLE :: nan_emask(:)
    LOGICAL, ALLOCATABLE :: inf_mask_r(:,:), inf_mask_i(:,:)
    LOGICAL :: has_inf_dbg
#endif
  
  
  INTEGER              :: lwork, info
  !
  REAL(DP)             :: abstol
  INTEGER, ALLOCATABLE :: ifail(:)
  INTEGER, VARTYPE     :: iwork(:)
  REAL(DP), VARTYPE    :: rwork(:)
  COMPLEX(DP), VARTYPE :: work(:)
  !
  COMPLEX(DP), VARTYPE :: v_h(:,:,:) !!!! M.Iovine : added a dimension to the array
  REAL(DP), VARTYPE    :: e_h(:,:)
#if (! defined(__USE_GLOBAL_BUFFER)) && defined(__CUDA)
  ATTRIBUTES( PINNED ) :: work, iwork, rwork, v_h, e_h
#endif
  !
  INTEGER              :: lwork_d, lrwork_d, liwork, lrwork
  REAL(DP), VARTYPE    :: rwork_d(:)
  COMPLEX(DP), VARTYPE :: work_d(:)
  !!!! M.Iovine - We add a line to define the info array corresponding to the Batched routine of the Cusolver NVIDIA library:
  INTEGER, ALLOCATABLE :: d_info(:)
  !!!! M.Iovine - We define host arrays of c_devpr pointers:
  TYPE(c_devptr), ALLOCATABLE :: arr_of_ptr_s(:)
  TYPE(c_devptr), ALLOCATABLE :: arr_of_ptr_h(:)
  !!!! M.Iovine - We define device arrays of c_devpr pointers:
  TYPE(c_devptr), ALLOCATABLE, DEVICE :: arr_of_ptr_h_d(:), arr_of_ptr_s_d(:)
  COMPLEX(DP) :: alpha=(1.D0, 0.D0) !! M.Iovine - we define the alpha coefficient for the triangular cublas linear solver 
  INTEGER :: b, ik,ind_min !! M.Iovine - indices for eigenvectors ordering 
  COMPLEX(DP) :: minim !! M.Iovine - min. value for padding
  COMPLEX(DP) :: min_temp !! M.Iovine - temp variable for swapping of arrays for ordering eigenvalues and eigenvect
  ! various work space
  !
  ! Temp arrays to save H and S.
  REAL(DP), VARTYPE    :: h_diag_d(:), s_diag_d(:)
#if defined(__CUDA)
  ATTRIBUTES( DEVICE ) :: work_d, rwork_d, h_diag_d, s_diag_d, d_info, arr_of_ptr_s_d, arr_of_ptr_h_d !!!! M.Iovine - added d_info and arr_of_ptr to the device attributes
  INTEGER                      :: devInfo_d, h_meig
  ATTRIBUTES( DEVICE )         :: devInfo_d
  !TYPE(cusolverDnHandle), SAVE :: cuSolverHandle !M.Iovine - commented line
  TYPE(cusolverDnHandle), SAVE :: cuSolverHandle_batched !! M.Iovine - we define a proper handle for the batched cdiaghg_gpu
  LOGICAL, SAVE                :: cuSolverInitialized = .FALSE.
  !
  !! Device arrays to save the old values of the Hamiltonian and Overlap matrices!!
  COMPLEX(DP), VARTYPE   :: h_bkp_d(:,:,:), s_bkp_d(:,:,:) !!!! M.Iovine : The array to store the old values of the Hamiltonian before appying the 
  ATTRIBUTES( DEVICE )   :: h_bkp_d, s_bkp_d

  !!!! M.Iovine - We instantiate a cusolverDnSyevjInfo variable that is needed for the Batched routine provided by NVIDIA:
  TYPE(cusolverDnSyevjInfo) :: syevj_params

  !! M.Iovine - we define the variables for the cublas handle:
  INTEGER(kind=cuda_stream_kind) :: mycudaStream
  type(cublasHandle) :: cublas_handle
  INTEGER :: istat_cublas

#endif
  INTEGER :: i, j, k !!!! M.IOvine - added index k for the third dimension of the arrays
  mycudaStream = laxlib_cuda_stream
  istat_cublas = cublasCreate(cublas_handle) !!! M.Iovine - introduced cublas handle
  istat_cublas = cublasSetStream(cublas_handle, mycudaStream) !!!M.IOvine - COMMENTED TO DEBUGG!!
  IF (istat_cublas /= 0) CALL lax_error__( ' cdiaghg_gpu ', 'cublasSetStream', ABS(istat_cublas) ) 
#undef VARTYPE

  !
  !
  !
  !
  CALL start_clock_gpu( 'cdiaghg' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
      print *, 'DEBUG: n =', n, ' m =', m, ' ldh =', ldh, ' n_k =', n_k
      !
      ! Keeping compatibility for both CUSolver and CustomEigensolver, CUSolver below
      !
#if defined(__CUDA)

#if ! defined(__USE_GLOBAL_BUFFER)
      ALLOCATE(h_bkp_d(n,n,n_k), s_bkp_d(n,n,n_k), STAT = info) !!!! M.Iovine : h_bkp_d is a 3 dimensional array now!!
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate h_bkp_d or s_bkp_d ', ABS( info ) )
      !! M.Iovine - We allocate the array d_info necessary for the cusolver diagonalization:
      ALLOCATE(d_info(n_k),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate d_info ', ABS( info ) )
      !! M.Iovine - we allocate the host arrays of pointers for the Cholesky factorization:
      ALLOCATE(arr_of_ptr_s(n_k),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate arr_of_ptr_s ', ABS( info ) )
      ALLOCATE(arr_of_ptr_h(n_k),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate arr_of_ptr_h ', ABS( info ) )
      !! M.Iovine - we allocate the device arrays of pointers for the Cholesky factorization:
      ALLOCATE(arr_of_ptr_s_d(n_k),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate arr_of_ptr_s ', ABS( info ) )
      ALLOCATE(arr_of_ptr_h_d(n_k),STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate arr_of_ptr_h ', ABS( info ) )

#else
      CALL dev%lock_buffer( h_bkp_d,  (/ n, n /), info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate h_bkp_d ', ABS( info ) )
      CALL dev%lock_buffer( s_bkp_d,  (/ n, n /), info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate s_bkp_d ', ABS( info ) )
#endif
      !
!$cuf kernel do(3) <<<*,*,0,laxlib_cuda_stream>>> !!!! M.Iovine - Modified the loops nested from 2 to 3 --> so we changed kernel do(2) to kenel do(3) :
      DO k=1,n_k
         DO j=1,n
            DO i=1,n
                h_bkp_d(i,j,k) = h_d(i,j,k)
                s_bkp_d(i,j,k) = s_d(i,j,k)
            ENDDO
         ENDDO
      ENDDO
      !i
!!!!DEBUG LINES:
IF (.NOT. ALLOCATED(nan_chk)) ALLOCATE(nan_chk(n,n), nan_mask_r(n,n), nan_mask_i(n,n), &
                                        inf_mask_r(n,n), inf_mask_i(n,n))
nan_chk = h_d(1:n,1:n,1)
nan_mask_r = ieee_is_nan(REAL(nan_chk))
nan_mask_i = ieee_is_nan(AIMAG(nan_chk))
inf_mask_r = (.NOT. ieee_is_finite(REAL(nan_chk))) .AND. (.NOT. nan_mask_r)
inf_mask_i = (.NOT. ieee_is_finite(AIMAG(nan_chk))) .AND. (.NOT. nan_mask_i)
has_nan_dbg = ANY(nan_mask_r) .OR. ANY(nan_mask_i)
has_inf_dbg = ANY(inf_mask_r) .OR. ANY(inf_mask_i)
print *, '[1] h_d input has_nan =', has_nan_dbg, ' has_inf =', has_inf_dbg, &
         ' maxabs =', MAXVAL(ABS(nan_chk), MASK = ieee_is_finite(REAL(nan_chk)) .AND. ieee_is_finite(AIMAG(nan_chk)))
!!!!

#if defined(_OPENMP)
      IF (omp_get_num_threads() > 1) CALL lax_error__( ' cdiaghg_gpu ', 'cdiaghg_gpu is not thread-safe',  ABS( info ) )
#endif
      !IF ( .NOT. cusolver_initialized(cusolver_thread) ) THEN
       ! info = cusolverDnCreate(cusolver_handle(cusolver_thread))
        !IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnCreate',  ABS( info ) )
        !cusolver_initialized(cusolver_thread) = .TRUE.
        !info = cusolverDnSetStream(cusolver_handle(cusolver_thread), laxlib_cuda_stream )
        !IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnSetStream',  ABS( info ) )   
      !ENDIF
      IF( .NOT. cuSolverInitialized ) THEN   !!! M.Iovine - we introduce a new Handle to avoid leaving changes to the next calls done
                                                    !through the kernel calls inside the iterative loop!!
        info = cusolverDnCreate(cuSolverHandle_batched)
        IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu_batched ', 'cusolverDnCreate',  ABS( info ) )
         cuSolverInitialized = .TRUE.
         info = cusolverDnSetStream(cuSolverHandle_batched, laxlib_cuda_stream) 
         IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu_batched ', 'cusolverDnSetStream',  ABS( info ) )
      ENDIF 


    !!!! M.Iovine - Cholesky factorization of the s overlap matrix:
    !! We assign to each element of arr_of_ptr the device pointer to each 2D slice in s_d :
    do k = 1, n_k
        arr_of_ptr_s(k) = c_devloc(s_d(:,:,k))
    end do
    
    !!! M.Iovine - we copy the c_devptr in the host array arr_of_ptr_s to the device array arr_of_ptr_s_d
    arr_of_ptr_s_d = arr_of_ptr_s
    
    
    !!!DEBUGG :
    istat_cublas = cudaGetLastError()
IF (istat_cublas /= 0) THEN
   print *, 'STICKY CUDA ERROR before ZpotrfBatched: ', &
             cudaGetErrorString(istat_cublas)
END IF
    !!!!

     
    !cuSolverHandle = cusolver_handle(cusolver_thread) !!M.Iovine - this line must before any cuSolver routine kernel call!
    info = cusolverDnZpotrfBatched(cuSolverHandle_batched, CUBLAS_FILL_MODE_LOWER, n, arr_of_ptr_s_d, ldh, d_info(1), n_k)
    IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnZpotrfBatched',  ABS( info ) )
    !!!!
    
    !!!M.Iovine - Cholesky check:
    info = cudaDeviceSynchronize()
    IF (info /= 0) CALL lax_error__(' cdiaghg_gpu ', 'sync after Cholesky', ABS(info))
    
    IF (.NOT. ALLOCATED(dinfo_host)) ALLOCATE(dinfo_host(n_k))
    dinfo_host = d_info(1:n_k)
    print *, '[CHOLESKY d_info] per-batch status =', dinfo_host
    !!!!!

!!! DEBUGGING LINES:
IF (.NOT. ALLOCATED(nan_chk)) ALLOCATE(nan_chk(n,n), nan_mask_r(n,n), nan_mask_i(n,n), &
                                        inf_mask_r(n,n), inf_mask_i(n,n))
nan_chk = h_d(1:n,1:n,1)
nan_mask_r = ieee_is_nan(REAL(nan_chk))
nan_mask_i = ieee_is_nan(AIMAG(nan_chk))
inf_mask_r = (.NOT. ieee_is_finite(REAL(nan_chk))) .AND. (.NOT. nan_mask_r)
inf_mask_i = (.NOT. ieee_is_finite(AIMAG(nan_chk))) .AND. (.NOT. nan_mask_i)
has_nan_dbg = ANY(nan_mask_r) .OR. ANY(nan_mask_i)
has_inf_dbg = ANY(inf_mask_r) .OR. ANY(inf_mask_i)
print *, '[2] h_d input has_nan =', has_nan_dbg, ' has_inf =', has_inf_dbg, &
         ' maxabs =', MAXVAL(ABS(nan_chk), MASK = ieee_is_finite(REAL(nan_chk)) .AND. ieee_is_finite(AIMAG(nan_chk)))
!!!!

    !! We assign to each element of arr_of_ptr the device pointer to each 2D slice in h_d :
    do k = 1, n_k
        arr_of_ptr_h(k) = c_devloc(h_d(:,:,k))
    end do
    
    !!! M.Iovine - we copy the c_devptr in the host array arr_of_ptr_h to the device array arr_of_ptr_h_d
    arr_of_ptr_h_d = arr_of_ptr_h


    !!!!M.Iovine - triagular pre and post multiplication of Hamiltonian:
    !! Ly = H :
    info = cublasZtrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, &
           CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, &
           n, n, alpha, arr_of_ptr_s_d, ldh, arr_of_ptr_h_d, ldh, n_k)
    IF ( info /= CUBLAS_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cublasZtrsmBatched-LEFT',  ABS( info ) )

    info = cudaDeviceSynchronize() !!M.Iovine - added a synchronize
    IF (info /= 0) CALL lax_error__(' cdiaghg_gpu ', 'sync after triang left', ABS(info))

!!! DEBUGGING LINES
IF (.NOT. ALLOCATED(nan_chk)) ALLOCATE(nan_chk(n,n), nan_mask_r(n,n), nan_mask_i(n,n), &
                                        inf_mask_r(n,n), inf_mask_i(n,n))
nan_chk = h_d(1:n,1:n,1)
nan_mask_r = ieee_is_nan(REAL(nan_chk))
nan_mask_i = ieee_is_nan(AIMAG(nan_chk))
inf_mask_r = (.NOT. ieee_is_finite(REAL(nan_chk))) .AND. (.NOT. nan_mask_r)
inf_mask_i = (.NOT. ieee_is_finite(AIMAG(nan_chk))) .AND. (.NOT. nan_mask_i)
has_nan_dbg = ANY(nan_mask_r) .OR. ANY(nan_mask_i)
has_inf_dbg = ANY(inf_mask_r) .OR. ANY(inf_mask_i)
print *, '[3] h_d input has_nan =', has_nan_dbg, ' has_inf =', has_inf_dbg, &
         ' maxabs =', MAXVAL(ABS(nan_chk), MASK = ieee_is_finite(REAL(nan_chk)) .AND. ieee_is_finite(AIMAG(nan_chk)))
!!!!

    !! y = w L(conj. transpose) :
    info = cublasZtrsmBatched(cublas_handle, CUBLAS_SIDE_RIGHT, &
           CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_C, CUBLAS_DIAG_NON_UNIT, &
           n, n, alpha, arr_of_ptr_s_d, ldh, arr_of_ptr_h_d, ldh, n_k)    
    IF ( info /= CUBLAS_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cublasZtrsmBatched-RIGT',  ABS( info ) )
      !
    
    info = cudaDeviceSynchronize() !!M.Iovine - added a synchronize
    IF (info /= 0) CALL lax_error__(' cdiaghg_gpu ', 'sync after triang right', ABS(info))

!!! DEBUGGING LINES:
IF (.NOT. ALLOCATED(nan_chk)) ALLOCATE(nan_chk(n,n), nan_mask_r(n,n), nan_mask_i(n,n), &
                                        inf_mask_r(n,n), inf_mask_i(n,n))
nan_chk = h_d(1:n,1:n,1)
nan_mask_r = ieee_is_nan(REAL(nan_chk))
nan_mask_i = ieee_is_nan(AIMAG(nan_chk))
inf_mask_r = (.NOT. ieee_is_finite(REAL(nan_chk))) .AND. (.NOT. nan_mask_r)
inf_mask_i = (.NOT. ieee_is_finite(AIMAG(nan_chk))) .AND. (.NOT. nan_mask_i)
has_nan_dbg = ANY(nan_mask_r) .OR. ANY(nan_mask_i)
has_inf_dbg = ANY(inf_mask_r) .OR. ANY(inf_mask_i)
print *, '[4] h_d input has_nan =', has_nan_dbg, ' has_inf =', has_inf_dbg, &
         ' maxabs =', MAXVAL(ABS(nan_chk), MASK = ieee_is_finite(REAL(nan_chk)) .AND. ieee_is_finite(AIMAG(nan_chk)))


      !!!! M.Iovine - we define the parameters for the Jacobi algorithm corresponding to the diagonalization done through the batched cuSolver routine:
      info = cusolverDnCreateSyevjInfo(syevj_params)
      IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnCreateSyevjInfo',  ABS( info ) )
      info = cusolverDnXsyevjSetTolerance(syevj_params, 0.D0) !! The tolerance is set to the default value (0)
      IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnXsyevjSetTolerance',  ABS( info ) )
      info = cusolverDnXsyevjSetMaxSweeps(syevj_params, 100) !!! The maximum swe    eps is the default value (100)
      IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnXsyevjSetMaxSweeps',  ABS( info ) )
      info = cusolverDnXsyevjSetSortEig(syevj_params, 1) !!M.Iovine - ordering in ascending order
      IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnXsyevjSetSortEig',  ABS( info ) )
      !!!!

      !!!! M.Iovine - We change the routine from the single kernel call to the batched routine of NVIDIA Cusolver:
      info = cusolverDnZheevjBatched_bufferSize(cuSolverHandle_batched, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, &
                                               n, h_d, ldh, e_d, lwork_d, syevj_params, n_k)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', ' cusolverDnZheevjBatched failed ', ABS( info ) )
      !
#if ! defined(__USE_GLOBAL_BUFFER)
      ALLOCATE(work_d(1*lwork_d), STAT = info)
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate work_d ', ABS( info ) )
      !!!! M.Iovine - We allocate the array d_info and check if it's successful :
      !ALLOCATE(d_info(n_k), STAT=info)
      !IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate d_info ', ABS( info ) )
#else
      CALL dev%lock_buffer( work_d,  lwork_d, info )
      IF( info /= 0 ) CALL lax_error__( ' cdiaghg_gpu ', ' cannot allocate work_d ', ABS( info ) )
#endif
      !
      !! Debugging lines:
      IF (.NOT. ALLOCATED(h_sym_chk)) ALLOCATE(h_sym_chk(n,n))
        h_sym_chk = h_d(1:n,1:n,1)
        max_asym = 0.0_DP
        DO jj_sym = 1, n
        DO ii_sym = 1, jj_sym-1
                max_asym = MAX(max_asym, ABS(h_sym_chk(ii_sym,jj_sym) - CONJG(h_sym_chk(jj_sym,ii_sym))))
          END DO
        END DO
        print *, '[SYM CHECK] max Hermitian asymmetry in h_d before eigensolver =', max_asym


      !!!! M.Iovine - We change the routine from the single kernel call to the batched routine of NVIDIA Cusolver:
      info = cusolverDnZheevjBatched(cuSolverHandle_batched, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, &
      n, h_d, ldh, e_d, work_d, lwork_d, d_info(1), syevj_params, n_k)
      IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', ' cusolverDnZheevjBatched failed ', ABS( info ) )
     
     info = cudaDeviceSynchronize() !!M.Iovine - added a synchronize
    IF (info /= 0) CALL lax_error__(' cdiaghg_gpu ', 'sync after batched diag.', ABS(info))
    dinfo_host = d_info(1:n_k)
    print *, '[DIAG CUSOLVER CHECK NEW d_info] per-batch status =', dinfo_host

     !!! DEBUGGING LINES: 
IF (.NOT. ALLOCATED(nan_chk)) ALLOCATE(nan_chk(n,n), nan_mask_r(n,n), nan_mask_i(n,n), &
                                        inf_mask_r(n,n), inf_mask_i(n,n))
nan_chk = h_d(1:n,1:n,1)
nan_mask_r = ieee_is_nan(REAL(nan_chk))
nan_mask_i = ieee_is_nan(AIMAG(nan_chk))
inf_mask_r = (.NOT. ieee_is_finite(REAL(nan_chk))) .AND. (.NOT. nan_mask_r)
inf_mask_i = (.NOT. ieee_is_finite(AIMAG(nan_chk))) .AND. (.NOT. nan_mask_i)
has_nan_dbg = ANY(nan_mask_r) .OR. ANY(nan_mask_i)
has_inf_dbg = ANY(inf_mask_r) .OR. ANY(inf_mask_i)
print *, '[5] h_d input has_nan =', has_nan_dbg, ' has_inf =', has_inf_dbg, &
         ' maxabs =', MAXVAL(ABS(nan_chk), MASK = ieee_is_finite(REAL(nan_chk)) .AND. ieee_is_finite(AIMAG(nan_chk)))


      !! M.Iovine - we destroy the SyevjInfo object:
      info = cusolverDnDestroySyevjInfo(syevj_params)
    !!! LAST PART ADDED FOR DEBUGGING:
    do k = 1, n_k
        arr_of_ptr_h(k) = c_devloc(h_d(:,:,k))
    end do

    do k = 1, n_k
        arr_of_ptr_s(k) = c_devloc(s_d(:,:,k))
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!! M.Iovine : we need to take into account the fact that the eigenvalues got after the factorization and the triangular matrix mult. are the same of the initial problem, but this is not true for the eigenvectors, so we need to solve a triangular system to get the effective eigenvectors:
    !! M.Iovine : it is important to observe that the current eigenvectors are stored in the h_d matrix: h_d = L(conj. transpose) * eigvect
    info = cublasZtrsmBatched(cublas_handle, CUBLAS_SIDE_LEFT, &
           CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_C, CUBLAS_DIAG_NON_UNIT, &
           n, n, alpha, arr_of_ptr_s_d, ldh, arr_of_ptr_h_d, ldh, n_k)
    IF ( info /= CUBLAS_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cublasZtrsmBatched-eigenvectors',  ABS( info ) )
   
   info = cudaDeviceSynchronize() !!M.Iovine - added a synchronize
    IF (info /= 0) CALL lax_error__(' cdiaghg_gpu ', 'sync after triang final', ABS(info))

   !!! DEBUGGING LINES:   
IF (.NOT. ALLOCATED(nan_chk)) ALLOCATE(nan_chk(n,n), nan_mask_r(n,n), nan_mask_i(n,n), &
                                        inf_mask_r(n,n), inf_mask_i(n,n))
nan_chk = h_d(1:n,1:n,1)
nan_mask_r = ieee_is_nan(REAL(nan_chk))
nan_mask_i = ieee_is_nan(AIMAG(nan_chk))
inf_mask_r = (.NOT. ieee_is_finite(REAL(nan_chk))) .AND. (.NOT. nan_mask_r)
inf_mask_i = (.NOT. ieee_is_finite(AIMAG(nan_chk))) .AND. (.NOT. nan_mask_i)
has_nan_dbg = ANY(nan_mask_r) .OR. ANY(nan_mask_i)
has_inf_dbg = ANY(inf_mask_r) .OR. ANY(inf_mask_i)
print *, '[6] h_d input has_nan =', has_nan_dbg, ' has_inf =', has_inf_dbg, &
         ' maxabs =', MAXVAL(ABS(nan_chk), MASK = ieee_is_finite(REAL(nan_chk)) .AND. ieee_is_finite(AIMAG(nan_chk)))

IF (.NOT. ALLOCATED(nan_echk)) ALLOCATE(nan_echk(n), nan_emask(n))
info = cudaDeviceSynchronize()
IF (info /= 0) CALL lax_error__('cdiaghg_gpu', &
                                 'sync before reading e_d', ABS(info))
nan_echk = e_d(1:n,1)
nan_emask = ieee_is_nan(nan_echk)
print *, '[E] e_d has_nan =', ANY(nan_emask), ' first vals =', nan_echk(1:5)

    !!!! M.Iovine - We need to order in acending way the eigenvalues
    !!!! and the corresponding eigenvectors:
    !! we add !$cuf kernel do(1) because we need to do the reordering
    !! on the device --> 
    !!$cuf kernel do(1) <<<*,*,0,laxlib_cuda_stream>>>
    !do ik = 1, n_k
     !   do b = 1, n-1
      !      ind_min = b
       !     do k = b,n
        !        if (e_d(k, ik) .le. e_d(ind_min, ik)) then
         !           ind_min = k
          !      end if
           ! end do
           ! minim = e_d(ind_min,ik)
           ! e_d(ind_min,ik) = e_d(b,ik)
           ! e_d(b,ik) = minim
            !!! We swap also the columns with indices equal to the
            !!! eigenvalues ones in order to reorder also the 
            !!! eigenvectors array!
            !do i = 1, n
             !   min_temp = h_d(i,ind_min,ik)
              !  h_d(i,ind_min,ik) = h_d(i,b,ik)
               ! h_d(i,b,ik) = min_temp
            !end do
        !end do
    !end do



!!!! M.Iovine - Modified the loops nested from 2 to 3 --> so we changed kernel do(2) to kernel do(3) :
!!$cuf kernel do(3) <<<*,*,0,laxlib_cuda_stream>>>
      !DO k=1,n_k
       !  DO j=1,n
        !    DO i=1,n
         !       IF(j <= m) v_d(i,j,k) = h_d(i,j,k) !!!!M.Iovine - the array becomes 3D and also m is an array because 
          !      h_d(i,j,k) = h_bkp_d(i,j,k)
           !     s_d(i,j,k) = s_bkp_d(i,j,k)
            !ENDDO
        ! ENDDO
      !ENDDO
!$cuf kernel do(3) <<<*,*,0,laxlib_cuda_stream>>>
DO k = 1, n_k
   DO j = 1, m
      DO i = 1, n
         v_d(i, j, k) = h_d(i, j, k)
      END DO
   END DO
END DO

!$cuf kernel do(3) <<<*,*,0,laxlib_cuda_stream>>>
DO k = 1, n_k
   DO j = 1, n
      DO i = 1, n
         h_d(i, j, k) = h_bkp_d(i, j, k)
         s_d(i, j, k) = s_bkp_d(i, j, k)
      END DO
   END DO
END DO


info = cudaDeviceSynchronize()
 IF (info /= 0) CALL lax_error__('cdiaghg_gpu', &
                                 'sync before reading e_d', ABS(info))

!
      !
      ! Do not destroy the handle to save the (re)creation time on each call.
      !
      !info = cusolverDnDestroy(cuSolverHandle_batched)
      !IF( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu_batched ', ' cusolverDnDestroy_batched failed ', ABS( info ) )
      !

 !!! DEBUGGING LINES:   
IF (.NOT. ALLOCATED(nan_chk)) ALLOCATE(nan_chk(n,n), nan_mask_r(n,n), nan_mask_i(n,n), &
                                        inf_mask_r(n,n), inf_mask_i(n,n))
nan_chk = v_d(1:n,1:n,1)
nan_mask_r = ieee_is_nan(REAL(nan_chk))
nan_mask_i = ieee_is_nan(AIMAG(nan_chk))
inf_mask_r = (.NOT. ieee_is_finite(REAL(nan_chk))) .AND. (.NOT. nan_mask_r)
inf_mask_i = (.NOT. ieee_is_finite(AIMAG(nan_chk))) .AND. (.NOT. nan_mask_i)
has_nan_dbg = ANY(nan_mask_r) .OR. ANY(nan_mask_i)
has_inf_dbg = ANY(inf_mask_r) .OR. ANY(inf_mask_i)
print *, '[7] v_d input has_nan =', has_nan_dbg, ' has_inf =', has_inf_dbg, &
         ' maxabs =', MAXVAL(ABS(nan_chk), MASK = ieee_is_finite(REAL(nan_chk)) .AND. ieee_is_finite(AIMAG(nan_chk)))
      
      
      
#if ! defined(__USE_GLOBAL_BUFFER)
      DEALLOCATE(work_d)
      DEALLOCATE(h_bkp_d, s_bkp_d)
      DEALLOCATE(d_info) !!!! M.Iovine - we deallocate d_info
      DEALLOCATE(arr_of_ptr_h) !! M.Iovine - we deallocate arr_of_ptr
      DEALLOCATE(arr_of_ptr_s) !! M.Iovine - we deallocate arr_of_ptr
      DEALLOCATE(arr_of_ptr_h_d) !! M.Iovine - we deallocate arr_of_ptr
      DEALLOCATE(arr_of_ptr_s_d) !! M.Iovine - we deallocate arr_of_ptr
#else
      CALL dev%release_buffer( work_d,  info )
      CALL dev%release_buffer( h_bkp_d, info )
      CALL dev%release_buffer( s_bkp_d, info )
#endif
!      IF (ALLOCATED(h_sym_chk)) DEALLOCATE(h_sym_chk) !!M.IOvine - debugging lines
!! We destroy the cublas handle:
!#if defined(__CUDA)
 !   istat_cublas = cublasDestroy(cublas_handle)
!#endif

      !
      ! Keeping compatibility for both CUSolver and CustomEigensolver, CustomEigensolver below
      !
#else
     CALL lax_error__( 'cdiaghg', 'Called GPU eigensolver without GPU support', 1 )
#endif

!! We destroy the cublas handle:
#if defined(__CUDA)
    istat_cublas = cublasDestroy(cublas_handle)
#endif

info = cudaDeviceSynchronize()
 IF (info /= 0) CALL lax_error__('cdiaghg_gpu', &
                                 'sync before reading e_d', ABS(info))

     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#if defined __MPI
#if defined __GPU_MPI
  info = cudaDeviceSynchronize()
  IF ( info /= 0 ) &
        CALL lax_error__( 'cdiaghg', 'error synchronizing device (first)', ABS( info ) )
  !!!! M.Iovine - We add a loop to take into account the k-points of the batch - we need to make it be executed on the device :
  !$cuf kernel do(3) <<<*,*,0,laxlib_cuda_stream>>>
  DO k=1,n_k 
    CALL MPI_BCAST( e_d(:,k), n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
    IF ( info /= 0 ) &
            CALL lax_error__( 'cdiaghg', 'error broadcasting array e_d', ABS( info ) )
    CALL MPI_BCAST( v_d(:,:,k), ldh*m, MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
    IF ( info /= 0 ) &
            CALL lax_error__( 'cdiaghg', 'error broadcasting array v_d', ABS( info ) )
    info = cudaDeviceSynchronize() ! this is probably redundant...
    IF ( info /= 0 ) &
            CALL lax_error__( 'cdiaghg', 'error synchronizing device (second)', ABS( info ) )
  END DO  
#else
  !!!! M.Iovine - We add a loop to take into account the k-points of th batch :
  ALLOCATE(e_h(n, n_k), v_h(ldh, m, n_k))
  DO k=1,n_k
    e_h(1:n, k) = e_d(1:n, k)
    v_h(1:ldh, 1:m, k) = v_d(1:ldh, 1:m, k)
    CALL MPI_BCAST( e_h(:,k), n, MPI_DOUBLE_PRECISION, root_bgrp, intra_bgrp_comm, info )
    IF ( info /= 0 ) &
            CALL lax_error__( 'cdiaghg', 'error broadcasting array e_d', ABS( info ) )
    CALL MPI_BCAST( v_h(:,:,k), ldh*m, MPI_DOUBLE_COMPLEX, root_bgrp, intra_bgrp_comm, info )
    IF ( info /= 0 ) &
            CALL lax_error__( 'cdiaghg', 'error broadcasting array v_d', ABS( info ) )
    e_d(1:n, k) = e_h(1:n, k)
    v_d(1:ldh, 1:m, k) = v_h(1:ldh, 1:m, k)
  END DO
  DEALLOCATE(e_h, v_h)
#endif
#endif
  !
  CALL stop_clock_gpu( 'cdiaghg' )
  !
  
  !!!DEBUGGING LINES:
  IF (.NOT. ALLOCATED(nan_echk)) ALLOCATE(nan_echk(n), nan_emask(n))
nan_echk = e_d(1:n,1)
nan_emask = ieee_is_nan(nan_echk)
print *, '[E2] e_d has_nan =', ANY(nan_emask), ' first vals =', nan_echk(1:5)

!!!! DEBUGGING LINES:  
#if defined(__CUDA)
        IF (ALLOCATED(nan_chk)) DEALLOCATE(nan_chk, nan_mask_r, nan_mask_i, inf_mask_r, inf_mask_i)
        IF (ALLOCATED(nan_echk)) DEALLOCATE(nan_echk, nan_emask)
#endif
  RETURN
  !
END SUBROUTINE laxlib_cdiaghg_gpu_batched

!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE laxlib_pcdiaghg( n, h, s, ldh, e, v, idesc )
  !----------------------------------------------------------------------------
  !
  !! Called by pdiaghg interface.
  !! Calculates eigenvalues and eigenvectors of the generalized problem.
  !! Solve Hv = eSv, with H symmetric matrix, S overlap matrix.
  !! complex matrices version.
  !! On output both matrix are unchanged.
  !!
  !! Parallel version with full data distribution
  !!
  !
  USE laxlib_parallel_include
  USE laxlib_descriptor,      ONLY : la_descriptor, laxlib_intarray_to_desc
  USE laxlib_processors_grid, ONLY : ortho_parent_comm
#if defined __SCALAPACK
  USE laxlib_processors_grid, ONLY : ortho_cntx, np_ortho, me_ortho, ortho_comm
  USE zhpev_module,     ONLY : pzheevd_drv
#endif
  !
  IMPLICIT NONE
  !
  include 'laxlib_kinds.fh'
  include 'laxlib_param.fh'
  include 'laxlib_mid.fh'
  include 'laxlib_low.fh'
  !
  INTEGER, INTENT(IN) :: n
  !! dimension of the matrix to be diagonalized and number of eigenstates to be calculated
  INTEGER, INTENT(IN) :: ldh
  !! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,ldh)
  !! matrix to be diagonalized
  COMPLEX(DP), INTENT(INOUT) :: s(ldh,ldh)
  !! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
  !! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,ldh)
  !! eigenvectors (column-wise)
  INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
  !! laxlib descriptor
  !  
  TYPE(la_descriptor) :: desc
  !
  INTEGER, PARAMETER  :: root = 0
  INTEGER             :: nx, info
#if defined __SCALAPACK
  INTEGER             :: descsca( 16 )
#endif
    ! local block size
  COMPLEX(DP), ALLOCATABLE :: ss(:,:), hh(:,:), tt(:,:)
    ! work space used only in parallel diagonalization
  !
  ! ... input s and h are copied so that they are not destroyed
  !
  CALL start_clock( 'cdiaghg' )
  !
  CALL laxlib_intarray_to_desc(desc,idesc)
  !
  IF( desc%active_node > 0 ) THEN
     !
     nx   = desc%nrcx
     !
     IF( nx /= ldh ) &
        CALL lax_error__(" pcdiaghg ", " inconsistent leading dimension ", ldh )
     !
     ALLOCATE( hh( nx, nx ) )
     ALLOCATE( ss( nx, nx ) )
     !
     hh(1:nx,1:nx) = h(1:nx,1:nx)
     ss(1:nx,1:nx) = s(1:nx,1:nx)
     !
  END IF

  CALL start_clock( 'cdiaghg:choldc' )
  !
  ! ... Cholesky decomposition of sl ( L is stored in sl )
  !
  IF( desc%active_node > 0 ) THEN
     !
#if defined __SCALAPACK
     CALL descinit( descsca, n, n, desc%nrcx, desc%nrcx, 0, 0, ortho_cntx, SIZE( ss, 1 ) , info )
     !
     IF( info /= 0 ) CALL lax_error__( ' cdiaghg ', ' desckinit ', ABS( info ) )
#endif
     !
#if defined __SCALAPACK

     CALL pzpotrf( 'L', n, ss, 1, 1, descsca, info )

     IF( info /= 0 ) CALL lax_error__( ' cdiaghg ', ' problems computing cholesky ', ABS( info ) )
#else
     CALL laxlib_pzpotrf( ss, nx, n, idesc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:choldc' )
  !
  ! ... L is inverted ( sl = L^-1 )
  !
  CALL start_clock( 'cdiaghg:inversion' )
  !
  IF( desc%active_node > 0 ) THEN
     !
#if defined __SCALAPACK
     !CALL clear_upper_tr( ss )
     ! set to zero the upper triangle of ss
     !
     CALL sqr_setmat( 'U', n, ZERO, ss, size(ss,1), idesc )
     !
     CALL pztrtri( 'L', 'N', n, ss, 1, 1, descsca, info )
     !
     IF( info /= 0 ) CALL lax_error__( ' cdiaghg ', ' problems computing inverse ', ABS( info ) )
#else
     CALL laxlib_pztrtri( ss, nx, n, idesc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:inversion' )
  !
  ! ... vl = L^-1*H
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, idesc )
     !
  END IF
  !
  ! ... hl = ( L^-1*H )*(L^-1)^T
  !
  IF( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'N', 'C', n, ONE, v, nx, ss, nx, ZERO, hh, nx, idesc )
     !
     ! ensure that "hh" is really Hermitian, it is sufficient to set the diagonal
     ! properly, because only the lower triangle of hh will be used
     ! 
     CALL sqr_setmat( 'H', n, ZERO, hh, size(hh,1), idesc )
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  !
  IF ( desc%active_node > 0 ) THEN
     ! 
#ifdef TEST_DIAG
     CALL test_drv_begin()
#endif

#if defined(__SCALAPACK)
     !
     CALL pzheevd_drv( .true., n, desc%nrcx, hh, e, ortho_cntx, ortho_comm )
     !
#else
     !
     CALL laxlib_pzheevd( .true., n, idesc, hh, SIZE( hh, 1 ), e )
     !
#endif
     !
#ifdef TEST_DIAG
     CALL test_drv_end()
#endif
     !
  END IF
  !
  ! ... v = (L^T)^-1 v
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc%active_node > 0 ) THEN
     !
     CALL sqr_mm_cannon( 'C', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, idesc )
     !
  END IF
  !
#if defined __MPI
  CALL MPI_BCAST( e, SIZE(e), MPI_DOUBLE_PRECISION, root, ortho_parent_comm, info )
  IF ( info /= 0 ) &
        CALL lax_error__( 'pcdiaghg', 'error broadcasting array e', ABS( info ) )
#endif
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc%active_node > 0 ) THEN
     DEALLOCATE( ss, hh )
  END IF
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE test_drv_begin()
     ALLOCATE( tt( n, n ) )
     CALL laxlib_zsqmcll_x( n, hh, nx, tt, n, desc, desc%comm )
     RETURN
  END SUBROUTINE test_drv_begin
  !
  SUBROUTINE test_drv_end()
     !
     INTEGER :: i, j, k
     COMPLEX(DP), ALLOCATABLE :: diag(:,:)
     !
     IF( desc%myc == 0 .AND. desc%myr == 0 ) THEN

        write( 100, fmt="(A20,2D18.10)" ) ' e code = ', e( 1 ), e( n )
        ALLOCATE( diag( n*(n+1)/2, 1 ) )
        k = 1
        !write( 100, fmt="(I5)" ) n
        DO j = 1, n
           DO i = j, n
              diag( k, 1 ) = tt( i, j )
              !write( 100, fmt="(2I5,2D18.10)" ) i, j, tt( i, j )
              k = k + 1
           END DO
        END DO
        call zhpev_drv( 'V', 'L', N, diag(:,1), e, tt, n )
        write( 100, fmt="(A20,2D18.10)" ) ' e test = ', e( 1 ), e( n )
        !write( 100, * ) 'eigenvalues and eigenvectors'
        DO j = 1, n
           !write( 100, fmt="(1I5,1D18.10,A)" ) j, e( j )
           DO i = 1, n
              !write( 100, fmt="(2I5,2D18.10)" ) i, j, tt( i, j )
           END DO
        END DO
        close(100)
        DEALLOCATE( diag )
     END IF
#if defined __MPI
     CALL MPI_BCAST( tt, SIZE(tt), MPI_DOUBLE_COMPLEX, 0, desc%comm, info )
     IF ( info /= 0 ) &
        CALL lax_error__( 'test_drv_end', 'error broadcasting array e', ABS( info ) )
#endif
     CALL laxlib_zsqmdst_x( n, tt, n, hh, nx, desc )
     DEALLOCATE( tt )
     CALL lax_error__('cdiaghg','stop serial',1)
     RETURN
  END SUBROUTINE test_drv_end
  !
END SUBROUTINE laxlib_pcdiaghg
!
