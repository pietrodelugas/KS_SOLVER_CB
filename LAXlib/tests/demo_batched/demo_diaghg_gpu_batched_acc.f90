!----------------------------------------------------------------------------
! OpenACC-based demonstrator for the LAXlib batched GPU generalized
! eigensolver (laxlib_cdiaghg_gpu_batched, see also
! demo_diaghg_gpu_batched.f90 for the plain CUDA-Fortran version of this
! same test and a full description of the algorithm and the padding
! scheme).
!
! This version never "USE cudafor" and never declares a ", DEVICE" array.
! Instead it stages data on the GPU purely through OpenACC directives, the
! same way KS_Solvers/Davidson/cegterg.f90 and cbToy/cb_davidson_main.f90
! actually drive this routine in production:
!   * host arrays are mapped to the device with a structured
!     "!$acc data copyin(...) copyout(...)" region;
!   * "!$acc host_data use_device(...)" then hands the CUDA-Fortran
!     "diaghg" interface the device addresses of those OpenACC-mapped
!     arrays (this interop between OpenACC-managed memory and a
!     CUDA-Fortran DEVICE dummy argument is exactly what cegterg.f90 does
!     around its own batched "diaghg" call);
!   * the CUDA stream that LAXlib's cusolver/cublas handles must be bound
!     to is obtained from the OpenACC extension "acc_get_cuda_stream",
!     not from "cudafor" -- exactly as cb_davidson_main.f90 does via
!     mytime's clock_cuda_stream (itself an "acc_handle_kind" integer).
!
! Problem set and verification are identical to demo_diaghg_gpu_batched.f90:
! n_k = 10 Hermitian generalized eigenproblems H_k v = e S_k v with true
! dimensions ranging from 32x32 to 64x64 (zero-padded up to a common
! n = 64, decoupled padding block with large diagonal so it never mixes
! with the true spectrum), half of them using a genuine non-identity,
! diagonally-dominant (positive-definite) tridiagonal overlap matrix S_k.
! For every batch entry k, restricted to its true d_k x d_k block, we check:
!   (1) max_j || H_k v_j - e_j S_k v_j ||   (residual, generalized problem)
!   (2) || V_k^H S_k V_k - I ||             (S_k-orthonormality)
!   (3) max || v_j(padded rows) ||          (no leakage into the padding)
!----------------------------------------------------------------------------
#if defined(__CUDA)
PROGRAM demo_diaghg_gpu_batched_acc
  !
  USE openacc, ONLY : acc_get_cuda_stream, acc_async_sync, acc_handle_kind
  USE laxlib_cusolver_handles, ONLY : initialize_cusolver_handles, initialize_cublas_handles, &
                                       initialize_laxlib_cuda_stream, finalize_cusolver_handles, &
                                       finalize_cublas_handles
  IMPLICIT NONE
  include 'laxlib.fh'
  include 'laxlib_kinds.fh'
  !
  INTEGER, PARAMETER :: n_k = 10           ! number of matrices in the batch
  !! true dimension of each matrix in the batch: spans 32x32 .. 64x64
  INTEGER, PARAMETER :: real_dim(n_k) = (/ 32, 36, 40, 44, 48, 52, 56, 60, 62, 64 /)
  !! batch entries that get a non-identity (tridiagonal, PD) overlap matrix
  LOGICAL, PARAMETER :: nonidentity_s(n_k) = &
       (/ .false., .true., .false., .true., .false., .true., .false., .true., .false., .true. /)
  INTEGER, PARAMETER :: n   = 64           ! padded common dimension (= max(real_dim))
  INTEGER, PARAMETER :: ldh = n            ! leading dimension
  INTEGER, PARAMETER :: m   = n            ! number of eigenpairs requested
  REAL(DP), PARAMETER :: pad_scale = 1.0d3 ! padding eigenvalues >> any real one
  REAL(DP), PARAMETER :: tol = 1.0d-8
  !
  ! plain host arrays: no ", DEVICE" attribute anywhere in this program.
  ! They are mapped to the GPU below purely through OpenACC "data" clauses.
  COMPLEX(DP), ALLOCATABLE :: h_h(:,:,:), s_h(:,:,:), v_h(:,:,:)
  REAL(DP),    ALLOCATABLE :: e_h(:,:)
  !
  INTEGER  :: i, j, l, p, k, d, me_bgrp, root_bgrp, intra_bgrp_comm
  INTEGER(acc_handle_kind) :: my_cuda_stream
  REAL(DP) :: max_res, max_orth, max_leak, res_norm
  COMPLEX(DP) :: ovlp, hv, sv
  LOGICAL :: all_ok
  !
  me_bgrp = 0; root_bgrp = 0; intra_bgrp_comm = 0
  !
  ALLOCATE( h_h(ldh,n,n_k), s_h(ldh,n,n_k), v_h(ldh,m,n_k), e_h(n,n_k) )
  !
  ! --- build n_k distinct Hermitian generalized eigenproblems of true size ---
  ! --- real_dim(k), zero-padded up to the common batched dimension n       ---
  DO k = 1, n_k
     d = real_dim(k)
     h_h(:,:,k) = (0.0_DP, 0.0_DP)
     s_h(:,:,k) = (0.0_DP, 0.0_DP)
     !
     ! --- true d x d block ---
     DO i = 1, d
        h_h(i,i,k) = CMPLX( DBLE(i) + 0.1_DP*k, 0.0_DP, KIND=DP )
        DO j = i+1, d
           h_h(i,j,k) = CMPLX( DBLE(i+j)/(10.0_DP*k), DBLE(i-j)/(20.0_DP*k), KIND=DP )
           h_h(j,i,k) = CONJG( h_h(i,j,k) )
        END DO
     END DO
     !
     IF ( nonidentity_s(k) ) THEN
        ! diagonally-dominant (hence positive-definite) Hermitian tridiagonal
        ! overlap matrix: |S(i,i)| = 1 > sum of off-diagonal magnitudes (~0.07)
        DO i = 1, d
           s_h(i,i,k) = (1.0_DP, 0.0_DP)
        END DO
        DO i = 1, d-1
           s_h(i,i+1,k) = CMPLX( 0.05_DP, 0.02_DP, KIND=DP )
           s_h(i+1,i,k) = CONJG( s_h(i,i+1,k) )
        END DO
     ELSE
        DO i = 1, d
           s_h(i,i,k) = (1.0_DP, 0.0_DP)
        END DO
     END IF
     !
     ! --- padding block (d+1:n): decoupled from the real block, S_pad = I, ---
     ! --- H_pad diagonal >> any real eigenvalue so padding modes sort last ---
     DO i = d+1, n
        h_h(i,i,k) = CMPLX( pad_scale*DBLE(i), 0.0_DP, KIND=DP )
        s_h(i,i,k) = (1.0_DP, 0.0_DP)
     END DO
  END DO
  !
  ! --- initialize the cuSOLVER / cuBLAS handles LAXlib needs (1 serial slot); ---
  ! --- the CUDA stream comes from OpenACC, not from "cudafor"                ---
  CALL init_clocks(.true.)
  CALL initialize_cusolver_handles(1)
  CALL initialize_cublas_handles(1)
  my_cuda_stream = acc_get_cuda_stream( acc_async_sync )
  CALL initialize_laxlib_cuda_stream( my_cuda_stream, 1 )
  !
  ! --- map h_h, s_h to the device, solve the whole batch in one call, and ---
  ! --- bring e_h, v_h back -- all through OpenACC "data"/"host_data"      ---
  !$acc data copyin(h_h, s_h) copyout(e_h, v_h)
  !$acc host_data use_device(h_h, s_h, e_h, v_h)
  CALL diaghg( n, m, h_h, s_h, ldh, e_h, v_h, n_k, me_bgrp, root_bgrp, intra_bgrp_comm )
  !$acc end host_data
  !$acc end data
  !
  CALL finalize_cublas_handles()
  CALL finalize_cusolver_handles()
  !
  ! --- verify each batch entry, restricted to its true d x d block ---
  ! --- (h_h, s_h here still hold the original host-side matrices:  ---
  ! --- they were mapped "copyin" only, never copied back)          ---
  all_ok = .TRUE.
  DO k = 1, n_k
     d = real_dim(k)
     !
     ! residual of the generalized eigenproblem: H v_j - e_j S v_j, over the
     ! first d (true, non-padded) modes, using rows 1..d of H and S only --
     ! padding rows are decoupled by construction and checked separately below
     max_res = 0.0_DP
     DO j = 1, d
        DO i = 1, d
           hv = (0.0_DP, 0.0_DP)
           sv = (0.0_DP, 0.0_DP)
           DO l = 1, d
              hv = hv + h_h(i,l,k) * v_h(l,j,k)
              sv = sv + s_h(i,l,k) * v_h(l,j,k)
           END DO
           res_norm = ABS( hv - e_h(j,k)*sv )
           max_res = MAX( max_res, res_norm )
        END DO
     END DO
     !
     ! S-orthonormality of the true modes: V^H S V == I (restricted to 1..d)
     max_orth = 0.0_DP
     DO i = 1, d
        DO j = 1, d
           ovlp = (0.0_DP, 0.0_DP)
           DO l = 1, d
              sv = (0.0_DP, 0.0_DP)
              DO p = 1, d
                 sv = sv + s_h(l,p,k) * v_h(p,j,k)
              END DO
              ovlp = ovlp + CONJG(v_h(l,i,k)) * sv
           END DO
           IF (i == j) ovlp = ovlp - (1.0_DP, 0.0_DP)
           max_orth = MAX( max_orth, ABS(ovlp) )
        END DO
     END DO
     !
     ! leakage check: true modes must have (numerically) zero weight on the
     ! padded rows d+1:n, confirming the padding block does not perturb them
     max_leak = 0.0_DP
     DO j = 1, d
        DO i = d+1, n
           max_leak = MAX( max_leak, ABS(v_h(i,j,k)) )
        END DO
     END DO
     !
     WRITE(*,'(A,I3,A,I3,A,L1,A,ES10.3,A,ES10.3,A,ES10.3)') &
        'batch ', k, ':  dim = ', d, '  nonidentity_S = ', nonidentity_s(k), &
        '   max|Hv-eSv| = ', max_res, '   max|VHSV-I| = ', max_orth, &
        '   max|leak| = ', max_leak
     !
     IF ( max_res > tol .OR. max_orth > tol .OR. max_leak > tol ) all_ok = .FALSE.
  END DO
  !
  IF (all_ok) THEN
     WRITE(*,*) 'PASSED: all ', n_k, ' padded batched eigenproblems verified (OpenACC).'
  ELSE
     WRITE(*,*) 'FAILED: see per-batch diagnostics above.'
     STOP 1
  END IF
  !
  DEALLOCATE(h_h, s_h, v_h, e_h)
  !
END PROGRAM demo_diaghg_gpu_batched_acc
#else
PROGRAM demo_diaghg_gpu_batched_acc
END PROGRAM demo_diaghg_gpu_batched_acc
#endif
