!----------------------------------------------------------------------------
! Demonstrator for the LAXlib batched GPU generalized eigensolver.
!
! LAXlib/cdiaghg.f90 recently gained laxlib_cdiaghg_gpu_batched, which solves
! n_k independent generalized Hermitian eigenproblems H_k v_k = e_k S_k v_k
! with a single call into NVIDIA cuSOLVER's *Batched routines (cusolverDn
! ZpotrfBatched for the Cholesky factorization of S_k, cublasZtrsmBatched for
! the congruence transformation, cusolverDnZheevjBatched for the Jacobi
! eigensolver). The routine is reached through the "diaghg" generic
! interface (see LAXlib/laxlib_hi.h) once 3D device arrays plus an n_k
! argument are supplied, exactly as KS_Solvers/Davidson/cegterg.f90 does.
!
! Before any cuSOLVER/cuBLAS call, LAXlib expects the per-thread handle
! arrays declared in LAXlib/laxlib_cusolver_handle_mod.f90 to be allocated
! and bound to a CUDA stream:
!   CALL initialize_cusolver_handles(nthreads)   ! allocate cusolver_handle(:)
!   CALL initialize_cublas_handles(nthreads)     ! allocate cublas_handle(:)
!   CALL initialize_laxlib_cuda_stream(stream,i) ! bind thread i to a stream
! (cb_davidson_main.f90 does this once per OpenMP thread; here we only need
! a single serial "thread slot", so nthreads = 1.) The actual
! cusolverDnCreate/cublasCreate calls happen lazily, once, inside
! laxlib_cdiaghg_gpu_batched itself.
!
! The batched cuSOLVER call requires every matrix in the batch to share the
! same nominal dimension n (fixed leading dimension ldh too). Real
! Davidson k-points don't all have the same basis size, so smaller problems
! are PADDED up to n = max(real_dim(:)) exactly the way
! KS_Solvers/Davidson/cegterg.f90 pads its batched arrays: for row/col
! indices beyond the true dimension d_k,
!    H_pad(i,i) = large_value(i)   (>> any real eigenvalue, pushes the
!                                    padding modes to the top of the sorted
!                                    spectrum)
!    S_pad(i,i) = 1
!    every other padding entry, and every entry coupling the padding block
!    to the real d_k x d_k block, is zero.
! Because the padding block is block-diagonal with zero coupling to the
! real block, it does not perturb the true eigenpairs at all: the first
! d_k eigenvalues/eigenvectors coming out of the batched solve are exactly
! the eigenpairs of the true (H_k, S_k) problem, and the corresponding
! eigenvectors have (numerically) zero weight on the padded rows.
!
! S_k is NOT the identity for every batch entry: half of the batch uses a
! genuine (diagonally dominant, hence positive-definite) tridiagonal
! overlap matrix, to exercise the cusolverDnZpotrfBatched Cholesky step
! for real.
!
! This program builds n_k = 10 Hermitian generalized eigenproblems with
! distinct real dimensions ranging from 32x32 to 64x64, padded up to a common
! n = 64, solves all of them in a single batched call, and verifies the
! results on the host for every batch entry k by checking, over its true
! d_k x d_k block only:
!   (1) max_j || H_k v_j - e_j S_k v_j ||   (residual, generalized problem)
!   (2) || V_k^H S_k V_k - I ||             (S_k-orthonormality)
!   (3) max || v_j(padded rows) ||          (no leakage into the padding)
!----------------------------------------------------------------------------
#if defined(__CUDA)
PROGRAM demo_diaghg_gpu_batched
  !
  USE cudafor
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
  COMPLEX(DP), ALLOCATABLE :: h_h(:,:,:), s_h(:,:,:), v_h(:,:,:)
  REAL(DP),    ALLOCATABLE :: e_h(:,:)
  COMPLEX(DP), ALLOCATABLE, DEVICE :: h_d(:,:,:), s_d(:,:,:), v_d(:,:,:)
  REAL(DP),    ALLOCATABLE, DEVICE :: e_d(:,:)
  !
  INTEGER  :: i, j, l, p, k, d, me_bgrp, root_bgrp, intra_bgrp_comm
  REAL(DP) :: max_res, max_orth, max_leak, res_norm
  COMPLEX(DP) :: acc, hv, sv
  LOGICAL :: all_ok
  !
  me_bgrp = 0; root_bgrp = 0; intra_bgrp_comm = 0
  !
  ALLOCATE( h_h(ldh,n,n_k), s_h(ldh,n,n_k), v_h(ldh,m,n_k), e_h(n,n_k) )
  ALLOCATE( h_d(ldh,n,n_k), s_d(ldh,n,n_k), v_d(ldh,m,n_k), e_d(n,n_k) )
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
  h_d = h_h
  s_d = s_h
  !
  ! --- initialize the cuSOLVER / cuBLAS handles LAXlib needs (1 serial slot) ---
  CALL init_clocks(.true.)
  CALL initialize_cusolver_handles(1)
  CALL initialize_cublas_handles(1)
  CALL initialize_laxlib_cuda_stream( INT(0, KIND=cuda_stream_kind), 1 )
  !
  ! --- one batched call solves all n_k padded generalized eigenproblems ---
  CALL diaghg( n, m, h_d, s_d, ldh, e_d, v_d, n_k, me_bgrp, root_bgrp, intra_bgrp_comm )
  !
  e_h = e_d
  v_h = v_d
  !
  CALL finalize_cublas_handles()
  CALL finalize_cusolver_handles()
  !
  ! --- verify each batch entry, restricted to its true d x d block ---
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
           acc = (0.0_DP, 0.0_DP)
           DO l = 1, d
              sv = (0.0_DP, 0.0_DP)
              DO p = 1, d
                 sv = sv + s_h(l,p,k) * v_h(p,j,k)
              END DO
              acc = acc + CONJG(v_h(l,i,k)) * sv
           END DO
           IF (i == j) acc = acc - (1.0_DP, 0.0_DP)
           max_orth = MAX( max_orth, ABS(acc) )
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
     WRITE(*,*) 'PASSED: all ', n_k, ' padded batched eigenproblems verified.'
  ELSE
     WRITE(*,*) 'FAILED: see per-batch diagnostics above.'
     STOP 1
  END IF
  !
  DEALLOCATE(h_h, s_h, v_h, e_h, h_d, s_d, v_d, e_d)
  !
END PROGRAM demo_diaghg_gpu_batched
#else
PROGRAM demo_diaghg_gpu_batched
END PROGRAM demo_diaghg_gpu_batched
#endif
