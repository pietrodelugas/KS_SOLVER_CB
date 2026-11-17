!:/
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! NOTE (Ivan Carnimeo, May, 05th, 2022): 
!   cegterg and regterg have been ported to GPU with OpenACC, 
!   the previous CUF versions (cegterg_gpu and regterg_gpu) have been removed, 
!   and now cegterg and regterg are used for both CPU and GPU execution.
!   If you want to see the previous code checkout to commit: df3080b231c5daf52295c23501fbcaa9bfc4bfcc (on Thu Apr 21 06:18:02 2022 +0000)
!
! Modified for batched k-point processing
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!

MODULE cegterg_mod
  IMPLICIT NONE
  CONTAINS
!----------------------------------------------------------------------------
SUBROUTINE cegterg( h_psi_ptr, s_psi_ptr,  uspp, g_psi_ptr, h_psi_ptr_iter, g_psi_ptr_iter, &
                    npw, npwx, nvec, nvecx, npol, evc, ethr, e, btype, notcnv, lrot, dav_iter, nhpsi, n_k )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
  ! ... The n_k k-points ("slots" or "batches") of a group are solved together by ONE host thread and ONE OpenACC
  ! ... queue (clock_thread, normally 1): every step of the Davidson algorithm is executed for all the k-points that
  ! ... still have to converge, one after the other or, where a library supports it, as a single batched call.
  ! ... All the work arrays carry the k-point (batch) index as last dimension.
  !
#if defined(__CUDA)
  use cublas
  use cudafor
#endif
  USE nvtx
  USE device_blas_scalars_m, ONLY : blas_cone, blas_czero, device_blas_scalars_init
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
                            nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp,            ONLY : mp_size
  USE mytime,          ONLY : clock_thread, clock_cuda_stream
  USE openacc,         ONLY : acc_get_cuda_stream, c_devptr, acc_deviceptr
  USE laxlib_cusolver_handles, ONLY :  laxlib_cuda_stream,cusolver_thread, cublas_handle, cublas_initialized
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: n_k
    ! number of k-points (batches) solved together
  INTEGER, INTENT(IN) :: npw(n_k)
    ! number of plane waves of every k-point
  INTEGER, INTENT(IN) :: npwx, nvec, nvecx, npol
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! number of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx*npol,nvec,n_k)
    !  evc contains the  refined estimates of the eigenvectors (on the device, present)
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty ), present on the device
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec,n_k)
    ! contains the estimated roots (device pointer, see the call)
  INTEGER, INTENT(OUT) :: dav_iter(n_k), notcnv(n_k)
    ! number of iterations performed, number of unconverged roots, per k-point
  INTEGER, INTENT(OUT) :: nhpsi(n_k)
    ! total number of individual hpsi, per k-point
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, np, kdmx, n, m, ipol, ib
  INTEGER :: nbase(n_k), nb1(n_k), my_n(n_k), kdim(n_k)
    ! dimension of the reduced basis, first new vector, number of basis vectors and size of the vectors, per k-point
  INTEGER :: n_start, n_end
  INTEGER :: ierr
  COMPLEX(DP), ALLOCATABLE :: hc(:,:,:), sc(:,:,:), vc(:,:,:), vc_tmp(:,:,:)
  !$acc declare device_resident(hc, sc, vc, vc_tmp)
    ! Hamiltonian, S matrix and eigenvectors on the reduced basis, per k-point
  LOGICAL :: restart_now(n_k)
  INTEGER, ALLOCATABLE :: src_idx(:,:) ! unconverged roots, in order, for the reordering of vc (host and device copy)
  INTEGER, ALLOCATABLE :: par(:,:)
    ! per k-point integer parameters, host and device copy, updated with one transfer: the kernels that run over all the
    ! k-points of the group in one launch read them (1 nbase, 2 nb1, 3 notcnv, 4 npw, 5 active = not converged yet)
  INTEGER, PARAMETER :: ip_nbase = 1, ip_nb1 = 2, ip_notcnv = 3, ip_npw = 4, ip_act = 5
  INTEGER :: ibb, jj, kk, ii, ntail, jmax
  INTEGER :: nb_slot_max, col
  REAL(DP) :: thr_conv
  REAL(DP) :: nrm
  INTEGER :: notcnv_slot(n_k), n_slot
  INTEGER :: nvecx_pad ! nvecx + nvec: padding so that batched Zgemm calls sized on the max over batches stay in bounds
  REAL(DP), ALLOCATABLE :: ew(:,:)
    ! eigenvalues of the reduced hamiltonian, per k-point
  COMPLEX(DP), ALLOCATABLE :: psi(:,:,:), hpsi(:,:,:), spsi(:,:,:)
    ! work space, contains psi, H*psi, S*psi of every k-point
  LOGICAL :: conv(nvec,n_k)
    ! true if the root is converged
  REAL(DP) :: empty_ethr
    ! threshold for empty bands
  ! work arrays shared by the k-points
  COMPLEX(DP), ALLOCATABLE :: psi_fft_comp(:,:,:), hpsi_fft_comp(:,:,:)
  INTEGER, ALLOCATABLE :: npw_comp(:), nb1_comp(:), slot2batch(:), conv_comp(:,:)
  REAL(DP), ALLOCATABLE :: eig_new_comp(:,:), ew_old_comp(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc_comp(:,:,:), sc_comp(:,:,:), vc_comp(:,:,:)
  REAL(DP), ALLOCATABLE :: ew_comp(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc_comp_itr(:,:,:), sc_comp_itr(:,:,:), vc_comp_itr(:,:,:)
  REAL(DP), ALLOCATABLE :: ew_comp_itr(:,:)
  TYPE(c_devptr), ALLOCATABLE :: ptr_hc(:), ptr_vc(:), ptr_sc(:), ptr_psi(:), ptr_hpsi(:), ptr_spsi(:), ptr_psi_result(:), &
                                 ptr_hpsi_nb1(:)
  INTEGER :: nbase_comp(n_k), notcnv_comp(n_k), my_n_comp(n_k)
  LOGICAL :: done_comp(n_k), my_done(n_k)
  INTEGER :: my_slot, n_active, nbase_max, prev_nbase_max, prev_n_active
  INTEGER :: n_active_zgem, notcnv_max, my_n_max, kdim_max
  INTEGER :: info
  INTEGER :: i, j, k
  INTEGER :: async_id
  INTEGER(kind=cuda_stream_kind) :: mycudaStream
  !
  EXTERNAL  h_psi_ptr,    s_psi_ptr,    g_psi_ptr,  h_psi_ptr_iter, g_psi_ptr_iter
    ! h_psi_ptr(npwx,npw,nvec,psi,hpsi,i_batch,n_k)           H|psi> of all the k-points (first call)
    ! s_psi_ptr(npwx,npw,nvec,psi,spsi,i_batch)              S|psi> of one k-point (if needed)
    ! h_psi_ptr_iter, g_psi_ptr_iter: H|psi> and g_psi of the correction vectors of the active slots
  !
  nhpsi = 0
  CALL start_clock( 'cegterg' )
  !
  ! one host thread, one queue: clock_thread is the index of the queue and of the cuBLAS/cuSOLVER handles
  async_id = clock_thread
  IF ( mp_size(intra_bgrp_comm) > 1 ) CALL errore( 'cegterg', &
       'batched g_psi / normalization not implemented with plane-wave parallelism', 1 )
  IF ( mp_size(inter_bgrp_comm) > 1 ) CALL errore( 'cegterg', &
       'batched cegterg is not implemented with band-group parallelism', 1 )
  !
#if defined(__CUDA)
  ! Get cuda stream and link cublas to it
  mycudaStream = clock_cuda_stream
  IF ( .NOT. cublas_initialized(async_id) ) THEN
     info = cublasCreate(cublas_handle(async_id))
     IF ( info /= CUBLAS_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cublasCreate',  ABS( info ) )
     cublas_initialized(async_id) = .TRUE.
     info = cublasSetStream(cublas_handle(async_id), mycudaStream )
     IF ( info /= CUBLAS_STATUS_SUCCESS ) CALL lax_error__( ' cdiaghg_gpu ', 'cusolverDnSetStream',  ABS( info ) )
  ENDIF
#endif
  !$acc data deviceptr(e) async(async_id)
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     kdim(:) = npw(:)
     kdmx = npwx
  ELSE
     kdim(:) = npwx*npol
     kdmx = npwx*npol
  END IF
  !
  nvecx_pad = nvecx + nvec
  ALLOCATE(  psi( npwx*npol, nvecx_pad, n_k ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' cegterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx*npol, nvecx_pad, n_k ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' cegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !$acc enter data async(async_id) create(psi, hpsi)
  ! zero the (padded) work space: the batched Zgemm calls also read the padding columns
  !$acc kernels async(async_id) present(psi, hpsi)
  psi = ZERO
  hpsi = ZERO
  !$acc end kernels
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx*npol, nvecx_pad, n_k ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' cegterg ',' cannot allocate spsi ', ABS(ierr) )
     !$acc enter data async(async_id) create(spsi)
     !$acc kernels async(async_id) present(spsi)
     spsi = ZERO
     !$acc end kernels
  END IF
  !
  ALLOCATE( sc( nvecx_pad, nvecx_pad, n_k ), hc( nvecx_pad, nvecx_pad, n_k ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' cegterg ',' cannot allocate hc, sc ', ABS(ierr) )
  ALLOCATE( vc( nvecx, nvecx, n_k ), vc_tmp( nvecx, nvec, n_k ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' cegterg ',' cannot allocate vc ', ABS(ierr) )
  ALLOCATE( ew( nvecx, n_k ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' cegterg ',' cannot allocate ew ', ABS(ierr) )
  !$acc enter data create(ew) async(async_id)
  ALLOCATE( par(5,n_k), src_idx(nvec,n_k) )
  par = 0
  src_idx = 0
  !$acc enter data create(par, src_idx) async(async_id)
  !
  ! work arrays shared by the k-points (device only)
  ALLOCATE( psi_fft_comp(npwx*npol,nvecx,n_k), hpsi_fft_comp(npwx*npol,nvecx,n_k) )
  ALLOCATE( npw_comp(n_k), nb1_comp(n_k), slot2batch(n_k) )
  ALLOCATE( conv_comp(nvec,n_k), ew_old_comp(nvec,n_k), eig_new_comp(nvec,n_k) )
  ALLOCATE( hc_comp(nvec,nvec,n_k), sc_comp(nvec,nvec,n_k), vc_comp(nvec,nvec,n_k), ew_comp(nvec,n_k) )
  ALLOCATE( ptr_psi(n_k), ptr_psi_result(n_k), ptr_hpsi(n_k), ptr_hpsi_nb1(n_k), ptr_spsi(n_k), ptr_hc(n_k), ptr_sc(n_k), ptr_vc(n_k) )
  !$acc enter data async(async_id) create(psi_fft_comp, hpsi_fft_comp, npw_comp, nb1_comp, slot2batch, conv_comp, ew_old_comp, &
  !$acc eig_new_comp, hc_comp, sc_comp, vc_comp, ew_comp, ptr_psi, ptr_psi_result, ptr_spsi, ptr_hpsi, ptr_hpsi_nb1, ptr_hc, ptr_sc, ptr_vc)
  prev_nbase_max = -1
  prev_n_active  = -1
  !
  notcnv(:) = nvec
  nbase(:)  = nvec
  conv(:,:) = .FALSE.
  my_done(:) = .FALSE.
  done_comp(:) = .FALSE.
  !
  ! psi = evc
  !$acc parallel loop collapse(3) async(async_id) present(psi, evc)
  DO ib = 1, n_k
     DO j = 1, nvec
        DO k = 1, npwx*npol
           psi(k,j,ib) = evc(k,j,ib)
        END DO
     END DO
  END DO
  !
  ! ... hpsi contains h times the basis vectors (all the k-points, one call)
  !
  CALL start_clock( 'cegterg:h_psi_ptr init' )
  !$acc kernels async(async_id)
  psi_fft_comp(:,:,:) = psi(1:(npwx*npol), 1:nvecx, :)
  hpsi_fft_comp(:,:,:) = hpsi(1:(npwx*npol), 1:nvecx, :)
  !$acc end kernels
  DO ib = 1, n_k
     npw_comp(ib) = npw(ib)
  END DO
  !$acc update device(npw_comp(1:n_k)) async(async_id)
  CALL h_psi_ptr( npwx, npw_comp, nvecx, nvec, psi_fft_comp, hpsi_fft_comp, 1, n_k )
  !$acc kernels async(async_id)
  hpsi(:,1:nvecx,:) = hpsi_fft_comp(1:(npwx*npol),1:nvecx,:)
  psi(:,1:nvecx,:) = psi_fft_comp(1:(npwx*npol),1:nvecx,:)
  !$acc end kernels
  CALL stop_clock( 'cegterg:h_psi_ptr init' )
  nhpsi(:) = nhpsi(:) + nvec
  !
  ! ... spsi contains s times the basis vectors
  !
  IF ( uspp ) THEN
     DO ib = 1, n_k
        CALL s_psi_ptr( npwx, npw(ib), nvec, psi(1,1,ib), spsi(1,1,ib), ib )
     END DO
  END IF
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced
  ! ... space vc contains the eigenvectors of hc
  !
  CALL start_clock( 'cegterg:init' )
  !
  CALL device_blas_scalars_init()
  kdim_max = MAXVAL(kdim(1:n_k))
  DO ib = 1, n_k
     ptr_hc(ib) = c_devloc(hc(1,1,ib))
     ptr_sc(ib) = c_devloc(sc(1,1,ib))
     ptr_psi(ib) = acc_deviceptr(psi(1,1,ib))
     ptr_psi_result(ib) = acc_deviceptr(psi(1,1,ib))
     ptr_hpsi(ib) = acc_deviceptr(hpsi(1,1,ib))
     IF ( uspp ) ptr_spsi(ib) = acc_deviceptr(spsi(1,1,ib))
  END DO
  CALL start_clock('First Zgemm call init')
  !$acc update device(ptr_psi, ptr_hpsi, ptr_hc)
  !$acc host_data use_device(ptr_psi, ptr_hpsi, ptr_hc)
  info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_C, CUBLAS_OP_N, nvec, nvec, kdim_max, blas_cone, &
                            ptr_psi, kdmx, ptr_hpsi, kdmx, blas_czero, ptr_hc, nvecx_pad, n_k)
  !$acc end host_data
  CALL stop_clock('First Zgemm call init')
  CALL start_clock('Second Zgemm call init')
  IF ( uspp ) THEN
     !$acc update device(ptr_psi, ptr_spsi, ptr_sc)
     !$acc host_data use_device(ptr_psi, ptr_spsi, ptr_sc)
     info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_C, CUBLAS_OP_N, nvec, nvec, kdim_max, blas_cone, &
                               ptr_psi, kdmx, ptr_spsi, kdmx, blas_czero, ptr_sc, nvecx_pad, n_k)
     !$acc end host_data
  ELSE
     !$acc update device(ptr_psi, ptr_psi_result, ptr_sc)
     !$acc host_data use_device(ptr_psi, ptr_psi_result, ptr_sc)
     info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_C, CUBLAS_OP_N, nvec, nvec, kdim_max, blas_cone, &
                               ptr_psi, kdmx, ptr_psi_result, kdmx, blas_czero, ptr_sc, nvecx_pad, n_k)
     !$acc end host_data
  END IF
  CALL stop_clock('Second Zgemm call init')
  !
  DO ib = 1, n_k
     !$acc parallel vector_length(64) async(async_id)
     !$acc loop gang
     DO n = 1, nbase(ib)
        ! ... the diagonal of hc and sc must be strictly real
        hc(n,n,ib) = CMPLX( REAL( hc(n,n,ib) ), 0.D0 ,kind=DP)
        sc(n,n,ib) = CMPLX( REAL( sc(n,n,ib) ), 0.D0 ,kind=DP)
        !$acc loop vector
        DO m = n + 1, nbase(ib)
           hc(n,m,ib) = CONJG( hc(m,n,ib) )
           sc(n,m,ib) = CONJG( sc(m,n,ib) )
        END DO
     END DO
     !$acc end parallel
  END DO
  !
  CALL stop_clock( 'cegterg:init' )
  !
  IF ( lrot ) THEN
     !
     DO ib = 1, n_k
        !$acc kernels async(async_id)
        vc(1:nbase(ib),1:nbase(ib),ib) = ZERO
        !$acc end kernels
        !$acc parallel loop async(async_id)
        DO n = 1, nbase(ib)
           e(n,ib) = REAL( hc(n,n,ib) )
           vc(n,n,ib) = ONE
        END DO
     END DO
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian (all the k-points with one batched call)
     !
     CALL start_clock( 'cegterg:diag' )
     CALL nvtxRangePush("cdiaghg batched 1")
     !$acc kernels async(async_id)
     hc_comp(:,:,:) = hc(1:nvec,1:nvec,:)
     sc_comp(:,:,:) = sc(1:nvec,1:nvec,:)
     !$acc end kernels
     !
     !$acc host_data use_device(hc_comp, sc_comp, ew_comp, vc_comp)
     CALL diaghg( nvec, nvec, hc_comp, sc_comp, nvec, ew_comp, vc_comp, n_k, me_bgrp, root_bgrp, intra_bgrp_comm )
     !$acc end host_data
     !
     ! ... retrieve the eigenvectors/eigenvalues
     !
     !$acc kernels async(async_id)
     vc(1:nvec,1:nvec,:) = vc_comp(:,:,:)
     ew(1:nvec,:)        = ew_comp(:,:)
     ew_old_comp(1:nvec,:) = ew_comp(:,:)
     e(1:nvec,:) = ew_comp(:,:)
     !$acc end kernels
     CALL stop_clock( 'cegterg:diag' )
     CALL nvtxRangePop()
  END IF
  !
  ! ... iterate
  !
  dav_iter(:) = 0
  nbase_comp(:) = nvec
  notcnv_comp(:) = nvec
  my_n_comp(:) = nvec
  my_n(:) = nvec

  iterate: DO kter = 1, maxter
     !
     ! ... reorder eigenvectors so that coefficients for unconverged roots come first. This allows to use quick
     ! ... matrix-matrix multiplications to set a new basis vector (see below). The list of the unconverged roots is
     ! ... built on the host and the reordering is done with two kernels (instead of two per root): through a
     ! ... temporary, because the destination column of a root can be the source of a later one.
     ! ... ew(nbase+1:nbase+notcnv) gets the eigenvalues of these roots, for use in g_psi_ptr.
     !
     CALL start_clock( 'cegterg:update' )
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        dav_iter(ib) = kter
        np = 0
        DO n = 1, nvec
           IF ( .NOT. conv(n,ib) ) THEN
              np = np + 1
              src_idx(np,ib) = n
           END IF
        END DO
        !
        nb1(ib) = nbase(ib) + 1
        !
        ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
        !
        my_n(ib) = nbase(ib)
     END DO
     DO ib = 1, n_k
        par(ip_nbase,ib) = nbase(ib)
        par(ip_nb1,ib) = nb1(ib)
        par(ip_notcnv,ib) = notcnv(ib)
        par(ip_npw,ib) = npw(ib)
        par(ip_act,ib) = MERGE( 0, 1, my_done(ib) )
     END DO
     !$acc update device(par, src_idx)
     !
     ! the reordering of all the active k-points: two kernels in total (np of a k-point is notcnv of that k-point)
     !$acc parallel loop collapse(3) async(async_id) present(vc, vc_tmp, src_idx, par)
     DO ib = 1, n_k
        DO j = 1, nvec
           DO i = 1, nvecx
              IF ( par(ip_act,ib) == 1 .AND. j <= par(ip_notcnv,ib) ) vc_tmp(i,j,ib) = vc(i,src_idx(j,ib),ib)
           END DO
        END DO
     END DO
     !$acc parallel loop collapse(3) async(async_id) present(vc, vc_tmp, ew, e, src_idx, par)
     DO ib = 1, n_k
        DO j = 1, nvec
           DO i = 1, nvecx
              IF ( par(ip_act,ib) == 1 .AND. j <= par(ip_notcnv,ib) ) THEN
                 vc(i,j,ib) = vc_tmp(i,j,ib)
                 IF ( i == 1 ) ew(par(ip_nbase,ib)+j,ib) = e(src_idx(j,ib),ib)
              END IF
           END DO
        END DO
     END DO
     !
     n_active_zgem = COUNT(.NOT. done_comp(1:n_k))
     DO ib = 1, n_k
        IF ( .NOT. my_done(ib) ) THEN
           notcnv_comp(ib) = notcnv(ib)
           my_n_comp(ib) = my_n(ib)
        END IF
     END DO
     !
     IF ( n_active_zgem > 0 ) THEN
        notcnv_max = MAXVAL(notcnv_comp(1:n_k), MASK=.NOT. done_comp(1:n_k))
        my_n_max = MAXVAL(my_n_comp(1:n_k), MASK=.NOT. done_comp(1:n_k))
     END IF
     !
     ! ... The batched Zgemm below contracts over my_n_max (the maximum over the active
     ! ... batches), but a batch only has my_n valid basis vectors (e.g. right after a
     ! ... restart, when another batch still has a larger basis). The contribution of the
     ! ... extra rows of vc must vanish: zero them (the stale values there are not zero).
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        IF ( my_n_max > my_n(ib) ) THEN
           !$acc kernels async(async_id) present(vc)
           vc(my_n(ib)+1:my_n_max, 1:notcnv_max, ib) = ZERO
           !$acc end kernels
        END IF
     END DO
     !
     ! pointers of the active slots
     DO ib = 1, n_k
        IF ( .NOT. my_done(ib) .AND. .NOT. done_comp(ib) ) THEN
           my_slot = 1 + COUNT(.NOT. done_comp(1:ib-1))
           ptr_vc(my_slot) = c_devloc(vc(1,1,ib))
           ptr_hc(my_slot) = c_devloc(hc(nb1(ib),1,ib))
           ptr_sc(my_slot) = c_devloc(sc(nb1(ib),1,ib))
           ptr_psi(my_slot) = acc_deviceptr(psi(1,1,ib))
           ptr_hpsi(my_slot) = acc_deviceptr(hpsi(1,1,ib))
           ptr_hpsi_nb1(my_slot) = acc_deviceptr(hpsi(1,nb1(ib),ib))
           ptr_psi_result(my_slot) = acc_deviceptr(psi(1,nb1(ib),ib))
           IF ( uspp ) ptr_spsi(my_slot) = acc_deviceptr(spsi(1,1,ib))
        END IF
     END DO
     !
     IF ( uspp ) THEN
        !$acc update device(ptr_vc, ptr_psi, ptr_spsi, ptr_psi_result)
        !$acc host_data use_device(ptr_psi_result, ptr_spsi, ptr_vc)
        info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_N, CUBLAS_OP_N, kdim_max, notcnv_max, my_n_max, blas_cone, &
                                  ptr_spsi, kdmx, ptr_vc, nvecx, blas_czero, ptr_psi_result, kdmx, n_active_zgem)
        !$acc end host_data
     ELSE
        !$acc update device(ptr_vc, ptr_psi, ptr_psi_result)
        !$acc host_data use_device(ptr_psi, ptr_psi_result, ptr_vc)
        info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_N, CUBLAS_OP_N, kdim_max, notcnv_max, my_n_max, blas_cone, &
                                  ptr_psi, kdmx, ptr_vc, nvecx, blas_czero, ptr_psi_result, kdmx, n_active_zgem)
        !$acc end host_data
     END IF
     !
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        !
        ! NB: no mp_sum over inter_bgrp_comm here (band groups are not supported)
        !$acc parallel loop collapse(3) async(async_id)
        DO np = 1, notcnv(ib)
           DO ipol = 1, npol
              DO k = 1, npwx
                 psi(k + (ipol-1)*npwx, nbase(ib)+np, ib) = - ew(nbase(ib)+np,ib)*psi(k + (ipol-1)*npwx, nbase(ib)+np, ib)
              END DO
           END DO
        END DO
        !$acc end parallel
     END DO
     !
     !$acc update device(ptr_hpsi, ptr_psi_result)
     !$acc host_data use_device(ptr_hpsi, ptr_vc, ptr_psi_result)
     info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_N, CUBLAS_OP_N, kdim_max, notcnv_max, my_n_max, blas_cone, &
                               ptr_hpsi, kdmx, ptr_vc, nvecx, blas_cone, ptr_psi_result, kdmx, n_active_zgem)
     !$acc end host_data
     !
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        ! clean up garbage if there is any
        IF ( npw(ib) < npwx ) THEN
           !$acc parallel loop collapse(2) async(async_id)
           DO j = nb1(ib), nbase(ib)+notcnv(ib)
              DO k = npw(ib)+1, npwx
                 psi(k,j,ib) = ZERO
              END DO
           END DO
        END IF
        IF ( npol == 2 ) THEN
           !$acc parallel loop collapse(2) async(async_id)
           DO j = nb1(ib), nbase(ib)+notcnv(ib)
              DO k = npwx+npw(ib)+1, 2*npwx
                 psi(k,j,ib) = ZERO
              END DO
           END DO
        END IF
     END DO
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration (g_psi) and the hpsi and spsi of the new functions
     !
     CALL start_clock( 'cegterg:h_psi_ptr iter' )
     !
     n_active = COUNT(.NOT. done_comp(1:n_k))
     ! slot tables on the host (one transfer), then one kernel for the copies of all the active slots
     n_slot = 0
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        n_slot = n_slot + 1
        npw_comp(n_slot) = npw(ib)
        nb1_comp(n_slot) = nb1(ib)
        slot2batch(n_slot) = ib
     END DO
     !$acc update device(npw_comp(1:n_active), nb1_comp(1:n_active), slot2batch(1:n_active))
     nb_slot_max = MAXVAL( notcnv(1:n_k), MASK=.NOT. my_done(1:n_k) )
     !$acc parallel loop collapse(3) async(async_id) present(psi_fft_comp, hpsi_fft_comp, psi, hpsi, slot2batch)
     DO n = 1, n_active
        DO j = 1, nvecx
           DO k = 1, npwx*npol
              psi_fft_comp(k,j,n) = psi(k,j,slot2batch(n))
              hpsi_fft_comp(k,j,n) = hpsi(k,j,slot2batch(n))
           END DO
        END DO
     END DO
     !$acc parallel loop collapse(2) async(async_id) present(eig_new_comp, ew, par, slot2batch)
     DO n = 1, n_active
        DO i = 1, nb_slot_max
           IF ( i <= par(ip_notcnv,slot2batch(n)) ) eig_new_comp(i,n) = ew(par(ip_nb1,slot2batch(n))+i-1,slot2batch(n))
        END DO
     END DO
     !
     ! notcnv_comp is indexed by batch on the host, while h_psi_ptr_iter works with
     ! slots (the active batches, packed): build the slot-indexed copy explicitly.
     ! (it is a plain host array: the routine copies it to the device by itself)
     n_slot = 0
     DO n = 1, n_k
        IF ( .NOT. done_comp(n) ) THEN
           n_slot = n_slot + 1
           notcnv_slot(n_slot) = notcnv_comp(n)
        END IF
     END DO
     ! approximate inverse iteration for all the slots (the callback works on psi_fft_comp)
     CALL g_psi_ptr_iter(npwx, npw_comp, nvecx, notcnv_slot, nb1_comp, psi_fft_comp, eig_new_comp, nvec, async_id, n_active, slot2batch)
     !
     ! "normalize" the correction vectors of all the slots in order to improve numerical stability of the subspace
     ! diagonalization: psi(:,nb1:nb1+notcnv-1) / sqrt(<psi|psi>)  (one gang for every (slot, band) pair)
     nb_slot_max = MAXVAL( notcnv_slot(1:n_active) )
     !$acc parallel loop gang collapse(2) async(async_id) copyin(notcnv_slot(1:n_active)) private(nrm, col)
     DO n = 1, n_active
        DO i = 1, nb_slot_max
           IF ( i <= notcnv_slot(n) ) THEN
              col = nb1_comp(n) + i - 1
              nrm = 0.D0
              !$acc loop vector collapse(2) reduction(+:nrm)
              DO ipol = 1, npol
                 DO k = 1, npw_comp(n)
                    nrm = nrm + DBLE( psi_fft_comp(k+(ipol-1)*npwx,col,n) )**2 + AIMAG( psi_fft_comp(k+(ipol-1)*npwx,col,n) )**2
                 END DO
              END DO
              !$acc loop vector collapse(2)
              DO ipol = 1, npol
                 DO k = 1, npw_comp(n)
                    psi_fft_comp(k+(ipol-1)*npwx,col,n) = psi_fft_comp(k+(ipol-1)*npwx,col,n) / SQRT( nrm )
                 END DO
              END DO
           END IF
        END DO
     END DO
     CALL h_psi_ptr_iter(npwx, npw_comp, nvecx, notcnv_slot, nb1_comp, psi_fft_comp, hpsi_fft_comp, async_id, n_active, slot2batch)
     CALL stop_clock( 'cegterg:h_psi_ptr iter' )
     !
     DO ib = 1, n_k
        IF ( .NOT. my_done(ib) ) nhpsi(ib) = nhpsi(ib) + notcnv(ib)
     END DO
     !$acc parallel loop collapse(3) async(async_id) present(psi_fft_comp, hpsi_fft_comp, psi, hpsi, slot2batch)
     DO n = 1, n_active
        DO j = 1, nvecx
           DO k = 1, npwx*npol
              hpsi(k,j,slot2batch(n)) = hpsi_fft_comp(k,j,n)
              psi(k,j,slot2batch(n)) = psi_fft_comp(k,j,n)
           END DO
        END DO
     END DO
     !
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        IF ( uspp ) CALL s_psi_ptr( npwx, npw(ib), notcnv(ib), psi(1,nb1(ib),ib), spsi(1,nb1(ib),ib), ib )
        my_n(ib) = nbase(ib) + notcnv(ib)
     END DO
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
     ! We need to recompute notcnv_comp and my_n_comp because they were re-computed since the first assignment:
     DO ib = 1, n_k
        IF ( .NOT. my_done(ib) ) THEN
           notcnv_comp(ib) = notcnv(ib)
           my_n_comp(ib) = my_n(ib)
        END IF
     END DO
     notcnv_max = MAXVAL(notcnv_comp(:))
     my_n_max = MAXVAL(my_n_comp(:))
     !
     !$acc update device(ptr_hpsi_nb1, ptr_psi, ptr_hc)
     !$acc host_data use_device(ptr_hpsi_nb1, ptr_psi, ptr_hc)
     info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_C, CUBLAS_OP_N, notcnv_max, my_n_max, kdim_max, blas_cone, &
                               ptr_hpsi_nb1, kdmx, ptr_psi, kdmx, blas_czero, ptr_hc, nvecx_pad, n_active_zgem)
     !$acc end host_data
     !
     ! The pointers to the matrices are re-assigned (as the dimensions changed)
     DO ib = 1, n_k
        IF ( .NOT. my_done(ib) .AND. .NOT. done_comp(ib) ) THEN
           my_slot = 1 + COUNT(.NOT. done_comp(1:ib-1))
           ptr_sc(my_slot) = c_devloc(sc(nb1(ib),1,ib))
           ptr_psi(my_slot) = acc_deviceptr(psi(1,1,ib))
           ptr_psi_result(my_slot) = acc_deviceptr(psi(1,nb1(ib),ib))
           IF ( uspp ) ptr_spsi(my_slot) = acc_deviceptr(spsi(1,nb1(ib),ib))
        END IF
     END DO
     !
     IF ( uspp ) THEN
        !$acc update device(ptr_spsi, ptr_psi, ptr_sc)
        !$acc host_data use_device(ptr_spsi, ptr_psi, ptr_sc)
        info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_C, CUBLAS_OP_N, notcnv_max, my_n_max, kdim_max, blas_cone, &
                                  ptr_spsi, kdmx, ptr_psi, kdmx, blas_czero, ptr_sc, nvecx_pad, n_active_zgem)
        !$acc end host_data
     ELSE
        !$acc update device(ptr_psi_result, ptr_psi, ptr_sc)
        !$acc host_data use_device(ptr_psi_result, ptr_psi, ptr_sc)
        info = cublasZgemmBatched(cublas_handle(async_id), CUBLAS_OP_C, CUBLAS_OP_N, notcnv_max, my_n_max, kdim_max, blas_cone, &
                                  ptr_psi_result, kdmx, ptr_psi, kdmx, blas_czero, ptr_sc, nvecx_pad, n_active_zgem)
        !$acc end host_data
     END IF
     !
     CALL stop_clock( 'cegterg:overlap' )
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        nbase(ib) = nbase(ib) + notcnv(ib)
        !
        !$acc parallel vector_length(64) async(async_id)
        !$acc loop gang
        DO n = 1, nbase(ib)
           ! ... the diagonal of hc and sc must be strictly real
           IF( n>=nb1(ib) ) THEN
              hc(n,n,ib) = CMPLX( REAL( hc(n,n,ib) ), 0.D0 ,kind=DP)
              sc(n,n,ib) = CMPLX( REAL( sc(n,n,ib) ), 0.D0 ,kind=DP)
           ENDIF
           !$acc loop vector
           DO m = MAX(n+1,nb1(ib)), nbase(ib)
              hc(n,m,ib) = CONJG( hc(m,n,ib) )
              sc(n,m,ib) = CONJG( sc(m,n,ib) )
           END DO
        END DO
        !$acc end parallel
        !
        ! The current nbase is stored only if the convergence has not been reached
        IF ( .NOT. done_comp(ib) ) nbase_comp(ib) = nbase(ib)
     END DO
     !
     ! ... diagonalize the reduced hamiltonian (all the active slots with one batched call)
     !
     CALL nvtxRangePush("Cdiaghg call 2")
     n_active = COUNT(.NOT. done_comp(1:n_k))
     IF ( n_active > 0 ) THEN
        nbase_max = MAXVAL(nbase_comp(1:n_k), MASK=.NOT. done_comp(1:n_k))
     END IF
     IF ( prev_nbase_max /= nbase_max .OR. prev_n_active /= n_active ) THEN
        IF ( prev_nbase_max > 0 ) THEN
           !$acc exit data delete(hc_comp_itr, sc_comp_itr, vc_comp_itr, ew_comp_itr)
           DEALLOCATE( hc_comp_itr, sc_comp_itr, vc_comp_itr, ew_comp_itr )
        END IF
        ALLOCATE( hc_comp_itr(nbase_max,nbase_max,n_active), sc_comp_itr(nbase_max,nbase_max,n_active), &
                  vc_comp_itr(nbase_max,nbase_max,n_active), ew_comp_itr(nbase_max,n_active) )
        !$acc enter data create(hc_comp_itr, sc_comp_itr, vc_comp_itr, ew_comp_itr)
        prev_nbase_max = nbase_max
        prev_n_active  = n_active
     END IF
     DO ib = 1, n_k
        IF ( my_done(ib) .OR. done_comp(ib) ) CYCLE
        my_slot = 1 + COUNT(.NOT. done_comp(1:ib-1))
        ! hc_comp_itr / sc_comp_itr of this slot, nbase_max x nbase_max, in one kernel:
        !  - the nbase x nbase block is copied from hc / sc
        !  - the padding (rows and columns beyond nbase) is zero, except for the diagonal,
        !    where sc = 1 and hc = i*1e3 so that the padded eigenpairs are decoupled and sorted last
        !$acc parallel loop collapse(2) async(async_id)
        DO j = 1, nbase_max
           DO i = 1, nbase_max
              IF ( i <= nbase(ib) .AND. j <= nbase(ib) ) THEN
                 hc_comp_itr(i,j,my_slot) = hc(i,j,ib)
                 sc_comp_itr(i,j,my_slot) = sc(i,j,ib)
              ELSE IF ( i == j ) THEN
                 hc_comp_itr(i,j,my_slot) = CMPLX( REAL(i,DP)*1.D3, 0.D0, kind=DP )
                 sc_comp_itr(i,j,my_slot) = CMPLX( 1.D0, 0.D0, kind=DP )
              ELSE
                 hc_comp_itr(i,j,my_slot) = CMPLX( 0.D0, 0.D0, kind=DP )
                 sc_comp_itr(i,j,my_slot) = CMPLX( 0.D0, 0.D0, kind=DP )
              END IF
           END DO
        END DO
     END DO
     !
     CALL start_clock( 'cegterg:diag' )
     IF ( n_active > 0 ) THEN
        !$acc host_data use_device(hc_comp_itr, sc_comp_itr, ew_comp_itr, vc_comp_itr)
        CALL diaghg( nbase_max, nvec, hc_comp_itr, sc_comp_itr, nbase_max, ew_comp_itr, vc_comp_itr, &
                     n_active, me_bgrp, root_bgrp, intra_bgrp_comm )
        !$acc end host_data
        !$acc wait(async_id)
        !
        ! ... test for convergence of all the slots at once: conv_comp(:,batch) is true for the roots whose
        ! ... eigenvalue changed by less than the threshold since the previous iteration
        !$acc parallel loop collapse(2) present(ew_comp_itr, ew_old_comp, conv_comp, btype, slot2batch)
        DO n = 1, n_active
           DO i = 1, nvec
              IF ( btype(i) == 1 ) THEN
                 thr_conv = ethr
              ELSE
                 thr_conv = empty_ethr
              END IF
              IF ( ABS( ew_comp_itr(i,n) - ew_old_comp(i,slot2batch(n)) ) < thr_conv ) THEN
                 conv_comp(i,slot2batch(n)) = 1
              ELSE
                 conv_comp(i,slot2batch(n)) = 0
              END IF
              ew_old_comp(i,slot2batch(n)) = ew_comp_itr(i,n)
           END DO
        END DO
        !$acc update self(conv_comp(1:nvec,1:n_k))
     END IF
     !
     ! We retrieve data:
     DO ib = 1, n_k
        IF ( my_done(ib) .OR. done_comp(ib) ) CYCLE
        my_slot = 1 + COUNT(.NOT. done_comp(1:ib-1))
        !$acc kernels async(async_id)
        vc(1:nbase_max,1:nbase_max,ib) = vc_comp_itr(1:nbase_max,1:nbase_max,my_slot)
        ew(1:nbase_max,ib)        = ew_comp_itr(1:nbase_max,my_slot)
        !$acc end kernels
     END DO
     CALL nvtxRangePop()
     CALL stop_clock( 'cegterg:diag' )
     !
     ! ... test for convergence
     !
     ! the test was done for all the slots right after the batched diagonalization: take the column of this batch
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        conv(1:nvec,ib) = ( conv_comp(1:nvec,ib) == 1 )
        notcnv(ib) = COUNT( .NOT. conv(:,ib) )
        !$acc kernels async(async_id)
        e(1:nvec,ib) = ew(1:nvec,ib)
        !$acc end kernels
     END DO
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     CALL start_clock( 'cegterg:last' )
     restart_now(:) = .FALSE.
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        restart_now(ib) = ( notcnv(ib) == 0 .OR. nbase(ib)+notcnv(ib) > nvecx .OR. dav_iter(ib) == maxter )
        IF ( restart_now(ib) ) THEN
           my_n(ib) = nbase(ib)
           !$acc host_data use_device(evc, psi, vc)
           ! the product is issued on the stream of the queue, so it is ordered with the kernels enqueued on it
           CALL MYZGEMM_H( cublas_handle(async_id), 'N','N', kdim(ib), nvec, my_n(ib), ONE, psi(1,1,ib), kdmx, &
                           vc(1,1,ib), nvecx, ZERO, evc(1,1,ib), kdmx )
           !$acc end host_data
           IF ( notcnv(ib) == 0 ) my_done(ib) = .TRUE.
        END IF
     END DO
     !
     ! Update of done_comp and exit test
     done_comp(1:n_k) = my_done(1:n_k)
     IF ( ALL(done_comp(1:n_k)) .OR. kter == maxter ) THEN
        CALL stop_clock( 'cegterg:last' )
        EXIT iterate
     END IF
     !
     DO ib = 1, n_k
        IF ( my_done(ib) ) CYCLE
        IF ( restart_now(ib) ) THEN
           !
           ! ... refresh psi, H*psi and S*psi
           !
           !$acc parallel loop collapse(2) async(async_id) present(psi, evc)
           DO j = 1, nvec
              DO k = 1, npwx*npol
                 psi(k,j,ib) = evc(k,j,ib)
              END DO
           END DO
           !
           IF ( uspp ) THEN
              !$acc host_data use_device(psi, spsi, vc)
              CALL MYZGEMM_H( cublas_handle(async_id), 'N','N', kdim(ib), nvec, my_n(ib), ONE, spsi(1,1,ib), kdmx, &
                              vc(1,1,ib), nvecx, ZERO, psi(1,nvec+1,ib), kdmx )
              !$acc end host_data
              !$acc parallel loop collapse(2) async(async_id) present(psi, spsi)
              DO j = 1, nvec
                 DO k = 1, npwx*npol
                    spsi(k,j,ib) = psi(k,nvec+j,ib)
                 END DO
              END DO
           END IF
           !
           !$acc host_data use_device(psi, hpsi, vc)
           CALL MYZGEMM_H( cublas_handle(async_id), 'N','N', kdim(ib), nvec, my_n(ib), ONE, hpsi(1,1,ib), kdmx, &
                           vc(1,1,ib), nvecx, ZERO, psi(1,nvec+1,ib), kdmx )
           !$acc end host_data
           !$acc parallel loop collapse(2) async(async_id) present(psi, hpsi)
           DO j = 1, nvec
              DO k = 1, npwx*npol
                 hpsi(k,j,ib) = psi(k,nvec+j,ib)
              END DO
           END DO
           !
           ! ... refresh the reduced hamiltonian
           !
           nbase(ib) = nvec
           !
           !$acc kernels async(async_id)
           DO n = 1, nbase(ib)
              hc(n,n,ib) = CMPLX( e(n,ib), 0.0_DP ,kind=DP)
              sc(n,n,ib) = ONE
              vc(n,n,ib) = ONE
              DO j = n+1, nbase(ib)
                 hc(j,n,ib) = ZERO
                 hc(n,j,ib) = ZERO
                 sc(j,n,ib) = ZERO
                 sc(n,j,ib) = ZERO
                 vc(j,n,ib) = ZERO
                 vc(n,j,ib) = ZERO
              END DO
           END DO
           !$acc end kernels
        END IF
     END DO
     CALL stop_clock( 'cegterg:last' )
     !
  END DO iterate
  !
  ! make sure that everything the caller will read has been computed
  !$acc wait(async_id)
  !
  IF ( prev_nbase_max > 0 ) THEN
     !$acc exit data delete(hc_comp_itr, sc_comp_itr, vc_comp_itr, ew_comp_itr)
     DEALLOCATE( hc_comp_itr, sc_comp_itr, vc_comp_itr, ew_comp_itr )
  END IF
  !$acc exit data delete(psi_fft_comp, hpsi_fft_comp, npw_comp, nb1_comp, slot2batch, conv_comp, ew_old_comp, eig_new_comp, &
  !$acc hc_comp, sc_comp, vc_comp, ew_comp, ptr_psi, ptr_psi_result, ptr_spsi, ptr_hpsi, ptr_hpsi_nb1, ptr_hc, ptr_sc, ptr_vc)
  DEALLOCATE( psi_fft_comp, hpsi_fft_comp, npw_comp, nb1_comp, slot2batch, conv_comp, ew_old_comp, eig_new_comp )
  DEALLOCATE( hc_comp, sc_comp, vc_comp, ew_comp )
  DEALLOCATE( ptr_psi, ptr_psi_result, ptr_hpsi, ptr_hpsi_nb1, ptr_spsi, ptr_hc, ptr_sc, ptr_vc )
  !
  !$acc exit data delete(ew, par, src_idx) async(async_id)
  DEALLOCATE( ew, par, src_idx )
  DEALLOCATE( vc )
  DEALLOCATE( vc_tmp )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  !
  IF ( uspp ) THEN
     !$acc exit data async(async_id) delete(spsi)
     DEALLOCATE( spsi )
  END IF
  !
  !$acc exit data async(async_id) delete(psi, hpsi)
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  !
  !$acc end data
  !
  CALL stop_clock( 'cegterg' )
  !
  RETURN
  !
END SUBROUTINE cegterg
!
!  Subroutine with distributed matrixes
!  (written by Carlo Cavazzoni)
!
!----------------------------------------------------------------------------
SUBROUTINE pcegterg(h_psi_ptr, s_psi_ptr, uspp, g_psi_ptr, &
                    npw, npwx, nvec, nvecx, npol, evc, ethr, &
                    e, btype, notcnv, lrot, dav_iter, nhpsi )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  !
  USE util_param,       ONLY : DP, stdout
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier, &
                               mp_size, mp_type_free, mp_allgather
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! number of spin polarizations
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
    ! chunking parameters
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx*npol,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : S|psi> not needed
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  INTEGER, INTENT(OUT) :: nhpsi
    ! total number of indivitual hpsi
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, ipol, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  INTEGER :: ierr
  REAL(DP), ALLOCATABLE :: ew(:)
  COMPLEX(DP), ALLOCATABLE :: hl(:,:), sl(:,:), vl(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr
    ! threshold for empty bands
  INTEGER :: idesc(LAX_DESC_SIZE), idesc_old(LAX_DESC_SIZE)
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: notcnv_ip( : )
  INTEGER, ALLOCATABLE :: ic_notcnv( : )
  !
  INTEGER :: np_ortho(2), ortho_parent_comm
  LOGICAL :: do_distr_diag_inside_bgrp
  !
  REAL(DP), EXTERNAL :: ddot
  !
  EXTERNAL  h_psi_ptr, s_psi_ptr, g_psi_ptr
    ! h_psi_ptr(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi_ptr(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  nhpsi = 0
  CALL start_clock( 'cegterg' )
  !
  CALL laxlib_getval( np_ortho = np_ortho, ortho_parent_comm = ortho_parent_comm, &
    do_distr_diag_inside_bgrp = do_distr_diag_inside_bgrp )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'pcegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  ! compute the number of chuncks
  numblock  = (npw+blocksize-1)/blocksize
  !
  ALLOCATE(  psi( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate psi ', ABS(ierr) )
  !
  ALLOCATE( hpsi( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx*npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ! ... Initialize the matrix descriptor
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ic_notcnv ', ABS(ierr) )
  !
  ALLOCATE( notcnv_ip( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate notcnv_ip ', ABS(ierr) )
  !
  CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
  !
  IF( la_proc ) THEN
     !
     ! only procs involved in the diagonalization need to allocate local
     ! matrix block.
     !
     ALLOCATE( vl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  ELSE
     !
     ALLOCATE( vl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ew ', ABS(ierr) )
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  CALL threaded_memcpy(psi, evc, nvec*npol*npwx*2)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_ptr( npwx, npw, nvec, psi, hpsi ) ; nhpsi = nhpsi + nvec
  !
  IF ( uspp ) CALL s_psi_ptr( npwx, npw, nvec, psi, spsi )
  !
  ! ... hl contains the projection of the hamiltonian onto the reduced
  ! ... space, vl contains the eigenvectors of hl. Remember hl, vl and sl
  ! ... are all distributed across processors, global replicated matrixes
  ! ... here are never allocated
  !
  CALL start_clock( 'cegterg:init' )

  CALL compute_distmat( hl, psi, hpsi )
  !
  IF ( uspp ) THEN
     !
     CALL compute_distmat( sl, psi, spsi )
     !
  ELSE
     !
     CALL compute_distmat( sl, psi, psi )
     !
  END IF
  CALL stop_clock( 'cegterg:init' )
  !
  !
  IF ( lrot ) THEN
     !
     CALL set_e_from_h()
     !
     CALL set_to_identity( vl, idesc )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !     Calling block parallel algorithm
     !
     CALL start_clock( 'cegterg:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
     END IF
     CALL stop_clock( 'cegterg:diag' )
     !
     e(1:nvec) = ew(1:nvec)
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter ; !write(*,*) kter, notcnv, conv
     !
     CALL start_clock( 'cegterg:update' )
     !
     CALL reorder_v()
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL hpsi_dot_v()
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi_ptr( npwx, npw, notcnv, npol, psi(1,nb1), ew(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
           ew(n) = ddot( 2*npw, psi(1,nbn), 1, psi(1,nbn), 1 )
           !
        ELSE
           !
           ew(n) = ddot( 2*npw, psi(1,nbn), 1, psi(1,nbn), 1 ) + &
                   ddot( 2*npw, psi(npwx+1,nbn), 1, psi(npwx+1,nbn), 1 )
           !
        END IF
        !
     END DO
     !
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !
     !$omp parallel do collapse(3)
     DO n = 1, notcnv
        DO ipol = 1, npol
           DO m = 1, numblock
              psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw, m*blocksize)+(ipol-1)*npwx,nbase+n) = &
              psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw, m*blocksize)+(ipol-1)*npwx,nbase+n) / &
                    SQRT( ew(n) )
           END DO
        END DO
     END DO
     !$omp end parallel do
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi_ptr( npwx, npw, notcnv, psi(1,nb1), hpsi(1,nb1) ) ; nhpsi = nhpsi + notcnv
     !
     IF ( uspp ) CALL s_psi_ptr( npwx, npw, notcnv, psi(1,nb1), spsi(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
     ! we need to save the old descriptor in order to redistribute matrices
     !
     idesc_old = idesc
     !
     ! ... RE-Initialize the matrix descriptor
     !
     CALL desc_init( nbase+notcnv, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
     !
     IF( la_proc ) THEN

        !  redistribute hl and sl (see dsqmred), since the dimension of the subspace has changed
        !
        vl = hl
        DEALLOCATE( hl )
        ALLOCATE( hl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )

        CALL laxlib_zsqmred( nbase, vl, idesc_old(LAX_DESC_NRCX), idesc_old, nbase+notcnv, hl, nx, idesc )

        vl = sl
        DEALLOCATE( sl )
        ALLOCATE( sl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )

        CALL laxlib_zsqmred( nbase, vl, idesc_old(LAX_DESC_NRCX), idesc_old, nbase+notcnv, sl, nx, idesc )

        DEALLOCATE( vl )
        ALLOCATE( vl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )

     END IF
     !
     !
     CALL update_distmat( hl, psi, hpsi )
     !
     IF ( uspp ) THEN
        !
        CALL update_distmat( sl, psi, spsi )
        !
     ELSE
        !
        CALL update_distmat( sl, psi, psi )
        !
     END IF
     !
     CALL stop_clock( 'cegterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     ! ... diagonalize the reduced hamiltonian
     !     Call block parallel algorithm
     !
     CALL start_clock( 'cegterg:diag' )
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pdiaghg ew and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id ) CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
           CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
           CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pdiaghg( nbase, hl, sl, nx, ew, vl, idesc )
     END IF
     CALL stop_clock( 'cegterg:diag' )
     !
     ! ... test for convergence
     !
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'cegterg:last' )
        !
        CALL refresh_evc()
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        CALL threaded_memcpy(psi, evc, nvec*npol*npwx*2)
        !
        IF ( uspp ) THEN
           !
           CALL refresh_spsi()
           !
        END IF
        !
        CALL refresh_hpsi()
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        CALL desc_init( nvec, nx, la_proc, idesc, rank_ip, irc_ip, nrc_ip )
        !
        IF( la_proc ) THEN
           !
           ! note that nx has been changed by desc_init
           ! we need to re-alloc with the new size.
           !
           DEALLOCATE( vl, hl, sl )
           ALLOCATE( vl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
           ALLOCATE( hl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
           ALLOCATE( sl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
           !
        END IF
        !
        CALL set_h_from_e( )
        !
        CALL set_to_identity( vl, idesc )
        CALL set_to_identity( sl, idesc )
        !
        CALL stop_clock( 'cegterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  !
  DEALLOCATE( vl, hl, sl )
  !
  DEALLOCATE( rank_ip )
  DEALLOCATE( irc_ip )
  DEALLOCATE( nrc_ip )
  DEALLOCATE( ic_notcnv )
  DEALLOCATE( notcnv_ip )
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )
  !
  CALL stop_clock( 'cegterg' )
  !call print_clock( 'cegterg' )
  !call print_clock( 'cegterg:init' )
  !call print_clock( 'cegterg:diag' )
  !call print_clock( 'cegterg:update' )
  !call print_clock( 'cegterg:overlap' )
  !call print_clock( 'cegterg:last' )
  !
  RETURN
  !
  !
CONTAINS
  !
  SUBROUTINE set_to_identity( distmat, idesc )
     INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
     COMPLEX(DP), INTENT(OUT) :: distmat(:,:)
     INTEGER :: i
     distmat = ( 0_DP , 0_DP )
     IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) .AND. idesc(LAX_DESC_ACTIVE_NODE) > 0 ) THEN
        DO i = 1, idesc(LAX_DESC_NC)
           distmat( i, i ) = ( 1_DP , 0_DP )
        END DO
     END IF
     RETURN
  END SUBROUTINE set_to_identity
  !
  !
  SUBROUTINE reorder_v()
     !
     INTEGER :: ipc
     INTEGER :: nc, ic
     INTEGER :: nl, npl
     !
     np = 0
     !
     notcnv_ip = 0
     !
     n = 0
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        npl = 0
        !
        IF( ic <= nvec ) THEN
           !
           DO nl = 1, min( nvec - ic + 1, nc )
              !
              n  = n  + 1
              !
              IF ( .NOT. conv(n) ) THEN
                 !
                 ! ... this root not yet converged ...
                 !
                 np  = np  + 1
                 npl = npl + 1
                 IF( npl == 1 ) ic_notcnv( ipc ) = np
                 !
                 ! ... reorder eigenvectors so that coefficients for unconverged
                 ! ... roots come first. This allows to use quick matrix-matrix
                 ! ... multiplications to set a new basis vector (see below)
                 !
                 notcnv_ip( ipc ) = notcnv_ip( ipc ) + 1
                 !
                 IF ( npl /= nl ) THEN
                    IF( la_proc .AND. idesc(LAX_DESC_MYC) == ipc-1 ) THEN
                       vl( :, npl) = vl( :, nl )
                    END IF
                 END IF
                 !
                 ! ... for use in g_psi_ptr
                 !
                 ew(nbase+np) = e(n)
                 !
              END IF
              !
           END DO
           !
        END IF
        !
     END DO
     !
  END SUBROUTINE reorder_v
  !
  !
  SUBROUTINE hpsi_dot_v()
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, ir, ic, notcl, root, np, ipol, ib
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP), ALLOCATABLE :: ptmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     ALLOCATE( ptmp( npwx*npol, nx ) )

     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        IF( notcnv_ip( ipc ) > 0 ) THEN

           notcl = notcnv_ip( ipc )
           ic    = ic_notcnv( ipc )

           beta = ZERO

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 vtmp(:,1:notcl) = vl(:,1:notcl)
              END IF

              CALL mp_bcast( vtmp(:,1:notcl), root, ortho_parent_comm )
              !
              IF ( uspp ) THEN
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    spsi(1, ir), kdmx, vtmp, nx, beta, psi(1,nb1+ic-1), kdmx )
                 !
              ELSE
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    psi(1, ir), kdmx, vtmp, nx, beta, psi(1,nb1+ic-1), kdmx )
                 !
              END IF
              !
              CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                      hpsi(1, ir), kdmx, vtmp, nx, beta, ptmp, kdmx )

              beta = ONE

           END DO

           !$omp parallel do collapse(3)
           DO np = 1, notcl
              DO ipol = 1, npol
                 DO ib = 1, numblock
                    !
                    psi( (ib-1)*blocksize+(ipol-1)*npwx+1: &
                         MIN(npw, ib*blocksize)+(ipol-1)*npwx,nbase+np+ic-1) = &
                    ptmp((ib-1)*blocksize+(ipol-1)*npwx+1: &
                         MIN(npw, ib*blocksize)+(ipol-1)*npwx,np) - &
                    ew(nbase+np+ic-1) * psi((ib-1)*blocksize+(ipol-1)*npwx+1:&
                       MIN(npw, ib*blocksize)+(ipol-1)*npwx,nbase+np+ic-1)
                    !
                 END DO
              END DO
           END DO
           !$omp end parallel do
           !
           ! clean up garbage if there is any
           IF (npw < npwx) psi(npw+1:npwx,nbase+ic:nbase+notcl+ic-1) = ZERO
           IF (npol == 2)  psi(npwx+npw+1:2*npwx,nbase+ic:nbase+notcl+ic-1) = ZERO
           !
        END IF
        !
     END DO

     DEALLOCATE( vtmp )
     DEALLOCATE( ptmp )

     RETURN
  END SUBROUTINE hpsi_dot_v
  !
  !
  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO

           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          psi(1,ir), kdmx, vl, nx, beta, evc(1,ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          psi(1,ir), kdmx, vtmp, nx, beta, evc(1,ic), kdmx )
              END IF
              !

              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_evc
  !
  !
  SUBROUTINE refresh_spsi( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO
           !
           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi(1,ir), kdmx, vl, nx, beta, psi(1,nvec+ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi(1,ir), kdmx, vtmp, nx, beta, psi(1,nvec+ic), kdmx )
              END IF
              !
              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     CALL threaded_memcpy(spsi, psi(1,nvec+1), nvec*npol*npwx*2)
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_spsi
  !
  !
  !
  SUBROUTINE refresh_hpsi( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = ZERO
           !
           DO ipr = 1, idesc(LAX_DESC_NPR)
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == idesc(LAX_DESC_MYR) .AND. ipc-1 == idesc(LAX_DESC_MYC) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi(1,ir), kdmx, vl, nx, beta, psi(1,nvec+ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi(1,ir), kdmx, vtmp, nx, beta, psi(1,nvec+ic), kdmx )
              END IF
              !
              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )
     !
     CALL threaded_memcpy(hpsi, psi(1,nvec+1), nvec*npol*npwx*2)
     !
     RETURN
  END SUBROUTINE refresh_hpsi
  !
  !
  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = ZERO
     !
     !  Only upper triangle is computed, then the matrix is hermitianized
     !
     DO ipc = 1, idesc(LAX_DESC_NPC) !  loop on column procs
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! idesc(LAX_DESC_NPR) ! ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE , &
                       v(1,ir), kdmx, w(1,ic), kdmx, ZERO, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     !  The matrix is hermitianized using upper triangle
     !
     CALL laxlib_zsqmher( nbase, dm, nx, idesc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE update_distmat( dm, v, w )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, icc, ii
     COMPLEX(DP) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )

     ALLOCATE( vtmp( nx, nx ) )
     !
     vtmp = ZERO
     !
     DO ipc = 1, idesc(LAX_DESC_NPC)
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic+nc-1 >= nb1 ) THEN
           !
           nc = MIN( nc, ic+nc-1 - nb1 + 1 )
           IF( ic >= nb1 ) THEN
              ii = ic
              icc = 1
           ELSE
              ii = nb1
              icc = nb1-ic+1
           END IF
           !
           ! icc to nc is the local index of the unconverged bands
           ! ii is the global index of the first unconverged bands
           !
           DO ipr = 1, ipc ! idesc(LAX_DESC_NPR) use symmetry
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE, v(1, ir), &
                          kdmx, w(1,ii), kdmx, ZERO, vtmp, nx )
              IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) vtmp = vtmp/nbgrp
              !
              IF(  (idesc(LAX_DESC_ACTIVE_NODE) > 0) .AND. &
                   (ipr-1 == idesc(LAX_DESC_MYR)) .AND. (ipc-1 == idesc(LAX_DESC_MYC)) ) THEN
                 CALL mp_root_sum( vtmp(:,1:nc), dm(:,icc:icc+nc-1), root, ortho_parent_comm )
              ELSE
                 CALL mp_root_sum( vtmp(:,1:nc), dm, root, ortho_parent_comm )
              END IF

           END DO
           !
        END IF
        !
     END DO
     !
     CALL laxlib_zsqmher( nbase+notcnv, dm, nx, idesc )
     !
     DEALLOCATE( vtmp )
     RETURN
  END SUBROUTINE update_distmat
  !
  !
  SUBROUTINE set_e_from_h()
     INTEGER :: nc, ic, i
     e(1:nbase) = 0_DP
     IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) .AND. la_proc ) THEN
        nc = idesc(LAX_DESC_NC)
        ic = idesc(LAX_DESC_IC)
        DO i = 1, nc
           e( i + ic - 1 ) = REAL( hl( i, i ) )
        END DO
     END IF
     CALL mp_sum( e(1:nbase), ortho_parent_comm )
     RETURN
  END SUBROUTINE set_e_from_h
  !
  SUBROUTINE set_h_from_e()
     INTEGER :: nc, ic, i
     IF( la_proc ) THEN
        hl = ZERO
        IF( idesc(LAX_DESC_MYC) == idesc(LAX_DESC_MYR) ) THEN
           nc = idesc(LAX_DESC_NC)
           ic = idesc(LAX_DESC_IC)
           DO i = 1, nc
              hl(i,i) = CMPLX( e( i + ic - 1 ), 0_DP ,kind=DP)
           END DO
        END IF
     END IF
     RETURN
  END SUBROUTINE set_h_from_e
  !
END SUBROUTINE pcegterg

END MODULE cegterg_mod
