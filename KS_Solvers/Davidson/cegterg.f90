
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
! M.Iovine : added the argument n_k as INTENT(IN)
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE cegterg( h_psi_ptr, s_psi_ptr, uspp, g_psi_ptr, &
                    npw, npwx, nvec, nvecx, npol, evc, ethr, &
                    e, btype, notcnv, lrot, dav_iter, nhpsi, i_batch, n_k )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
#if defined(__CUDA)
  use cublas
  use cudafor
  USE ieee_arithmetic, ONLY: ieee_is_nan, ieee_is_finite !! M.Iovine - added for debugging
#endif
  USE util_param,    ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id,&
                            nbgrp, my_bgrp_id, me_bgrp, root_bgrp
  USE mp,            ONLY : mp_sum, mp_gather, mp_bcast, mp_size,&
                            mp_type_create_column_section, mp_type_free
  USE device_memcpy_m, ONLY : dev_memcpy, dev_memset, dev_memcpy_async, &
                              dev_memset_async !Fixed: import device memory management routines
  USE mytime,          ONLY : clock_thread, clock_cuda_stream, cegterg_locker !Fixed: import thread private varibale
  USE openacc,         ONLY : acc_get_cuda_stream
#if defined(_OPENMP) 
  USE omp_lib, only:  omp_set_lock, omp_unset_lock
#endif
  !
  IMPLICIT NONE
  !
  include 'laxlib.fh'
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! number of spin polarizations
  INTEGER, INTENT(IN) :: i_batch
    ! batch index for k-point batching
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx*npol,nvec)
    !  evc contains the  refined estimates of the eigenvectors  
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  INTEGER, INTENT(OUT) :: nhpsi
    ! total number of individual hpsi
  INTEGER, INTENT(IN) :: n_k !! M.Iovine : added the number of k-points as argument in input for the allocation of arrays
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
    ! adapted npw and npwx
    ! do-loop counters
  INTEGER :: n_start, n_end, my_n
  INTEGER :: column_section_type
    ! defines a column section for communication
  INTEGER :: ierr
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
  !$acc declare device_resident(hc, sc, vc)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: hc_batched(:,:,:), sc_batched(:,:,:), vc_batched(:,:,:) !! M.Iovine - we define the corresponding 3D arrays for batched calls directly on the device
  !$acc declare device_resident(hc_batched, sc_batched, vc_batched)
  COMPLEX(DP), ALLOCATABLE :: hc_comp(:,:,:), sc_comp(:,:,:), vc_comp(:,:,:) !! M.Iovine - we define the corresponding 3D arrays
  !$acc declare device_resident(hc_comp, sc_comp, vc_comp)
  COMPLEX(DP), ALLOCATABLE :: hc_batched_old(:,:,:), sc_batched_old(:,:,:) !! M.Iovine added for debugging
  !$acc declare device_resident(hc_batched_old, sc_batched_old )
  REAL(DP), ALLOCATABLE :: ew_batched(:,:) !! M.Iovine - eigenvalues of the reduced Hamiltonian considered for batch on k-points directly on the device
  !$acc declare device_resident(ew_batched)
  REAL(DP), ALLOCATABLE, DEVICE :: ew_comp(:,:) !! M.Iovine - eigenvalues of the reduced Hamiltonian
  !$acc declare device_resident(ew_comp)
  REAL(DP), ALLOCATABLE :: ew_full_check(:) !! Added for debugging!!
  LOGICAL, ALLOCATABLE :: ew_full_nan(:), ew_full_inf(:) !! Added for debugging!!
  COMPLEX(DP), ALLOCATABLE :: vc_full_check(:,:,:) !! Added for debugging!!
  LOGICAL, ALLOCATABLE :: vc_full_nan(:,:), vc_full_inf(:,:) !! Added for debugging!!
  INTEGER :: first_nan_idx, ii_check !! M.IOvine - Added for debugging!!
  LOGICAL :: psi_has_nan !! M.IOvine - Added for debugging!!
  INTEGER :: ig, ibnd !! M.IOvine - Added for debugging!!
  COMPLEX(DP) :: maxv_hc, maxv_sc, maxv_vc !! M.Iovine - we define the variables for the maximum elements of the submatrix of dim. nbase x nbase 
  INTEGER, ALLOCATABLE :: nbase_batched(:) !! M.Iovine - definition of array to store the nbase value of each k_point -> of each thread
  INTEGER :: nbase_max, l !! M.Iovine - variable to store the maximum nbase found among the k-points and index l for the padding
  INTEGER :: info !!! M.Iovine - introduced for debugging
  !REAL(DP), ALLOCATABLE :: ew_host(:) !!! M.Iovine - introduced for debugging
  !COMPLEX(DP), ALLOCATABLE :: vc_host(:,:) !!! M.Iovine - introduced for debugging
  COMPLEX(DP), ALLOCATABLE :: vc_host(:,:) !!! M.Iovine - introduced for debugging
  REAL(DP), ALLOCATABLE :: ew_host(:) !!! M.Iovine - introduced for debugging
  COMPLEX(DP), ALLOCATABLE :: sc_host(:,:) !!! M.Iovine - introduced for debugging
  COMPLEX(DP), ALLOCATABLE :: hc_host(:,:) !!! M.Iovine - introduced for debugging
  REAL(DP), ALLOCATABLE :: e_check(:) !!! M.Iovine - introduced for debugging
  COMPLEX(DP), ALLOCATABLE :: hc_entry_check(:,:) !!! M.Iovine - degugging line
  LOGICAL, ALLOCATABLE :: hc_entry_nan(:,:) !!! M.Iovine - debugging line
  COMPLEX(DP), ALLOCATABLE :: hc_pre_check(:,:) !!! M.Iovine - debugging line
  LOGICAL :: has_nan_v, has_inf_v !!! M.Iovine - debugging line
  REAL(DP) :: max_abs_v !!! M.Iovine - debugging line
  COMPLEX(DP), ALLOCATABLE :: vc_check(:,:) !!M.Iovine - debugging line
  COMPLEX(DP), ALLOCATABLE :: ew_check(:) !!M.Iovine - debugging line
  REAL(DP) :: ef_check(nvec) !!M.Iovine - debugging line
  INTEGER :: istat !! M.Iovine - added for synchronization
  COMPLEX(DP), DEVICE, ALLOCATABLE :: d_res_chk(:), d_Sv_chk(:) !! M.Iovine - added for synchronization
  COMPLEX(DP), ALLOCATABLE         :: h_res_chk(:), h_Sv_chk(:) !! M.Iovine - added for synchronization
  REAL(DP)                         :: e1_chk, ortho_val_chk !! M.Iovine - added for synchronization
  INTEGER                          :: i_chk, j_chk, n_dim_chk !! M.Iovine - added for synchronization
  COMPLEX(DP), ALLOCATABLE :: Sv_chk(:) !! M.Iovine - added for synchronization
  COMPLEX(DP), ALLOCATABLE :: res_chk(:) !! M.Iovine - added for synchronization
  COMPLEX(DP), ALLOCATABLE :: Hv_chk(:) !! M.Iovine - added for synchronization
  REAL(DP) :: res_norm_chk !! M.Iovine - added for synchronization
  REAL(DP) :: s_norm_chk !! M.Iovine - added for synchronization
  REAL(DP) :: lambda_chk !! M.Iovine - added for synchronization
  REAL(DP) :: max_res_2   = 0.0_DP !! M.Iovine - added for synchronization
  REAL(DP) :: max_res_inf = 0.0_DP !! M.Iovine - added for synchronization
  REAL(DP) :: s_norm !! M.Iovine - added for synchronization
  REAL(DP) :: max_snorm_err = 0.0_DP !! M.Iovine - added for synchronization
  REAL(DP) :: res_norm_2 !!M.Iovine - added for synchronization
  REAL(DP) :: res_norm_inf !!M.Iovine - added for synchronization
  REAL(DP) :: max_ortho_err = 0.0_DP !!M.Iovine - added for synchronization
  REAL(DP) :: ortho_err !!M.Iovine - added for synchronization
  REAL(DP), ALLOCATABLE :: ew(:)
  !!$acc declare device_resident(ew)
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  INTEGER, ALLOCATABLE :: recv_counts(:), displs(:)
    ! receive counts and memory offsets
  INTEGER, PARAMETER :: blocksize = 256
  INTEGER :: numblock
    ! chunking parameters
  INTEGER :: i,j,k
  ! GPU stream management
  INTEGER :: async_id
  INTEGER(kind=cuda_stream_kind) :: mycudaStream, prova
#if defined(__CUDA)
  type(cublasHandle) :: myblasHandle(20) 
  INTEGER :: istat_cublas
#endif

  !
  REAL(DP), EXTERNAL :: MYDDOT_VECTOR_GPU
  !$acc routine(MYDDOT_VECTOR_GPU) vector
  !
  EXTERNAL  h_psi_ptr,    s_psi_ptr,    g_psi_ptr
    ! h_psi_ptr(npwx,npw,nvec,psi,hpsi,i_batch)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nvec,spsi,i_batch)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx*npol,nvec)
    ! g_psi_ptr(npwx,npw,notcnv,psi,e,i_batch)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  nhpsi = 0
  CALL start_clock( 'cegterg' ); !write(*,*) 'start cegterg' ; FLUSH(6)
  !
  ! Setup GPU stream using clock_thread (Fixed)
  async_id = clock_thread

#if defined(__CUDA)
  ! Get cuda stream and link cublas to it
  mycudaStream = clock_cuda_stream
  istat_cublas = cublasCreate(myblasHandle(i_batch))
  istat_cublas = cublasSetStream(myblasHandle(i_batch), mycudaStream)
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
#if ! defined(__CUDA)
  ! compute the number of chunks
  numblock  = (npw+blocksize-1)/blocksize
#endif
  !
  ALLOCATE(  psi( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx*npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !$acc enter data async(async_id) create(psi, hpsi)
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx*npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi ', ABS(ierr) )
     !$acc enter data async(async_id) create(spsi)
  END IF
  !
  ALLOCATE( sc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate sc ', ABS(ierr) )
  ALLOCATE( hc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hc ', ABS(ierr) )
  ALLOCATE( vc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc ', ABS(ierr) )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  !$acc enter data create(ew) async(async_id) 
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew ', ABS(ierr) )
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate conv ', ABS(ierr) )
  ALLOCATE( recv_counts(mp_size(inter_bgrp_comm)), displs(mp_size(inter_bgrp_comm)) )
  !!M.Iovine - allocation on the device of 3D arrays for batching :
  ALLOCATE(hc_batched(nvecx, nvecx, n_k))
  ALLOCATE(sc_batched(nvecx, nvecx, n_k))
  ALLOCATE(vc_batched(nvecx, nvecx, n_k))
  ALLOCATE(hc_batched_old(nvecx, nvecx, n_k)) !!M.Iovine - debugging
  ALLOCATE(sc_batched_old(nvecx, nvecx, n_k)) !!M.Iovine - debugging
  ALLOCATE(ew_batched(nvecx, n_k))
  ALLOCATE(nbase_batched(n_k))
  !ALLOCATE(ew_host(n_k)) !!! introduced for debugging
  !ALLOCATE(vc_host(nvecx,n_k)) !!! introduced for debugging
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  
  !!!DEBUGG LINES:
  ALLOCATE(hc_host(nbase, nbase))
  ALLOCATE(vc_host(nbase,nbase))
  ALLOCATE(ew_host(nbase))
  ALLOCATE(sc_host(nbase,nbase))
  ALLOCATE(Sv_chk(nbase))
  ALLOCATE(Hv_chk(nbase))
  ALLOCATE(res_chk(nbase))
  !!
  
  !
  !$acc host_data use_device(evc, psi)
  CALL dev_memcpy_async(psi, evc, mycudaStream, (/ 1 , npwx*npol /), 1, &
                            (/ 1 , nvec /), 1)
  !$acc end host_data
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi_ptr( npwx, npw, nvec, psi, hpsi, i_batch ) ; nhpsi = nhpsi + nvec
  !
  ! ... spsi contains s times the basis vectors
  !
  IF ( uspp ) CALL s_psi_ptr( npwx, npw, nvec, psi, spsi, i_batch )
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
  CALL start_clock( 'cegterg:init' )
  !
  !$acc host_data use_device(evc, psi, hpsi, spsi, hc, sc)
  CALL divide_all(inter_bgrp_comm,nbase,n_start,n_end,recv_counts,displs)
  CALL mp_type_create_column_section(sc(1,1), 0, nbase, nvecx, column_section_type)
  my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
  !
  if (n_start .le. n_end) &
  CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi, kdmx, hpsi(1,n_start), kdmx, ZERO, hc(1,n_start), nvecx )
  !
  if (n_start .le. n_end) & 
#if defined(__CUDA)
        CALL mp_sum( hc, 1, nbase, n_start, n_end , intra_bgrp_comm )
#else
        CALL mp_sum( hc(1:nbase, n_start:n_end), intra_bgrp_comm )
#endif
  CALL mp_gather( hc, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
  !
  IF ( uspp ) THEN
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi, kdmx, spsi(1,n_start), kdmx, &
                 ZERO, sc(1,n_start), nvecx )
     !
  ELSE
     !
     if (n_start .le. n_end) &
     CALL ZGEMM( 'C','N', nbase, my_n, kdim, ONE, psi, kdmx, psi(1,n_start), kdmx, &
                 ZERO, sc(1,n_start), nvecx )
     !
  END IF
  !
  if (n_start .le. n_end) & 
#if defined(__CUDA)
         CALL mp_sum( sc, 1, nbase, n_start, n_end , intra_bgrp_comm )
#else
         CALL mp_sum( sc(1:nbase, n_start:n_end), intra_bgrp_comm )
#endif
  CALL mp_gather( sc, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
  !$acc end host_data
  !
  CALL mp_type_free( column_section_type )
  !
  !$acc parallel vector_length(64) async(async_id) 
  !$acc loop gang 
  DO n = 1, nbase
     !
     ! ... the diagonal of hc and sc must be strictly real
     !
     hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 ,kind=DP)
     sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 ,kind=DP)
     !
     !$acc loop vector 
     DO m = n + 1, nbase
        !
        hc(n,m) = CONJG( hc(m,n) )
        sc(n,m) = CONJG( sc(m,n) )
        !
     END DO
     !
  END DO
  !$acc end parallel
  !
  CALL stop_clock( 'cegterg:init' )
  !
  IF ( lrot ) THEN
     !
     !$acc host_data use_device(vc)
     CALL dev_memset_async(vc, ZERO, mycudaStream, (/1, nbase/), 1, (/1, nbase/), 1)
     !$acc end host_data
     !
     !$acc parallel loop async(async_id)
     DO n = 1, nbase
        !
        e(n) = REAL( hc(n,n) )
        !
        vc(n,n) = ONE
        !
     END DO
     !
     CALL mp_bcast( e, root_bgrp_id, inter_bgrp_comm )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
     !!$acc host_data use_device(hc, sc, vc, ew, hc_batched, sc_batched, vc_batched, ew_batched)
     !CALL start_clock( 'cegterg:diag' )
     !$acc wait(async_id) !! M.Iovine - we add a wait before the assignment to the 3d arrays
     !call omp_set_lock(cegterg_locker) !!!! M.Iovine - We comment the line
     !!!! M.IOvine - Changes for 3D arrays: 
     nbase_batched(i_batch) = nbase !!!M.Iovine : we store the dimension for the current thread
     !$omp barrier
     !!! M.Iovine : we find the maximum for nbase_batched and nvecx_batched :
     nbase_max = maxval(nbase_batched)
     
     

     !$acc parallel loop collapse(2) present(hc, sc, hc_batched, sc_batched, vc_batched)
     do i = 1, nvecx
        do j = 1, nvecx
                hc_batched(j,i,i_batch) = CMPLX(hc(j,i), kind=DP)
                sc_batched(j,i,i_batch) = CMPLX(sc(j,i), kind=DP)
                vc_batched(j,i,i_batch) = CMPLX(vc(j,i), kind=DP)
        end do
     end do
     

 !!$acc parallel present(ew, ew_batched)
 !do j = 1, nvecx
  !    ew_batched(j,i_batch) = ew(j)
 !end do
      
     !!!M.Iovine - for each thread, we write the corresponding matricesin the batched arrays and we add the Padding:
     !if (nbase < nbase_max) then
       !!! We find the maximum element of the nbasexnbase arrays:
       !maxv_hc = maxval(abs(hc_batched(1:nbase,1:nbase,i_batch)))
       !do l=(nbase+1), nbase_max
            !hc_batched(:,l:nbase_max,i_batch) = CMPLX(0.D0,0.D0,kind=DP)
            !hc_batched(l:nbase_max,:,i_batch) = CMPLX(0.D0,0.D0,kind=DP)
            !sc_batched(:,l:nbase_max,i_batch) = CMPLX(0.D0,0.D0,kind=DP)
            !sc_batched(l:nbase_max,:,i_batch) = CMPLX(0.D0,0.D0,kind=DP)
            !vc_batched(:,l:nbase_max,i_batch) = CMPLX(0.D0,0.D0,kind=DP)
            !vc_batched(l:nbase_max,:,i_batch) = CMPLX(0.D0,0.D0,kind=DP)
  !!! Diagonal elements (we base them on the maximum element      of the reduced matrices --> we multiply for 10^5 for faster convergen     ce
            !hc_batched(l,l,i_batch) = CMPLX(l*1e3,0.D0,kind=DP)
            !sc_batched(l,l,i_batch) = CMPLX(1.D0,0.D0,kind=DP)
            !vc_batched(l,l,i_batch) = CMPLX(1.D0,0.D0,kind=DP)
       !end do
       !print *, 'I am inside the Padding loop'
     !end if

     !!$omp master
     
     !!!M.Iovine - we copy part of the 3d arrays into the computation arrays of dimension : nbase x nbase x n_k
     IF (.NOT. ALLOCATED(hc_comp)) ALLOCATE(hc_comp(nbase_max,nbase_max,n_k))
     IF (.NOT. ALLOCATED(sc_comp)) ALLOCATE(sc_comp(nbase_max,nbase_max,n_k))
     IF (.NOT. ALLOCATED(vc_comp)) ALLOCATE(vc_comp(nbase_max,nbase_max,n_k))
     IF (.NOT. ALLOCATED(ew_comp)) ALLOCATE(ew_comp(nbase_max,n_k))
     !!$acc enter data create(hc_comp, sc_comp, vc_comp, ew_comp)

     !$acc parallel loop collapse(2) present(hc, sc, vc, hc_comp, sc_comp, vc_comp)
     do i = 1, nbase_max
         do j = 1, nbase_max
                 hc_comp(j,i,i_batch) = hc(j,i)
                 sc_comp(j,i,i_batch) = sc(j,i)
                 vc_comp(j,i,i_batch) = vc(j,i)
         end do
     end do


     CALL start_clock( 'cegterg:diag' )
     !$omp master
     !$acc host_data use_device(hc_comp, sc_comp, ew_comp, vc_comp)
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase_max, nvec, hc_comp, sc_comp, nbase_max, ew_comp, vc_comp, n_k, me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     !$acc end host_data

     !$omp end master
     !$omp barrier 
     
     !$acc wait

     !$acc parallel loop present(ew)
     DO i = 1, nbase_max
          ew(i) = ew_comp(i, i_batch)
     END DO
     
     !$acc update self(ew)
     
     print *, 'Debugging', ew


     !$acc parallel loop collapse(2) present(vc)
     DO j = 1, nbase_max
          DO i = 1, nbase_max
             vc(i, j) = vc_comp(i, j, i_batch)
          END DO
     END DO
     
     !$acc update self(vc)
     print *, 'Debugging Vc', vc(:,1)
     
     !$acc wait(async_id)
     !!!!M.Iovine - deallocation of new arrays: 
     !!$acc exit data delete(hc_comp, vc_comp, ew_comp, sc_comp)
     DEALLOCATE(hc_comp)
     DEALLOCATE(vc_comp)
     DEALLOCATE(sc_comp)
     DEALLOCATE(ew_comp)
     !!!!
!!!! M.Iovine : copy to host of ew, vc, hc, sc for residual norm check :
     !$acc parallel loop present(ew) copyout(ew_host)
     DO i = 1, nvec
          ew_host(i) = ew(i)
     END DO
 
     !$acc parallel loop collapse(2) present(vc) copyout(vc_host)
     DO i = 1, nvec
       DO j = 1, nvec
        vc_host(j,i) = vc(j,i)
       END DO
     END DO
     
     print *, 'DEbugging vc host', vc_host(:,1)

     !$acc parallel loop collapse(2) present(hc) copyout(hc_host)
     DO i = 1, nbase
       DO j = 1, nbase
          hc_host(j,i) = hc(j,i)
       END DO
     END DO
     
     !$acc parallel loop collapse(2) present(sc) copyout(sc_host)
     DO i = 1, nbase
        DO j = 1, nbase
          sc_host(j,i) = sc(j,i)
        END DO
     END DO
!!!! M.Iovine - Residual check :
!DO k = 1, nbase
    !Hv_chk(:) = (0.0_DP, 0.0_DP)
    !Sv_chk(:) = (0.0_DP, 0.0_DP)

   ! DO i = 1, nbase
        !DO j = 1, nbase

        !    Hv_chk(i) = Hv_chk(i) + &
         !               hc_host(i,j) * vc_host(j,k)

       !     Sv_chk(i) = Sv_chk(i) + &
      !                  sc_host(i,j) * vc_host(j,k)

     !   END DO
    !END DO

    ! r_k = H v_k - lambda_k S v_k
    !DO i = 1, nbase
     !   res_chk(i) = Hv_chk(i) - &
      !               ew_host(k) * Sv_chk(i)
    !END DO

    ! ||r||_2
    !res_norm_2 = SQRT(SUM(ABS(res_chk(:))**2))

    ! ||r||_inf
    !res_norm_inf = MAXVAL(ABS(res_chk(:)))

    !max_res_2   = MAX(max_res_2,   res_norm_2)
    !max_res_inf = MAX(max_res_inf, res_norm_inf)

    !PRINT '(A,I6,A,ES16.8,A,ES16.8)', &
   !       'Eigenpair ', k, &
  !        '  residual L2 = ', res_norm_2, &
 !         '  residual Linf = ', res_norm_inf

! END DO 

!!!! M.Iovine - S-normlization:
!DO k = 1, nbase

 !   Sv_chk(:) = (0.0_DP, 0.0_DP)

  !  DO i = 1, nbase
   !     DO j = 1, nbase
    !        Sv_chk(i) = Sv_chk(i) + &
     !                   sc_host(i,j) * vc_host(j,k)
      !  END DO
    !END DO

    !s_norm = REAL( SUM(CONJG(vc_host(:,k)) * Sv_chk(:)), KIND=DP )

    !max_snorm_err = MAX(max_snorm_err, ABS(s_norm - 1.0_DP))

   ! PRINT '(A,I6,A,ES16.8,A,ES16.8)', &
    !      'Eigenvector ', k, &
     !     '  S-norm = ', s_norm, &
      !    '  error = ', ABS(s_norm - 1.0_DP)

!END DO

!!!!

!!!! M.Iovine S-orthogonality:
!DO k = 1, nbase

 ! DO l = 1, nbase

      ! Compute v_k^H S v_l

  !    s_norm = 0.0_DP

   !   DO i = 1, nbase

    !       Sv_chk(i) = (0.0_DP, 0.0_DP)

     !      DO j = 1, nbase
      !          Sv_chk(i) = Sv_chk(i) + &
       !                     sc_host(i,j) * vc_host(j,l)
        !   END DO

         !  s_norm = s_norm + &
          !          REAL(CONJG(vc_host(i,k)) * Sv_chk(i), KIND=DP)

      !END DO

      !IF (k == l) THEN
       ! ortho_err = ABS(s_norm - 1.0_DP)
      !ELSE
       ! ortho_err = ABS(s_norm)
      !END IF

     ! max_ortho_err = MAX(max_ortho_err, ortho_err)

    !END DO

!END DO

!PRINT *, 'Maximum S-orthogonality error = ', max_ortho_err

!!!!
     IF( nbgrp > 1 ) THEN
        !$acc host_data use_device(vc, ew)
        CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
        !$acc end host_data
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !
     CALL dev_memcpy_async(e, ew, mycudaStream, (/ 1, nvec /), 1 )
     !


!!!! : M.Iovine - debugging on vc:
IF (.NOT. ALLOCATED(vc_check)) ALLOCATE(vc_check(nbase,nbase))   !! size to whatever region you actually want to check
!$acc parallel loop collapse(2) present(vc) copyout(vc_check)
do i = 1, nbase
   do j = 1, nbase
      vc_check(j,i) = vc(j,i)
   end do
end do

!!! Debugging - M.Iovine - vc print in a file:
   open(unit=10, file='vc_print_new.dat', status ='replace', action='write')
     do i = 1, nbase_max
           write (10, '(10("(",F10.6,",",F10.6,") "))') (vc_check(i,j), j = 1, nbase_max)
     end do
   close(10)

!!!!

!!!! : M.Iovine - debugging on ew:
IF (.NOT. ALLOCATED(ew_check)) ALLOCATE(ew_check(nbase))   !! size to whatever region you actually want to check
!$acc parallel loop  present(ew) copyout(ew_check)
do i = 1, nbase
      ew_check(i) = ew(i)
end do

!!! Debugging - M.Iovine - vc print in a file:
open(unit=10, file='ew_print_new.dat', status ='replace', action='write')
do i = 1, nbase
        write(10, '(F6.2)') ew_check(i)
end do
close(10)

   !!! Debugging - M.Iovine - e print in a file:
   !!!! : M.Iovine - debugging on ew:
    !$acc parallel loop present(e) copyout(e_check)
    do i = 1, nvec
          ef_check(i) = e(i)
    end do

    open(unit=10, file='e_print_new.dat', status ='replace', action='write')
     do i = 1, nvec
             write(10, '(F6.2)') ef_check(i)
     end do
    close(10)
!!!!!!!


has_nan_v = ANY(IEEE_IS_NAN(REAL(vc_check))) .OR. ANY(IEEE_IS_NAN(AIMAG(vc_check)))
has_inf_v = .NOT. ALL(IEEE_IS_FINITE(REAL(vc_check))) .OR. .NOT. ALL(IEEE_IS_FINITE(AIMAG(vc_check)))
max_abs_v = MAXVAL(ABS(vc_check))

        IF (has_nan_v) THEN
               CALL lax_error__( 'cdiaghg', 'CRITICAL: NaN detected in vc / v_d!', 1 )
        ELSE IF (has_inf_v) THEN
               CALL lax_error__( 'cdiaghg', 'CRITICAL: Inf detected in vc / v_d!', 1 )
        ELSE IF (max_abs_v < 1.0d-12) THEN
               WRITE(*,*) 'WARNING: vc / v_d appears to be completely ZERO or uninitialized!'
        ELSE
               WRITE(*,'(A,E16.8)') ' [VC CHECK] v_d input valid. maxabs = ', max_abs_v
        END IF

END IF
 !!!!!!  

  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter ; !write(*,*) kter, notcnv, conv
     !
     CALL start_clock( 'cegterg:update' )
     !
     np = 0
     !
     DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
           !
           ! ... this root not yet converged ... 
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix 
           ! ... multiplications to set a new basis vector (see below)
           !
           IF ( np /= n ) THEN
             !$acc parallel loop async(async_id)  
             DO i = 1, nvecx
               vc(i,np) = vc(i,n)
             END DO 
           END IF 
           !
           ! ... for use in g_psi_ptr
           !
           !$acc kernels async(async_id) 
           ew(nbase+np) = e(n)
           !$acc end kernels 
           !
        END IF
        !
     END DO
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
     !
     !$acc host_data use_device(psi, spsi, vc)
     IF ( uspp ) THEN
        !
        if (n_start .le. n_end) &
        CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, spsi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                    ZERO, psi(1,nb1), kdmx )
        !     
     ELSE
        !
        if (n_start .le. n_end) &
        CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, psi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                    ZERO, psi(1,nb1), kdmx )
        !
     END IF
     !$acc end host_data
! NB: must not call mp_sum over inter_bgrp_comm here because it is done later to the full correction
     !
#if defined(__CUDA)
     !$acc parallel loop collapse(3) async(async_id) 
     DO np=1,notcnv
        DO ipol = 1, npol
           DO k=1,npwx
             psi(k + (ipol-1)*npwx,nbase+np) = - ew(nbase+np)*psi(k + (ipol-1)*npwx,nbase+np)
           END DO
        END DO
     END DO
     !$acc end parallel 
#else
     !$omp parallel do collapse(3)
     DO n = 1, notcnv
        DO ipol = 1, npol
           DO m = 1, numblock
              psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw, m*blocksize)+(ipol-1)*npwx,nbase+n) = &
                         - ew(nbase+n) * &
              psi( (m-1)*blocksize+(ipol-1)*npwx+1: &
                    MIN(npw, m*blocksize)+(ipol-1)*npwx,nbase+n)
           END DO
        END DO
     END DO
     !$omp end parallel do
#endif
     !
     !$acc host_data use_device(psi, hpsi, vc, ew)
     if (n_start .le. n_end) &
     CALL ZGEMM( 'N','N', kdim, notcnv, my_n, ONE, hpsi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                 ONE, psi(1,nb1), kdmx )
     CALL mp_sum( psi(:,nb1:nbase+notcnv), inter_bgrp_comm )
     !
     ! clean up garbage if there is any
     IF (npw < npwx) CALL dev_memset_async(psi, ZERO, mycudaStream, [npw+1,npwx], 1, [nb1, nbase+notcnv])
     IF (npol == 2)  CALL dev_memset_async(psi, ZERO, mycudaStream, [npwx+npw+1,2*npwx], 1, [nb1, nbase+notcnv])
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi_ptr( npwx, npw, notcnv, npol, psi(1,nb1), ew(nb1), i_batch )
     !$acc end host_data
     !
     !!!! M.Iovine - debugging check for NaN prints in output after preconditioning:
     ! --- Right after preconditioning (e.g., after CALL g_psi or psi update) ---
     !psi_has_nan = .FALSE.

     !!$acc parallel loop collapse(2) reduction(.OR. : psi_has_nan) present(psi)
     !DO ibnd = 1, notcnv
      !  DO ig = 1, npw
       !         IF (IEEE_IS_NAN(REAL(psi(ig, ibnd))) .OR. IEEE_IS_NAN(AIMAG(psi(ig, ibnd)))) THEN
        !           psi_has_nan = .TRUE.
         !       END IF
       ! END DO
     !END DO

    ! IF (psi_has_nan) THEN
     !   WRITE(6, '(A, I0)') 'CRITICAL: NaN detected in psi after preconditioning at i_batch = ', i_batch
      !  FLUSH(6)
       ! STOP 'cegterg: psi infected with NaN'
    ! END IF
     !!!!
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     !$acc parallel vector_length(96) async(async_id)  
     !$acc loop gang private(nbn)
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        ew(n) = MYDDOT_VECTOR_GPU( 2*npw, psi(1,nbn), psi(1,nbn) )
        !
     END DO
     !
     IF(npol.ne.1) THEN 
       !$acc loop gang private(nbn)
       DO n = 1, notcnv 
         nbn = nbase + n
         ew(n) = ew(n)  + MYDDOT_VECTOR_GPU( 2*npw, psi(npwx+1,nbn), psi(npwx+1,nbn) ) 
       END DO 
     END IF 
     !$acc end parallel 
     !
     !$acc host_data use_device(ew)
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !$acc end host_data
     !
#if defined(__CUDA)
     !$acc parallel loop collapse(3) async(async_id) 
     DO i = 1,notcnv
        DO ipol = 1,npol
           DO k=1,npw
             psi(k + (ipol-1)*npwx,nbase+i) = psi(k+(ipol-1)*npwx,nbase+i)/SQRT( ew(i) )
           END DO
        END DO
     END DO
#else
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
#endif
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi_ptr( npwx, npw, notcnv, psi(1,nb1), hpsi(1,nb1), i_batch ) ; nhpsi = nhpsi + notcnv
     !
     IF ( uspp ) CALL s_psi_ptr( npwx, npw, notcnv, psi(1,nb1), spsi(1,nb1), i_batch )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
     !$acc host_data use_device(psi, hpsi, spsi, hc, sc)
     CALL divide_all(inter_bgrp_comm,nbase+notcnv,n_start,n_end,recv_counts,displs)
     CALL mp_type_create_column_section(sc(1,1), nbase, notcnv, nvecx, column_section_type)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     !
     CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, hpsi(1,nb1), kdmx, psi(1,n_start), kdmx, &
                 ZERO, hc(nb1,n_start), nvecx )
     !
     if (n_start .le. n_end) &
#if defined(__CUDA)
       CALL mp_sum( hc, nb1, nbase+notcnv, n_start, n_end , intra_bgrp_comm )
#else
       CALL mp_sum( hc(nb1:nbase+notcnv, n_start:n_end) , intra_bgrp_comm )
#endif
     CALL mp_gather( hc, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
     !
     CALL divide(inter_bgrp_comm,nbase+notcnv,n_start,n_end)
     my_n = n_end - n_start + 1; !write (*,*) nbase+notcnv,n_start,n_end
     IF ( uspp ) THEN
        !
        CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, spsi(1,nb1), kdmx, psi(1,n_start), kdmx, &
                    ZERO, sc(nb1,n_start), nvecx )
        !     
     ELSE
        !
        CALL ZGEMM( 'C','N', notcnv, my_n, kdim, ONE, psi(1,nb1), kdmx, psi(1,n_start), kdmx, &
                    ZERO, sc(nb1,n_start), nvecx )
        !
     END IF
     !
     if (n_start .le. n_end) & 
#if defined(__CUDA)
         CALL mp_sum( sc, nb1, nbase+notcnv, n_start, n_end , intra_bgrp_comm )
#else
         CALL mp_sum( sc(nb1:nbase+notcnv, n_start:n_end) , intra_bgrp_comm )
#endif
     CALL mp_gather( sc, column_section_type, recv_counts, displs, root_bgrp_id, inter_bgrp_comm )
     !$acc end host_data
     !
     CALL mp_type_free( column_section_type )
     !
     CALL stop_clock( 'cegterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     !$acc parallel vector_length(64) async(async_id) 
     !$acc loop gang
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real
        !
        IF( n>=nb1 ) THEN
           hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 ,kind=DP)
           sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 ,kind=DP)
        ENDIF
        !
        !$acc loop vector
        DO m = MAX(n+1,nb1), nbase
           !
           hc(n,m) = CONJG( hc(m,n) )
           sc(n,m) = CONJG( sc(m,n) )
           !
        END DO
        !
     END DO
     !$acc end parallel 
     !
     ! ... diagonalize the reduced hamiltonian
     !
     call omp_set_lock(cegterg_locker) 
     !$acc host_data use_device(hc, sc, vc, ew)
     CALL start_clock( 'cegterg:diag' )
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL diaghg( nbase, nvec, hc, sc, nvecx, ew, vc, me_bgrp, root_bgrp, intra_bgrp_comm )
     END IF
     !$acc wait(async_id)
     call omp_unset_lock(cegterg_locker) 
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     CALL stop_clock( 'cegterg:diag' )
     !$acc end host_data
     !
     ! ... test for convergence
     !
     !$acc parallel loop copy(conv(1:nvec)) copyin(btype(1:nvec)) async(async_id)
     DO i = 1, nvec
       IF(btype(i) == 1) THEN
         conv(i) = ( ( ABS( ew(i) - e(i) ) < ethr ) )
       ELSE
         conv(i) = ( ( ABS( ew(i) - e(i) ) < empty_ethr ) )
       END IF 
     END DO
     !
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     !$acc host_data use_device(ew)
     CALL dev_memcpy_async (e, ew, mycudaStream, (/ 1, nvec /) )
     !$acc end host_data
     !
     
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. &
          nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'cegterg:last' )
        !
        CALL divide(inter_bgrp_comm,nbase,n_start,n_end)
        my_n = n_end - n_start + 1; !write (*,*) nbase,n_start,n_end
        !$acc host_data use_device(evc, psi, vc)
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, psi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                    ZERO, evc, kdmx )
        CALL mp_sum( evc, inter_bgrp_comm )
        !$acc end host_data
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
        !$acc host_data use_device(evc, psi, hpsi, spsi, vc)
        CALL dev_memcpy_async(psi, evc, mycudaStream, (/ 1, npwx*npol /), 1, &
                                      (/ 1, nvec /), 1)
        !
        IF ( uspp ) THEN
           !
           CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, spsi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                       ZERO, psi(1,nvec+1), kdmx)
           CALL dev_memcpy_async(spsi, psi(:,nvec+1:), mycudaStream,&
                                        (/1, npwx*npol/), 1, &
                                        (/1, nvec/), 1)
           CALL mp_sum( spsi(:,1:nvec), inter_bgrp_comm )
           !
        END IF
        !
        CALL ZGEMM( 'N','N', kdim, nvec, my_n, ONE, hpsi(1,n_start), kdmx, vc(n_start,1), nvecx, &
                    ZERO, psi(1,nvec+1), kdmx )
        CALL dev_memcpy_async(hpsi, psi(:,nvec+1:), mycudaStream,&
                                        (/1, npwx*npol/), 1, &
                                        (/1, nvec/), 1)
        CALL mp_sum( hpsi(:,1:nvec), inter_bgrp_comm )
        !$acc end host_data
        !
        ! ... refresh the reduced hamiltonian 
        !
        nbase = nvec
        !
        ! These variables are set to ZERO in the CUF Kernel below
        !hc(1:nbase,1:nbase) = ZERO
        !sc(1:nbase,1:nbase) = ZERO
        !vc(1:nbase,1:nbase) = ZERO
        !
        !$acc kernels  async(async_id) 
        DO n = 1, nbase
           hc(n,n) = CMPLX( e(n), 0.0_DP ,kind=DP)
           sc(n,n) = ONE
           vc(n,n) = ONE
           DO j = n+1, nbase
              hc(j,n) = ZERO
              hc(n,j) = ZERO
              sc(j,n) = ZERO
              sc(n,j) = ZERO
              vc(j,n) = ZERO
              vc(n,j) = ZERO
           END DO
           !
        END DO
        !$acc end kernels
        !
        CALL stop_clock( 'cegterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  !$acc exit data delete(ew) async(async_id) 
  DEALLOCATE( recv_counts )
  DEALLOCATE( displs )
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  DEALLOCATE( ew_batched ) !!M.Iovine - Deallocation of the arrays for the bacthed NVIDIA routines
  DEALLOCATE( vc_batched ) 
  DEALLOCATE( hc_batched ) 
  DEALLOCATE( sc_batched )
  !!! M.Iovine - added for debugging:
  IF (ALLOCATED(ew_full_check)) DEALLOCATE(ew_full_check, ew_full_nan, ew_full_inf)
  IF (ALLOCATED(vc_full_check)) DEALLOCATE(vc_full_check, vc_full_nan, vc_full_inf)
  DEALLOCATE(hc_host)
  DEALLOCATE(vc_host)
  DEALLOCATE(ew_host)
  DEALLOCATE(sc_host)
  DEALLOCATE(Sv_chk)
  DEALLOCATE(Hv_chk)
  DEALLOCATE(res_chk)
  DEALLOCATE(hc_batched_old) !!M.Iovine - debugging
  DEALLOCATE(sc_batched_old) !!M.Iovine - debugging
  DEALLOCATE(vc_check)
  DEALLOCATE(ew_check)
  !!!!!!!!!!
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
  ! Cleanup (Fixed)
#if defined(__CUDA)
  istat_cublas = cublasDestroy(myblasHandle(i_batch))  
#endif
  !
  CALL stop_clock( 'cegterg' ); !write(*,*) 'stop cegterg' ; FLUSH(6)
  !call print_clock( 'cegterg' )
  !call print_clock( 'cegterg:init' )
  !call print_clock( 'cegterg:diag' )
  !call print_clock( 'cegterg:update' )
  !call print_clock( 'cegterg:overlap' )
  !call print_clock( 'cegterg:last' )
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
