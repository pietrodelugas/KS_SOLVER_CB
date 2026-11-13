!
! qe_batch_module: device-resident batch arrays of qeToy (one column per batch = k-point of the group being solved),
! the FFT descriptor built from the dump and the loader of a group of k-points.
!
MODULE qe_batch_module
  USE qe_ham_module
  USE qe_kpoint_module, ONLY : qe_make_kpoint
  USE fft_types, ONLY : fft_type_descriptor
  USE openacc, ONLY : c_devptr
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  TYPE(fft_type_descriptor) :: dfft                       ! only the grid dimensions are set (full 3-D FFT, 'Rho' kind)
  INTEGER :: nb_group = 0                                 ! batches (k-points) per group
  INTEGER :: fft_cols_eff = 1                             ! columns of one batched FFT block
  INTEGER, ALLOCATABLE :: igk_batched(:,:), spin_batched(:)
  REAL(dp), ALLOCATABLE :: g2kin_batched(:,:), hdiag_batched(:,:), sdiag_batched(:,:)
  COMPLEX(dp), ALLOCATABLE :: vkb_batched(:,:,:), fft_array_batched(:,:)
  COMPLEX(dp), ALLOCATABLE :: becp_batched(:,:,:), ps_batched(:,:,:)   ! non local work arrays of H|psi> (indexed by slot)
  COMPLEX(dp), ALLOCATABLE :: becp_s(:,:,:), ps_s(:,:,:)               ! the same for S|psi> (indexed by batch/thread)
  ! S|psi> of the new vectors, computed together with the non local part of H|psi> (shared becp) and indexed by batch:
  COMPLEX(dp), ALLOCATABLE :: spsi_batch(:,:,:), ps2_batched(:,:,:)
  INTEGER, ALLOCATABLE :: sp_nb1(:), sp_m(:)          ! first column and number of columns of the valid S|psi> of every batch (host)
  LOGICAL, ALLOCATABLE :: sp_valid(:)                  ! S|psi> of batch ib is already computed (set by H|psi>, cleared by qe_s_psi_batched)
  TYPE(c_devptr), ALLOCATABLE :: ptr_ps2(:), ptr_sp_w(:)
  LOGICAL :: nonlocal_batched = .TRUE.                                 ! one batched GEMM pair for all the slots (else a GEMM pair per slot)
  TYPE(c_devptr), ALLOCATABLE :: ptr_vkb(:), ptr_psi_w(:), ptr_becp(:), ptr_ps(:), ptr_hpsi_w(:)   ! pointer arrays of the batched GEMMs (device)
  INTEGER, ALLOCATABLE :: kb2atom(:), kb2ih(:)                          ! atom and projector index within the atom of every vkb column
CONTAINS
  SUBROUTINE qe_setup_grid()
    ! The grid of the descriptor is the one of the dump (nothing to check against: it is the dump)
    dfft%nr1 = nr1; dfft%nr2 = nr2; dfft%nr3 = nr3
    dfft%nr1x = nr1x; dfft%nr2x = nr2x; dfft%nr3x = nr3x; dfft%nnr = nnr
    dfft%lpara = .FALSE.; dfft%lgamma = .FALSE.
    dfft%rho_clock_label = 'ffts'; dfft%wave_clock_label = 'fftw'
    IF (nr1x /= nr1 .OR. nr2x /= nr2 .OR. nr3x /= nr3) STOP 'qe_setup_grid: batched cuFFT needs nrNx = nrN'
    IF (nnr /= nr1*nr2*nr3) STOP 'qe_setup_grid: nnr /= nr1*nr2*nr3'
  END SUBROUTINE qe_setup_grid
  !
  SUBROUTINE qe_alloc_batch(nb, cols, nbndx)
    INTEGER, INTENT(IN) :: nb, cols, nbndx
    INTEGER :: na, ih
    nb_group = nb; fft_cols_eff = cols
    ALLOCATE(igk_batched(npwx,nb), spin_batched(nb), g2kin_batched(npwx,nb), hdiag_batched(npwx,nb), sdiag_batched(npwx,nb))
    ALLOCATE(vkb_batched(npwx,MAX(nkb,1),nb), fft_array_batched(nnr,cols))
    ALLOCATE(becp_batched(MAX(nkb,1),nbndx,nb), ps_batched(MAX(nkb,1),nbndx,nb), becp_s(MAX(nkb,1),nbndx,nb), ps_s(MAX(nkb,1),nbndx,nb))
    ALLOCATE(kb2atom(MAX(nkb,1)), kb2ih(MAX(nkb,1))); kb2atom = 1; kb2ih = 1
    IF (okvan) THEN
       ALLOCATE(spsi_batch(npwx,nbndx,nb), ps2_batched(MAX(nkb,1),nbndx,nb), ptr_ps2(nb), ptr_sp_w(nb))
       spsi_batch = 0; ps2_batched = 0
       !$acc enter data create(spsi_batch, ps2_batched, ptr_ps2, ptr_sp_w)
    END IF
    ALLOCATE(sp_nb1(nb), sp_m(nb), sp_valid(nb)); sp_nb1 = 1; sp_m = 0; sp_valid = .FALSE.
    ALLOCATE(ptr_vkb(nb), ptr_psi_w(nb), ptr_becp(nb), ptr_ps(nb), ptr_hpsi_w(nb))
    !$acc enter data create(ptr_vkb, ptr_psi_w, ptr_becp, ptr_ps, ptr_hpsi_w)
    DO na = 1, nat; DO ih = 1, nh(ityp(na))
       kb2atom(ofsbeta(na)+ih) = na; kb2ih(ofsbeta(na)+ih) = ih
    END DO; END DO
    igk_batched = 1; spin_batched = 1; g2kin_batched = 0; hdiag_batched = 1; sdiag_batched = 1; vkb_batched = 0
    !$acc enter data create(igk_batched, spin_batched, g2kin_batched, hdiag_batched, sdiag_batched, vkb_batched, fft_array_batched, &
    !$acc becp_batched, ps_batched, becp_s, ps_s) copyin(nl, vrs, deeq, qq_at, ityp, nh, ofsbeta, kb2atom, kb2ih)
  END SUBROUTINE qe_alloc_batch
  !
  SUBROUTINE qe_load_group(k0, n)
    ! k-points k0+1..k0+n of the run into batches 1..n (zero padded beyond npw), then to the device
    INTEGER, INTENT(IN) :: k0, n
    INTEGER :: ib, ik, npw
    REAL(dp), ALLOCATABLE :: et(:)
    ALLOCATE(et(nbnd))
    igk_batched = 1; g2kin_batched = 0; hdiag_batched = 1; sdiag_batched = 1; vkb_batched = 0
    DO ib = 1, n
       ik = k0 + ib; npw = ngk(ik)
       IF (kcompute) THEN
          CALL qe_make_kpoint(xk(:,ik), isk(ik), npwx, npw, igk_batched(:,ib), g2kin_batched(:,ib), hdiag_batched(:,ib), &
                              sdiag_batched(:,ib), vkb_batched(:,:,ib))
          IF (npw /= ngk(ik)) STOP 'qe_load_group: inconsistent number of plane waves'
       ELSE
          CALL qe_read_k(ik, igk_batched(:,ib), g2kin_batched(:,ib), hdiag_batched(:,ib), sdiag_batched(:,ib), vkb_batched(:,:,ib), et)
       END IF
       spin_batched(ib) = isk(ik)
    END DO
    !$acc update device(igk_batched, spin_batched, g2kin_batched, hdiag_batched, sdiag_batched, vkb_batched)
  END SUBROUTINE qe_load_group
  !
  SUBROUTINE qe_cublas_prepare(i_batch)
    ! cuBLAS handle of this thread, bound to the stream of its OpenACC queue (cegterg does the same at its start)
    USE laxlib_cusolver_handles, ONLY : cublas_handle, cublas_initialized
    USE mytime, ONLY : clock_cuda_stream
    USE cublas
    INTEGER, INTENT(IN) :: i_batch
    INTEGER :: info
    IF (.NOT. cublas_initialized(i_batch)) THEN
       info = cublasCreate(cublas_handle(i_batch))
       IF (info /= CUBLAS_STATUS_SUCCESS) STOP 'qe_cublas_prepare: cublasCreate'
       cublas_initialized(i_batch) = .TRUE.
       info = cublasSetStream(cublas_handle(i_batch), clock_cuda_stream)
       IF (info /= CUBLAS_STATUS_SUCCESS) STOP 'qe_cublas_prepare: cublasSetStream'
    END IF
  END SUBROUTINE qe_cublas_prepare
END MODULE qe_batch_module
