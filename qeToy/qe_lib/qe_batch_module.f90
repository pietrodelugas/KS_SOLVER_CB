!
! qe_batch_module: device-resident batch arrays of qeToy (one column per batch = k-point of the group being solved),
! the FFT descriptor built from the dump and the loader of a group of k-points.
!
MODULE qe_batch_module
  USE qe_ham_module
  USE fft_types, ONLY : fft_type_descriptor
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  TYPE(fft_type_descriptor) :: dfft                       ! only the grid dimensions are set (full 3-D FFT, 'Rho' kind)
  INTEGER :: nb_group = 0                                 ! batches (k-points) per group
  INTEGER :: fft_cols_eff = 1                             ! columns of one batched FFT block
  INTEGER, ALLOCATABLE :: igk_batched(:,:), spin_batched(:)
  REAL(dp), ALLOCATABLE :: g2kin_batched(:,:), hdiag_batched(:,:), sdiag_batched(:,:)
  COMPLEX(dp), ALLOCATABLE :: vkb_batched(:,:,:), fft_array_batched(:,:)
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
  SUBROUTINE qe_alloc_batch(nb, cols)
    INTEGER, INTENT(IN) :: nb, cols
    nb_group = nb; fft_cols_eff = cols
    ALLOCATE(igk_batched(npwx,nb), spin_batched(nb), g2kin_batched(npwx,nb), hdiag_batched(npwx,nb), sdiag_batched(npwx,nb))
    ALLOCATE(vkb_batched(npwx,MAX(nkb,1),nb), fft_array_batched(nnr,cols))
    igk_batched = 1; spin_batched = 1; g2kin_batched = 0; hdiag_batched = 1; sdiag_batched = 1; vkb_batched = 0
    !$acc enter data create(igk_batched, spin_batched, g2kin_batched, hdiag_batched, sdiag_batched, vkb_batched, fft_array_batched) &
    !$acc copyin(nl, vrs, deeq, qq_at, ityp, nh, ofsbeta)
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
       CALL qe_read_k(ik, igk_batched(:,ib), g2kin_batched(:,ib), hdiag_batched(:,ib), sdiag_batched(:,ib), vkb_batched(:,:,ib), et)
       spin_batched(ib) = isk(ik)
    END DO
    !$acc update device(igk_batched, spin_batched, g2kin_batched, hdiag_batched, sdiag_batched, vkb_batched)
  END SUBROUTINE qe_load_group
END MODULE qe_batch_module
