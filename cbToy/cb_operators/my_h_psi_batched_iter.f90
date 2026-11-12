subroutine my_h_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_comp, nb1_comp, psi, hpsi, i_batch, nk_batches, slot2batch)
  !! H*psi for the new correction vectors of all the k-points (batches) still active in the Davidson iterations.
  !!
  !! Slot ik (1..nk_batches) is the position of a batch among the active ones; slot2batch(ik) is the batch it belongs to
  !! (it selects igk_batched / ekin_batched). For every slot, nbnd_comp(ik) bands starting at column nb1_comp(ik) of
  !! psi(:,:,ik) are transformed.
  !!
  !! The FFT part (the local potential) works on chunks of bands and on all the slots at once:
  !!
  !!      column of fft_array_batched  col = jb + (ik-1)*chunk        jb = 1..chunk (band inside the chunk), ik = slot
  !!
  !! so that one chunk costs one zeroing kernel, one gather kernel, one batched invfft (howmany = chunk*nk_batches),
  !! one vloc kernel, one batched fwfft and one scatter kernel, with no host synchronization inside the loop.
  !! A term that needs the same (slot, band) columns (e.g. a non local part built from becp) can reuse the layout.
  !!
  !! chunk is a power of two <= band_chunk_eff, chosen for every chunk as the smallest one covering the bands still
  !! to be done, so that the same few cuFFT plans (howmany = chunk*n, created at start-up) are always enough.
  !! Columns of bands beyond nbnd_comp(ik) are left at zero.
  use iso_fortran_env, only: dp=> real64
  use cb_module, only: ekin_batched, vloc, igk_batched, dfft, fft_array_batched, band_chunk_eff
  use fft_interfaces, only: fwfft, invfft
  use mytime, only: clock_thread
  implicit none
  integer,intent(in) :: nk_batches
  integer,intent(in) :: npwx, i_batch, nbndx
  integer,intent(in) :: npw_comp(nk_batches)
  integer,intent(in) :: nbnd_comp(nk_batches) ! bands to be done for every slot (slot-indexed, host array)
  integer,intent(in) :: nb1_comp(nk_batches)  ! first column of the new vectors for every slot (slot-indexed, on the device)
  integer,intent(in) :: slot2batch(nk_batches)
  complex(dp) :: psi(npwx, nbndx, nk_batches), hpsi(npwx, nbndx, nk_batches)
  integer :: ig, ir, ik, jb, ib, ibr, col, ic
  integer :: b0, nb_max, nrem, chunk, ncols, nnr

  nnr    = dfft%nnr
  nb_max = maxval(nbnd_comp(1:nk_batches))

  !$acc data present(psi, hpsi, npw_comp, nb1_comp) copyin(nbnd_comp) async(clock_thread)
  b0 = 1
  do while (b0 <= nb_max)
     ! smallest power of two covering the remaining bands, not larger than band_chunk_eff
     nrem  = nb_max - b0 + 1
     chunk = 1
     do while (chunk < min(nrem, band_chunk_eff))
        chunk = 2*chunk
     end do
     ncols = chunk*nk_batches

     !$acc kernels async(clock_thread)
     fft_array_batched(:,1:ncols) = cmplx(0._dp, 0._dp, kind=dp)
     !$acc end kernels

     ! G-sphere -> FFT grid, all the (slot, band) columns of the chunk
     !$acc parallel loop collapse(3) async(clock_thread) private(ib, ibr, col)
     do ik = 1, nk_batches
        do jb = 1, chunk
           do ig = 1, npwx
              ib = b0 + jb - 1
              if (ig .le. npw_comp(ik) .and. ib .le. nbnd_comp(ik)) then
                 col = jb + (ik-1)*chunk
                 ibr = ib + nb1_comp(ik) - 1
                 fft_array_batched(dfft%nl(igk_batched(ig,slot2batch(ik))), col) = psi(ig,ibr,ik)
              end if
           end do
        end do
     end do

     !$acc host_data use_device(fft_array_batched)
     call invfft( 'Wave', fft_array_batched(:,:), dfft, ncols )
     !$acc end host_data

     !$acc parallel loop collapse(2) async(clock_thread)
     do ic = 1, ncols
        do ir = 1, nnr
           fft_array_batched(ir,ic) = fft_array_batched(ir,ic) * vloc(ir)
        end do
     end do

     !$acc host_data use_device(fft_array_batched)
     call fwfft( 'Wave', fft_array_batched(:,:), dfft, ncols )
     !$acc end host_data

     ! FFT grid -> G-sphere, plus the kinetic energy
     !$acc parallel loop collapse(3) async(clock_thread) private(ib, ibr, col)
     do ik = 1, nk_batches
        do jb = 1, chunk
           do ig = 1, npwx
              ib = b0 + jb - 1
              if (ig .le. npw_comp(ik) .and. ib .le. nbnd_comp(ik)) then
                 col = jb + (ik-1)*chunk
                 ibr = ib + nb1_comp(ik) - 1
                 hpsi(ig,ibr,ik) = ekin_batched(ig,slot2batch(ik)) * psi(ig,ibr,ik) +  &
                                   fft_array_batched(dfft%nl(igk_batched(ig,slot2batch(ik))), col)
              end if
           end do
        end do
     end do

     b0 = b0 + chunk
  end do
  !$acc wait(clock_thread)
  !$acc end data

end subroutine my_h_psi_batched_iter
