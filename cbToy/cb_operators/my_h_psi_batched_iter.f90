subroutine my_h_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_comp, nb1_comp, psi, hpsi, i_batch, nk_batches, slot2batch)
  !! H*psi for the new correction vectors of all the k-points (batches) still active in the Davidson iterations.
  !!
  !! Slot ik (1..nk_batches) is the position of a batch among the active ones; slot2batch(ik) is the batch it belongs to
  !! (it selects igk_batched / ekin_batched). For every slot, nbnd_comp(ik) bands starting at column nb1_comp(ik) of
  !! psi(:,:,ik) are transformed.
  !!
  !! The FFT part (the local potential) works on columns of fft_array_batched, one for every (slot, band) pair, packed
  !! one after the other (no padding, slots with fewer bands simply contribute fewer columns):
  !!
  !!      column ic  <->  (col_slot(ic), col_band(ic))
  !!
  !! The columns are processed in blocks of at most fft_cols_eff. For a block of nc columns: one zeroing kernel, one
  !! gather kernel, the batched invfft, one vloc kernel, the batched fwfft and one scatter kernel, with no host
  !! synchronization inside the loop (the FFT runs on the same stream as the OpenACC queue). The batched FFT of a block
  !! is split in pieces whose size is a power of two (binary expansion of nc, at most log2(fft_cols_eff)+1 calls), so
  !! that only the plans for 1, 2, 4, ... fft_cols_eff columns (created at start-up) are needed.
  !! A term that needs the same (slot, band) columns (e.g. a non local part built from becp) can reuse this layout.
  use iso_fortran_env, only: dp=> real64
  use cb_module, only: ekin_batched, vloc, igk_batched, dfft, fft_array_batched, fft_cols_eff
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
  integer :: col_slot(nk_batches*nbndx), col_band(nk_batches*nbndx)
  integer :: ig, ir, ik, ib, ibr, ic, jc
  integer :: ncol, c0, nc, pos, rem, k, nnr

  nnr = dfft%nnr

  ! list of the (slot, band) pairs
  ncol = 0
  do ik = 1, nk_batches
     do ib = 1, nbnd_comp(ik)
        ncol = ncol + 1
        col_slot(ncol) = ik
        col_band(ncol) = ib
     end do
  end do
  if (ncol == 0) return

  !$acc data present(psi, hpsi, npw_comp, nb1_comp) copyin(col_slot(1:ncol), col_band(1:ncol))
  do c0 = 0, ncol-1, fft_cols_eff
     nc = min(fft_cols_eff, ncol-c0)

     !$acc kernels async(clock_thread)
     fft_array_batched(:,1:nc) = cmplx(0._dp, 0._dp, kind=dp)
     !$acc end kernels

     ! G-sphere -> FFT grid
     !$acc parallel loop collapse(2) async(clock_thread) private(jc, ik, ibr)
     do ic = 1, nc
        do ig = 1, npwx
           jc = c0 + ic
           ik = col_slot(jc)
           if (ig .le. npw_comp(ik)) then
              ibr = col_band(jc) + nb1_comp(ik) - 1
              fft_array_batched(dfft%nl(igk_batched(ig,slot2batch(ik))), ic) = psi(ig,ibr,ik)
           end if
        end do
     end do

     pos = 1; rem = nc
     do while (rem > 0)
        k = 1
        do while (2*k <= rem)
           k = 2*k
        end do
        !$acc host_data use_device(fft_array_batched)
        call invfft( 'Wave', fft_array_batched(:,pos:pos+k-1), dfft, k )
        !$acc end host_data
        pos = pos + k; rem = rem - k
     end do

     !$acc parallel loop collapse(2) async(clock_thread)
     do ic = 1, nc
        do ir = 1, nnr
           fft_array_batched(ir,ic) = fft_array_batched(ir,ic) * vloc(ir)
        end do
     end do

     pos = 1; rem = nc
     do while (rem > 0)
        k = 1
        do while (2*k <= rem)
           k = 2*k
        end do
        !$acc host_data use_device(fft_array_batched)
        call fwfft( 'Wave', fft_array_batched(:,pos:pos+k-1), dfft, k )
        !$acc end host_data
        pos = pos + k; rem = rem - k
     end do

     ! FFT grid -> G-sphere, plus the kinetic energy
     !$acc parallel loop collapse(2) async(clock_thread) private(jc, ik, ibr)
     do ic = 1, nc
        do ig = 1, npwx
           jc = c0 + ic
           ik = col_slot(jc)
           if (ig .le. npw_comp(ik)) then
              ibr = col_band(jc) + nb1_comp(ik) - 1
              hpsi(ig,ibr,ik) = ekin_batched(ig,slot2batch(ik)) * psi(ig,ibr,ik) +  &
                                fft_array_batched(dfft%nl(igk_batched(ig,slot2batch(ik))), ic)
           end if
        end do
     end do
  end do
  !$acc wait(clock_thread)
  !$acc end data

end subroutine my_h_psi_batched_iter
