subroutine qe_h_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_comp, nb1_comp, psi, hpsi, i_batch, nk_batches, slot2batch)
  !! H*psi for the new correction vectors of all the k-points (batches) still active in the Davidson iterations.
  !!
  !! Slot ik (1..nk_batches) is the position of a batch among the active ones; slot2batch(ik) is the batch it belongs to
  !! (it selects igk_batched / g2kin_batched / spin_batched). For every slot, nbnd_comp(ik) bands starting at column nb1_comp(ik) of
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
  use qe_batch_module, only: g2kin_batched, igk_batched, spin_batched, dfft, fft_array_batched, fft_cols_eff, &
                             vkb_batched, becp_batched, ps_batched, kb2atom, kb2ih
  use qe_ham_module, only: nl, vrs, nkb, deeq, nh, ityp, ofsbeta
  use laxlib_cusolver_handles, only: cublas_handle
  use qe_batch_module, only: qe_cublas_prepare
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
  real(dp) :: dfac
  complex(dp) :: aux, acc
  integer :: m, npw, nb1, b, spin, ikb, na, nt, ih, jh, ofs, info

  nnr = dfft%nnr

  ! npw_comp, nb1_comp and slot2batch are set on the device by the caller (kernels): the non local GEMMs need them on the host
  !$acc update host(npw_comp(1:nk_batches), nb1_comp(1:nk_batches), slot2batch(1:nk_batches))
  ! NOTE: 'Rho' = full 3-D FFT of the grid (no stick pruning yet)

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
     !$acc parallel loop collapse(2) async(clock_thread) private(jc, ik, ibr, aux, dfac)
     do ic = 1, nc
        do ig = 1, npwx
           jc = c0 + ic
           ik = col_slot(jc)
           if (ig .le. npw_comp(ik)) then
              ibr = col_band(jc) + nb1_comp(ik) - 1
              aux = psi(ig,ibr,ik)
              fft_array_batched(nl(igk_batched(ig,slot2batch(ik))), ic) = aux
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
        call invfft( 'Rho', fft_array_batched(:,pos:pos+k-1), dfft, k )
        !$acc end host_data
        pos = pos + k; rem = rem - k
     end do

     !$acc parallel loop collapse(2) async(clock_thread) private(jc)
     do ic = 1, nc
        do ir = 1, nnr
           jc = c0 + ic
           fft_array_batched(ir,ic) = fft_array_batched(ir,ic) * vrs(ir,spin_batched(slot2batch(col_slot(jc))))
        end do
     end do

     pos = 1; rem = nc
     do while (rem > 0)
        k = 1
        do while (2*k <= rem)
           k = 2*k
        end do
        !$acc host_data use_device(fft_array_batched)
        call fwfft( 'Rho', fft_array_batched(:,pos:pos+k-1), dfft, k )
        !$acc end host_data
        pos = pos + k; rem = rem - k
     end do

     ! FFT grid -> G-sphere, plus the kinetic energy
     !$acc parallel loop collapse(2) async(clock_thread) private(jc, ik, ibr, aux, dfac)
     do ic = 1, nc
        do ig = 1, npwx
           jc = c0 + ic
           ik = col_slot(jc)
           if (ig .le. npw_comp(ik)) then
              ibr = col_band(jc) + nb1_comp(ik) - 1
              hpsi(ig,ibr,ik) = g2kin_batched(ig,slot2batch(ik)) * psi(ig,ibr,ik) + &
                                fft_array_batched(nl(igk_batched(ig,slot2batch(ik))), ic)
           end if
        end do
     end do
  end do
  !$acc wait(clock_thread)

  ! Non local part  hpsi += vkb D vkb^H psi, one GEMM pair per active slot on the stream of this thread:
  ! becp = vkb^H psi_new ; ps = D becp (kernel) ; hpsi_new += vkb ps
  if (nkb > 0) then
     call qe_cublas_prepare(i_batch)
     do ik = 1, nk_batches
        m = nbnd_comp(ik); if (m <= 0) cycle
        npw = npw_comp(ik); nb1 = nb1_comp(ik); b = slot2batch(ik)
        !$acc host_data use_device(vkb_batched, psi, becp_batched)
        call myzgemm_h( cublas_handle(i_batch), 'C', 'N', nkb, m, npw, (1._dp,0._dp), vkb_batched(1,1,b), npwx, &
                        psi(1,nb1,ik), npwx, (0._dp,0._dp), becp_batched(1,1,ik), nkb )
        !$acc end host_data
        !$acc parallel loop collapse(2) async(clock_thread) present(ps_batched, becp_batched, deeq, nh, ityp, ofsbeta, kb2atom, kb2ih, spin_batched) private(na, nt, ofs, spin, acc)
        do ic = 1, m
           do ikb = 1, nkb
              na = kb2atom(ikb); nt = ityp(na); ofs = ofsbeta(na); spin = spin_batched(b)
              acc = (0._dp,0._dp)
              do jh = 1, nh(nt)
                 acc = acc + deeq(kb2ih(ikb),jh,na,spin) * becp_batched(ofs+jh,ic,ik)
              end do
              ps_batched(ikb,ic,ik) = acc
           end do
        end do
        !$acc host_data use_device(vkb_batched, ps_batched, hpsi)
        call myzgemm_h( cublas_handle(i_batch), 'N', 'N', npw, m, nkb, (1._dp,0._dp), vkb_batched(1,1,b), npwx, &
                        ps_batched(1,1,ik), nkb, (1._dp,0._dp), hpsi(1,nb1,ik), npwx )
        !$acc end host_data
     end do
     !$acc wait(clock_thread)
  end if
  !$acc end data

end subroutine qe_h_psi_batched_iter
