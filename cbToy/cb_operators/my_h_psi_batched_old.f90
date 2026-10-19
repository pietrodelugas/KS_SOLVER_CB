subroutine my_h_psi_batched_old(npwx, npw, nbnd, psi, hpsi, i_batch)
  use iso_fortran_env, only: dp=> real64
  use cb_module, only: ekin_batched, vloc, igk_batched, dfft, fft_array_batched, aux_batched
  use fft_interfaces, only: fwfft, invfft 
  use mytime, only: clock_thread, fftw_locker
  use omp_lib, only: omp_set_lock, omp_unset_lock
  implicit none 
  integer,intent(in) :: npwx, npw, nbnd, i_batch
  complex(dp) :: psi(npwx, nbnd), hpsi(npwx, nbnd) 
  !
  integer :: ig, ibnd, ir
  !$acc data present(psi, hpsi) async(clock_thread)
  do ibnd =1, nbnd
    !$acc kernels async(i_batch)
    fft_array_batched(:,i_batch)=cmplx(0._dp, 0._dp, kind=dp)
    !$acc end kernels 
    !$acc parallel loop async(clock_thread)
    do ig = 1, npw
      ! if (igk_batched(ig,i_batch) <= 0 .or. igk_batched(ig,i_batch) > size(dfft%nl)) stop "igk_batched index out of range"
      fft_array_batched(dfft%nl(igk_batched(ig,i_batch)),i_batch) = psi(ig,ibnd) 
    end do 
    call omp_set_lock(fftw_locker)
    !$acc host_data use_device(fft_array_batched) 
    call invfft( 'Wave', fft_array_batched(:,i_batch), dfft) 
    !$acc end host_data 
    !call omp_unset_lock(fftw_locker) 
    !$acc parallel loop async(clock_thread)
    do ir = 1, dfft%nnr 
       fft_array_batched(ir, i_batch) = fft_array_batched(ir,i_batch) * vloc(ir) 
    end do 
    !call omp_set_lock(fftw_locker) 
    !$acc host_data use_device(fft_array_batched) 
    call fwfft('Wave', fft_array_batched(:,i_batch), dfft) 
    !$acc end host_data
    call omp_unset_lock(fftw_locker) 
    !$acc parallel loop async(clock_thread)
    do ig = 1, npw
       hpsi(ig,ibnd) = ekin_batched(ig,i_batch) * psi(ig,ibnd)   + &
                       fft_array_batched(dfft%nl(igk_batched(ig,i_batch)),i_batch)  
    end do 
  end do
  !$acc end data
end subroutine my_h_psi_batched_old
