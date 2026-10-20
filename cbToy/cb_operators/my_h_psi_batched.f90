subroutine my_h_psi_batched(npwx, npw_comp, nbndx, nbnd, psi, hpsi, i_batch, nk_batches)
  use iso_fortran_env, only: dp=> real64
  use cb_module, only: ekin_batched, vloc, igk_batched, dfft, fft_array_batched, aux_batched
  use fft_interfaces, only: fwfft, invfft 
  use mytime, only: clock_thread, fftw_locker
  use omp_lib, only: omp_set_lock, omp_unset_lock
  implicit none 
  integer,intent(in) :: npwx, nbnd, i_batch, nbndx !,npw
  integer,intent(in) :: npw_comp(nk_batches)
  complex(dp) :: psi(npwx, nbndx, nk_batches), hpsi(npwx, nbndx, nk_batches) !psi(npwx, nbnd), hpsi(npwx, nbnd) ! We make the arrays as assumed sized arguments
  !complex(dp), optional :: psi2D(npw,nbnd), hpsi2D(npw,nbnd) !psi(npwx, nbnd), hpsi(npwx, nbnd) ! We make the arrays as assumed sized arguments
  integer, intent(in) :: nk_batches !! We include an argument for the the size of the batch
  !
  integer :: ig, ibnd, ir, ik ! ik is the index for the element in the batch
  
  print *, "First my_h_psi_batched point : 1"

  !$acc data present(psi, hpsi, npw_comp) async(clock_thread)
  do ibnd =1, nbnd
        !print *, "First my_h_psi_batched point: 2"
        !print *, "Ibatch ! : ", i_batch

        !$acc kernels async(clock_thread)
        fft_array_batched(:,1:nk_batches)=cmplx(0._dp, 0._dp, kind=dp)
        !$acc end kernels
    
        !print *, "my_h_psi_batched point: 3"
    
        !print *, "Ibatch ! : ", i_batch
        !print *, "NK Batches :", nk_batches
        !$acc parallel loop async(clock_thread)
        do ik = 1, nk_batches
            do ig = 1, npw_comp(ik)
               !if (igk_batched(ig,i_batch) <= 0 .or. igk_batched(ig,i_batch) > size(dfft%nl)) stop "igk_batched index out of range"
               !do ik = 1, nk_batches
               fft_array_batched(dfft%nl(igk_batched(ig,ik)),ik) = psi(ig,ibnd,ik)
               !end do
            end do
        end do
        !if (.not. present(nk_batches)) call omp_set_lock(fftw_locker)
        !!$acc host_data use_device(fft_array_batched) 
        !print *, "my_h_psi_batched point: 5"
        !print *, "Ibatch ! : ", i_batch

        !$acc host_data use_device(fft_array_batched)
        call invfft( 'Wave', fft_array_batched(:,:), dfft, nk_batches)
        !$acc end host_data
    
        !print *, "my_h_psi_batched point: 6"
        !print *, "Ibatch ! : ", i_batch
    
        !call omp_unset_lock(fftw_locker) 
        !$acc parallel loop async(clock_thread)
        do ir = 1, dfft%nnr 
          !do ik = 1, nk_batches
          fft_array_batched(ir, 1:nk_batches) = fft_array_batched(ir,1:nk_batches) * vloc(ir) 
          !end do 
        end do 
    
        !print *, "my_h_psi_batched point: 7"

        !print *, "Ibatch ! : ", i_batch

        !call omp_set_lock(fftw_locker) 
        !!$acc host_data use_device(fft_array_batched) 
        !$acc host_data use_device(fft_array_batched) 
        call fwfft('Wave', fft_array_batched(:,:), dfft, nk_batches) 
        !$acc end host_data
    
        !print *, "my_h_psi_batched point: 8"
        !print *, "Ibatch ! : ", i_batch
    
        !$acc wait(clock_thread)

        !if (.not. present(nk_batches)) call omp_unset_lock(fftw_locker)
        !$acc parallel loop async(clock_thread)
        do ik = 1, nk_batches
           do ig = 1, npw_comp(ik)
           !do ik = 1, nk_batches
              hpsi(ig,ibnd,ik) = ekin_batched(ig,ik) * psi(ig,ibnd,ik)   + &
                                 fft_array_batched(dfft%nl(igk_batched(ig,ik)),ik)
           end do
        end do 
    
  end do
  !$acc end data
  
  !$acc wait(clock_thread)

  !print *, "my_h_psi_batched point: 9"

  !print *, "Ibatch ! : ", i_batch

end subroutine my_h_psi_batched
