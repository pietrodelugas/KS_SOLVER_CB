subroutine my_h_psi_batched(npwx, npw, nbnd, psi, hpsi, i_batch, nk_batches)
  use iso_fortran_env, only: dp=> real64
  use cb_module, only: ekin_batched, vloc, igk_batched, dfft, fft_array_batched, aux_batched
  use fft_interfaces, only: fwfft, invfft 
  use mytime, only: clock_thread, fftw_locker
  use omp_lib, only: omp_set_lock, omp_unset_lock
  implicit none 
  integer,intent(in) :: npwx, npw, nbnd, i_batch
  complex(dp) :: psi(npw,nbnd, nk_batches), hpsi(npw,nbnd, nk_batches) !psi(npwx, nbnd), hpsi(npwx, nbnd) ! We make the arrays as assumed sized arguments
  !complex(dp), optional :: psi_3d(npw,nbnd,nk_batches), hpsi_3d(npw,nbnd,nk_batches) !psi(npwx, nbnd), hpsi(npwx, nbnd) ! We make the arrays as assumed sized arguments
  integer, optional, intent(in) :: nk_batches !! We include an argument for the the size of the batch
  !
  integer :: ig, ibnd, ir, ik ! ik is the index for the element in the batch
  
  print *, "First my_h_psi_batched point : 1"

  !$acc data present(psi, hpsi, fft_array_batched) async(clock_thread)
  do ibnd =1, nbnd
    !if (.not. present(nk_batches)) then
    !   !$acc kernels async(i_batch)
    !   fft_array_batched(:,i_batch)=cmplx(0._dp, 0._dp, kind=dp)
    !   !$acc end kernels 
    !else
    
    print *, "First my_h_psi_batched point: 2"

    !$acc kernels async(clock_thread)
    fft_array_batched(:,1:nk_batches)=cmplx(0._dp, 0._dp, kind=dp)
    !$acc end kernels
    
    print *, "my_h_psi_batched point: 3"

    !$acc parallel loop async(clock_thread)
    do ig = 1, npw
      ! if (igk_batched(ig,i_batch) <= 0 .or. igk_batched(ig,i_batch) > size(dfft%nl)) stop "igk_batched index out of range"
      !if (.not. present(nk_batches)) then
       !  fft_array_batched(dfft%nl(igk_batched(ig,i_batch)),i_batch) = psi(ig,ibnd) 
      !else
      do ik = 1, nk_batches
         fft_array_batched(dfft%nl(igk_batched(ig,ik)),ik) = psi(ig,ibnd,ik)
      end do
    end do
    !if (.not. present(nk_batches)) call omp_set_lock(fftw_locker)
    !!$acc host_data use_device(fft_array_batched) 
    !if (.not. present(nk_batches)) then 
     !   !$acc host_data use_device(fft_array_batched)
     !  call invfft( 'Wave', fft_array_batched(:,i_batch), dfft)
     !   !$acc end host_data
    !else

    print *, "my_h_psi_batched point: 5"
        
        !$acc host_data use_device(fft_array_batched)
        call invfft( 'Wave', fft_array_batched(:,:), dfft, nk_batches)
        !$acc end host_data
    
    print *, "my_h_psi_batched point: 6"

    !end if
    !call omp_unset_lock(fftw_locker) 
    !$acc parallel loop async(clock_thread)
    do ir = 1, dfft%nnr 
       !if (.not. present(nk_batches)) then
       !   fft_array_batched(ir, i_batch) = fft_array_batched(ir,i_batch) * vloc(ir) 
       !else
       do ik = 1, nk_batches
          fft_array_batched(ir, ik) = fft_array_batched(ir,ik) * vloc(ir) 
       end do 
       !end if
    end do 
    
    print *, "my_h_psi_batched point: 7"

    !call omp_set_lock(fftw_locker) 
    !!$acc host_data use_device(fft_array_batched) 
    !if (.not. present(nk_batches)) then
    !    call fwfft('Wave', fft_array_batched(:,i_batch), dfft) 
    !else
    !$acc host_data use_device(fft_array_batched) 
    call fwfft('Wave', fft_array_batched(:,:), dfft, nk_batches) 
    !$acc end host_data

    print *, "my_h_psi_batched point: 8"

    !if (.not. present(nk_batches)) call omp_unset_lock(fftw_locker)
    !$acc parallel loop async(clock_thread)
    do ig = 1, npw
       !if (.not. present(nk_batches)) then
       !   hpsi(ig,ibnd) = ekin_batched(ig,i_batch) * psi(ig,ibnd)   + &
       !                   fft_array_batched(dfft%nl(igk_batched(ig,i_batch)),i_batch)  
       !else
          do ik = 1, nk_batches
             hpsi(ig,ibnd,ik) = ekin_batched(ig,ik) * psi(ig,ibnd,ik)   + &
                             fft_array_batched(dfft%nl(igk_batched(ig,ik)),ik)
          end do
    end do 
  end do
  !$acc end data

  print *, "my_h_psi_batched point: 9"

end subroutine my_h_psi_batched
        
  
   
