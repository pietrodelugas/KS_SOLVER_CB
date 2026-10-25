subroutine my_h_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_comp, nb1_comp, psi, hpsi, i_batch, nk_batches, slot2batch)
  use iso_fortran_env, only: dp=> real64
  use cb_module, only: ekin_batched, vloc, igk_batched, dfft, fft_array_batched, aux_batched
  use fft_interfaces, only: fwfft, invfft 
  use mytime, only: clock_thread, fftw_locker
  use omp_lib, only: omp_set_lock, omp_unset_lock
  implicit none 
  integer,intent(in) :: nk_batches
  integer,intent(in) :: npwx, i_batch, nbndx !,npw
  integer,intent(in) :: npw_comp(nk_batches)
  integer,intent(in) :: nbnd_comp(nk_batches) ! We need this vector to take into account the various values related to each omp thread
  integer,intent(in) :: nb1_comp(nk_batches)
  integer,intent(in) :: slot2batch(nk_batches)
  complex(dp) :: psi(npwx, nbndx, nk_batches), hpsi(npwx, nbndx, nk_batches) ! We make the arrays as assumed sized arguments
  integer :: ibnd_real_comp(nk_batches)
  integer :: ig, ibnd, ir, ik, nbnd_max, ibnd_real, n_active, ind_active_thr, init_ind ! ik is the index for the element in the batch
  integer :: n_active_comp(nk_batches)
  !print *, "First my_h_psi_batched point : 1"
  
  nbnd_max = maxval(nbnd_comp)
  n_active = 0 ! We initialize the counter for the effective number of elements in the batch to be computed 
  ibnd_real = 0
  ind_active_thr = 0
  init_ind = 0

  !print *, "NPW_COMP ELEMENTS: ", npw_comp
  !print *, "NB1_COMP ELEMENTS: ", nb1_comp

  !$acc data present(psi, hpsi, npw_comp, nb1_comp) async(clock_thread)
  do ibnd = 1, nbnd_max
        !$acc kernels async(clock_thread)
        fft_array_batched(:,1:nk_batches)=cmplx(0._dp, 0._dp, kind=dp)
        !$acc end kernels
        
        n_active = 0 ! We initialize the counter for the threads still active
        
        ! We update on the host the nb1_comp array in order to make possible to compute the true indices in the next loop 
        ! for both unified and distributed architercture hardwares :
        !$acc update host(nb1_comp)
        
        !!$acc kernels async(clock_thread)
        do ik = 1, nk_batches
           if (ibnd .le. nbnd_comp(ik)) then
              n_active = n_active + 1
              ibnd_real = ibnd + nb1_comp(ik) - 1
              !print *, "NB1 FIRST DO LOOP: ", nb1_comp(ik)
              !print *, "IBND REAL FIRST DO LOOP: ", ibnd_real
              n_active_comp(ik) = n_active
              ibnd_real_comp(ik) = ibnd_real
              ind_active_thr = ik
           end if
        end do
        !!$acc end kernels

        !$acc wait(clock_thread)

        n_active = 0
        
        !$acc parallel loop async(clock_thread) private(init_ind)
        do ik = 1, nk_batches
            if (ibnd .le. nbnd_comp(ik)) then
               !print *, "IK VALUE: ", ik
               !print *, "NB1 VALUE: ", nb1_comp(ik)
               !print *, "IBND_REAL CURRENT ik: ", ibnd_real_comp(ik)
               !print *, "NACTIVE CURRENT ik: ", n_active_comp(ik)
               init_ind = slot2batch(ik) ! We take into account the true k index for the computing of the igk and ekin arrays 
               do ig = 1, npw_comp(ik)
                  fft_array_batched(dfft%nl(igk_batched(ig,init_ind)), n_active_comp(ik)) = psi(ig,ibnd_real_comp(ik),ik)
               end do
            end if
        end do

        !print *, "THREAD: ", clock_thread
        !print *, "POINT AFTER THE SETTING BEFORE FFT CALL ITERATIVE SUBROUTINE"
        !print *, "NACTIVE BEFORE THE INVFFT CALL: ", ind_active_thr
        !print *, "HOMANY VALUE NEW SUBROUTINE: ", n_active_comp(ind_active_thr)
        !$acc host_data use_device(fft_array_batched)
        call invfft( 'Wave', fft_array_batched(:,:), dfft, n_active_comp(ind_active_thr))
        !$acc end host_data
    
        !$acc parallel loop async(clock_thread)
        do ir = 1, dfft%nnr 
          fft_array_batched(ir, 1:n_active_comp(ind_active_thr)) = fft_array_batched(ir,1:n_active_comp(ind_active_thr)) * vloc(ir) 
        end do
    
        !$acc host_data use_device(fft_array_batched) 
        call fwfft('Wave', fft_array_batched(:,:), dfft, n_active_comp(ind_active_thr)) 
        !$acc end host_data
    
        !!We apply the same logics that were applied for the first do loop:
        n_active = 0 !We reset the counter of the active threads 

        !$acc parallel loop async(clock_thread) private(init_ind)
        do ik = 1, nk_batches
           if (ibnd .le. nbnd_comp(ik)) then
              do ig = 1, npw_comp(ik)
                 init_ind = slot2batch(ik)
                 hpsi(ig,ibnd_real_comp(ik),ik) = ekin_batched(ig,init_ind) * psi(ig,ibnd_real_comp(ik),ik) +  &
                                                  fft_array_batched(dfft%nl(igk_batched(ig,init_ind)), n_active_comp(ik))
              end do
           end if
        end do 
  
  end do
  !$acc end data

end subroutine my_h_psi_batched_iter
