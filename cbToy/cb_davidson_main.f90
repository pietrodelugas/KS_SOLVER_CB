program cb_davidson_main

! global variables
   USE cb_module
#if defined(__MPI)
   use mp_global,            ONLY : mp_startup, mp_global_end
   use mp_world,             ONLY : world_comm
   use mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
#endif
   use mytime,               only: t0cpu, clock_label, nclock, clock_thread              
#if defined(__OPENMP)
   use omp_lib,               only: omp_get_thread_num
#endif
   !!use nvpl_lapack,         only: nvpl_lapack_set_num_threads
   implicit none
   !
   !include 'laxlib.fh'
   !
#if defined(__PAPI) 
   include 'f90papi.h'
   integer EventSet
   integer :: events(3)
   integer(8) :: values(3), values_cegterg_start(3), values_cegterg_stop(3), values_cegterg_tot(3) 
   integer, external :: PAPIF_library_init
   integer :: retval, inval
#endif
! local variables (used in the call to cegterg )
   logical, parameter :: gamma_only = .false. ! general k-point version
   complex(DP), allocatable :: evc(:,:), evc_batched(:,:,:) 
   real(dp), allocatable :: eig(:), eig_batched(:,:) 
   integer, parameter :: npol=1
   integer :: notcnv, dav_iter, nhpsi
   integer, allocatable :: notcnv_batched(:), dav_iter_batched(:), nhpsi_batched(:)
   logical :: overlap = .false. , lrot =.false.
! additional local variables
   real(dp) :: ref=0.d0
   integer :: i_batch, ik
#if defined(__MPI)
! local paralelization variables
   integer :: ndiag     ! input value of processors in the diagonalization group
   logical :: do_distr_diag_in_band_group = .false. ! whether or not the parallel diagonalization is performed inside the
                                                    ! band group or at the parallelization level above.
#endif
!------------------------------------------------------------------------
   external my_h_psi_batched, cb_h_psi, cb_s_psi_batched, cb_g_psi_batched
   external cb_s_psi, cb_g_psi
!  subroutine cb_h_psi(npwx,npw,nvec,psi,hpsi)  computes H*psi
!  subroutine cb_s_psi(npwx,npw,nvec,psi,spsi)  computes S*psi (if needed)
!  subroutine cb_g_psi(npwx,npw,nvec,psi,eig)   computes G*psi -> psi

#if defined(__PAPI)
  EventSet = PAPI_NULL
  events(1) = PAPI_TOT_CYC
  events(2) = PAPI_TOT_INS
  events(3) = PAPI_VEC_INS
  values = 0_8
  values_cegterg_tot = 0_8
  inval = PAPI_VER_CURRENT
  retval = PAPIF_library_init(inval) 
  if (retval .ne. PAPI_VER_CURRENT)  print *, 'PAPI init error, retval=', retval
  call PAPIF_create_eventset(EventSet, retval)
  call PAPIF_add_events(EventSet, events, 3, retval) 
#endif


#if defined(__MPI)
! this call creates the parallel communicators in the MAIN code 
  call mp_startup ( ndiag, diag_in_band_group = do_distr_diag_in_band_group )   
!--- THIS PART IS RELEVANT FOR THE PARALLEL USE OF THE ROUTINE IN KS_Solvers/Davidson -------------------------!
! this set the mpi communicators used internally by the routines in the Davidson library
! it passes 1) the top parent level communicator (could be different from world_comm)
!           2) the sub-communicator of the band group 
!           3) the communicator used across band groups
!           4) whether the distributed diagonalization is performed inside the band group or at the top level
 call set_mpi_comm_4_solvers( world_comm, intra_bgrp_comm, inter_bgrp_comm )

!--------------------------------------------------------------------------------------------------------------!
#endif

   nk_batches = 1 
   !$omp parallel num_threads(nk_batches) default(shared)  shared(t0cpu, nclock, clock_label) 
   call init_clocks(.true.)
   !$omp end parallel


 
   call start_clock('global') 
#if defined(__PAPI)
   call PAPIF_start(EventSet, retval)
#endif
   allocate(npw_batched(nk_batches)) 
   allocate(notcnv_batched(nk_batches), dav_iter_batched(nk_batches), nhpsi_batched(nk_batches))
   call input(gamma_only)
   call ggen(gamma_only)
   call set_cb_potential

   if (use_overlap) write(*,*) '** TEST:  CB hamiltonian modified so as to need an overlap matrix **'
   overlap = use_overlap

   allocate( evc_batched(npwx,nbnd,nk_batches), eig_batched(nbnd,nk_batches) )
   allocate( fft_array_batched(dfft%nnr, nk_batches), aux_batched(dfft%nnr, nk_batches) )
   allocate (evc(npwx, nbnd), eig(nbnd)) 
   !$acc enter data create(evc_batched, eig_batched, fft_array_batched, aux_batched)

   do ik =1,nks, nk_batches
     call start_clock('davidson')
     !$omp parallel num_threads(nk_batches) default(shared) private(i_batch) shared(t0cpu, nclock, clock_label) 
     !$omp do
     do i_batch = 1, min(nk_batches, nks - ik +1) 
       !clock thread is declared threadprivate in the module 
       clock_thread = i_batch  
#if defined(__OPENMP) 
       print '("First loop, batch ",3I5)', i_batch, clock_thread, omp_get_thread_num()  
#else
       print '("First loop, batch ",2I5)', i_batch, clock_thread
#endif
       current_k = ik + i_batch -1   
       call start_clock('init_data') 
       call init_k(current_k, i_batch) 
       !$acc update device(igk_batched(:,i_batch)) 
       call init_random_wfcs(npw_batched(i_batch), npwx, nbnd, evc_batched(1,1,i_batch),i_batch)  
       !$acc update device(evc_batched(:,:,i_batch)) 
       call stop_clock('init_data')     
#if defined(__PAPI) 
       call PAPIF_read(EventSet, values_cegterg_start, retval) 
#endif
       !$acc host_data use_device(eig_batched) 
       call cegterg( my_h_psi_batched, cb_s_psi_batched, overlap, cb_g_psi_batched, &
                      npw_batched(i_batch), npwx, nbnd, nbndx, npol, evc_batched(1,1,i_batch), ethr, &
                      eig_batched(1,i_batch), btype, notcnv_batched(i_batch), lrot, dav_iter_batched(i_batch), & 
                      nhpsi_batched(i_batch), i_batch )
       !$acc end host_data 
#if defined(__PAPI) 
     CALL PAPIF_read(EventSet, values_cegterg_stop, retval)
     values_cegterg_tot = values_cegterg_tot + values_cegterg_stop - values_cegterg_start
#endif 
     end do 
     !$omp end parallel 
     call stop_clock('davidson') 
     !$omp barrier   

     !$acc update self(eig_batched)
     ! Second loop: Process batches sequentially
     do i_batch =1, min(nk_batches, nks - ik +1 )  
        print '("Second loop, batch ",I5)', i_batch 
        current_k = ik + i_batch -1   
        
!--- THIS IS THE RELEVANT CALL TO THE ROUTINE IN KS_Solvers/Davidson ------------------------------------------!
!--------------------------------------------------------------------------------------------------------------!

        
        
        if (energy_shift .and. current_k==1) ref=eig_batched(4*ncell**3,i_batch)
     
        call write_bands(eig_batched(1,i_batch),ref)
        write (stdout,*) 'batch', i_batch, 'dav_iter, nhpsi, notcnv, ethr ', &
                         dav_iter, nhpsi, notcnv, ethr
     end do 
   end do
   
   !$acc exit data delete(evc, eig, fft_array_batched, aux_batched)
   !$acc exit data delete(dfft, dfft%nl, dfft%nnr, igk, vloc) 
#if defined(__PAPI) 
   call PAPIF_stop(EventSet, values, retval) 
#endif
   call stop_clock('global') 
   deallocate( eig )
   deallocate( evc )
   deallocate( evc_batched, eig_batched )
   deallocate( fft_array_batched, aux_batched )
   deallocate( notcnv_batched, dav_iter_batched, nhpsi_batched )
   call print_clock('davidson')
   call print_clock( 'cegterg' )
   call print_clock( 'cegterg:init' )
   call print_clock( 'cegterg:diag' )
   call print_clock( 'cegterg:update' )
   call print_clock( 'cegterg:overlap' )
   call print_clock( 'cegterg:last' )

   call print_clock('h_psi')
   call print_clock('s_psi')
   call print_clock('g_psi')
! 
  write (6,*)
  write (6,*) ' general FFT  routines'
  call print_clock('fftw')
  call print_clock('ffts')
  call print_clock('global') 
#if defined(__PAPI)
  print '(A, I0)', 'Total cycles:  ', values(1) 
  print '(A, I0)', 'Total instructions:', values(2) 
  print '(A, I0)', 'Vector instructions:', values(3) 
  print '(A, F12.5)', 'IPC:   ', real(values(2))/real(values(1)+1) 
  print '(A, I0)', 'Cegterg cycles:  ', values_cegterg_tot(1) 
  print '(A, I0)', 'Cegterg instructions:', values_cegterg_tot(2) 
  print '(A, I0)', 'Cegterg vec_ins:', values_cegterg_tot(3) 
  print '(A, F12.5)', 'IPC:   ', real(values_cegterg_tot(2))/real(values_cegterg_tot(1)+1)  
#endif

#if defined(__MPI)
   call mp_global_end( )
#endif
   call f_print_rmss() 
   end program cb_davidson_main
