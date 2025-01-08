program h_psi_checker
  use iso_fortran_env, only: dp => real64
  use cb_module, except=>dp 
#if defined(__MPI)
  use mp_global,            ONLY : mp_startup, mp_global_end
  use mp_world,             ONLY : world_comm, mpime
  use mp, only:   mp_sum
#endif
  implicit none 
  integer :: diag
  logical :: diag2
  real(dp) :: eig
  complex(dp), allocatable :: psi(:,:), hpsi(:,:)  
  integer :: iun, ibnd, ib, ik
  character,external :: int_to_char 
#if defined(__MPI)
  call mp_startup(diag, diag2) 
#endif
  
  call input(.false.)
  call ggen(.false.)
  call set_cb_potential() 
  allocate(psi(npwx, nbnd), hpsi(npwx, nbnd))  
  do ik = 1, nks
    call init_k()
    open(newunit=iun, file='wfcdat.'//int_to_char(ik),form='unformatted') 
    do ibnd=1, nbnd 
      read(iun) psi(1:npw,ibnd) 
    end do 
    close(iun) 
    !call init_random_wfcs(npw,npwx,nbnd,psi)
    call my_h_psi(npwx, npw, nbnd, psi, hpsi) 
    do ib = 1, min(nbnd,20)
      eig =  13.5*dot_product(psi(:,ib),hpsi(:,ib)) 
      !call mp_sum(eig, world_comm)  
      if (mpime == 0) print *, ik, ib, eig
    end do
  end do   
end program 
