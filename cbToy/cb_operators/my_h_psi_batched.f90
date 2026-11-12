subroutine my_h_psi_batched(npwx, npw_comp, nbndx, nbnd, psi, hpsi, i_batch, nk_batches)
  !! H*psi of the first nbnd vectors of the nk_batches k-points of the group (first H*psi of the Davidson, called once
  !! per group). It is the iterative routine with every k-point active, all the bands to be done and the first column
  !! of the new vectors equal to 1, so the (k-point, band) pairs are packed into the columns of the batched FFT array
  !! exactly as there.
  use iso_fortran_env, only: dp=> real64
  implicit none
  integer,intent(in) :: nk_batches
  integer,intent(in) :: npwx, nbnd, i_batch, nbndx
  integer,intent(in) :: npw_comp(nk_batches)
  complex(dp) :: psi(npwx, nbndx, nk_batches), hpsi(npwx, nbndx, nk_batches)
  integer :: nbnd_loc(nk_batches), nb1_loc(nk_batches), s2b_loc(nk_batches)
  integer :: ik

  do ik = 1, nk_batches
     nbnd_loc(ik) = nbnd   ! bands to be done for every slot
     nb1_loc(ik)  = 1      ! first column
     s2b_loc(ik)  = ik     ! slot ik is batch ik
  end do

  !$acc data copyin(nb1_loc, s2b_loc)
  call my_h_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_loc, nb1_loc, psi, hpsi, i_batch, nk_batches, s2b_loc)
  !$acc end data

end subroutine my_h_psi_batched
