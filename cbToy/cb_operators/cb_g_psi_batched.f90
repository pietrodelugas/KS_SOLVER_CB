 subroutine cb_g_psi_batched(npwx, npw, nvec, npol, psi, eig, i_batch)

! global variables
  USE mytime, ONLY: clock_thread
  USE cb_module, ONLY : DP
  USE cb_module, ONLY : ekin_batched
  implicit none
! input variables
  integer, intent(IN) :: npwx, npw, nvec, npol, i_batch
  complex(DP), intent(INOUT) :: psi(npwx,nvec)
  real(DP), intent(IN) :: eig(nvec)
! local variables
  integer :: ivec, ig
  real(DP) :: x, denm
  !$acc  data deviceptr(eig,psi) async(clock_thread) 
  call start_clock('g_psi')
  !$acc parallel loop collapse(2) async(clock_thread) 
  do ivec = 1, nvec
     do ig = 1, npw
        x = (ekin_batched(ig, i_batch) - eig(ivec))
        denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))
        psi (ig, ivec) = psi (ig, ivec) / denm
     enddo
  enddo
  !$acc end data 
  call stop_clock('g_psi')

 end subroutine cb_g_psi_batched




 subroutine cb_g_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_comp, nb1_comp, psi, eig, nvec, i_batch, nk_batches, slot2batch)
!! g_psi for the correction vectors of all the active slots at once (called from the iterative part of cegterg):
!! slot ik is the position of a batch among the active ones, slot2batch(ik) the batch it belongs to; for every slot
!! nbnd_comp(ik) columns starting at nb1_comp(ik) of psi(:,:,ik) are divided by the preconditioner, with the
!! eigenvalues eig(1:nbnd_comp(ik),ik).
  USE mytime, ONLY: clock_thread
  USE cb_module, ONLY : DP
  USE cb_module, ONLY : ekin_batched
  implicit none
  integer, intent(IN) :: nk_batches, npwx, nbndx, nvec, i_batch
  integer, intent(IN) :: npw_comp(nk_batches), nbnd_comp(nk_batches), nb1_comp(nk_batches), slot2batch(nk_batches)
  complex(DP), intent(INOUT) :: psi(npwx, nbndx, nk_batches)
  real(DP), intent(IN) :: eig(nvec, nk_batches)
  integer :: ik, ivec, ig, nb_max
  real(DP) :: x, denm
  nb_max = maxval(nbnd_comp(1:nk_batches))
  call start_clock('g_psi')
  !$acc parallel loop collapse(3) async(clock_thread) copyin(nbnd_comp) present(psi, eig, npw_comp, nb1_comp, slot2batch) private(x, denm)
  do ik = 1, nk_batches
     do ivec = 1, nb_max
        do ig = 1, npwx
           if (ivec .le. nbnd_comp(ik) .and. ig .le. npw_comp(ik)) then
              x = ekin_batched(ig, slot2batch(ik)) - eig(ivec, ik)
              denm = 0.5_dp*(1.d0+x+sqrt(1.d0+(x-1)*(x-1.d0)))
              psi(ig, nb1_comp(ik)+ivec-1, ik) = psi(ig, nb1_comp(ik)+ivec-1, ik) / denm
           end if
        enddo
     enddo
  enddo
  call stop_clock('g_psi')

 end subroutine cb_g_psi_batched_iter
