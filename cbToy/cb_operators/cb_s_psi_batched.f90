
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE cb_s_psi_batched( lda, n, m, psi, spsi, i_batch )
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the product of the S
  ! ... matrix with m wavefunctions contained in psi
  !
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  ! ...    spsi  S*psi
  !
  ! --- Wrapper routine: performs bgrp parallelization on non-distributed bands
  ! --- if suitable and required, calls old S\psi routine as cb_s_psi_
  !
  USE cb_module,        ONLY : DP
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: lda, n, m, i_batch
  COMPLEX(DP), INTENT(IN)  :: psi(lda,m)
  COMPLEX(DP), INTENT(OUT) :: spsi(lda,m)
  !
  INTEGER     :: m_start, m_end
  !
  ! band parallelization with non-distributed bands is performed if
  ! 1. enabled (variable use_bgrp_in_hpsi must be set to .T.)
  ! 2. there is more than one band, otherwise there is nothing to parallelize
  !
  call start_clock('s_psi')

  ! don't use band parallelization here
  CALL cb_s_psi_batched_( lda, n, m, psi, spsi , i_batch)
  CALL stop_clock('s_psi')

  RETURN
  !
END SUBROUTINE cb_s_psi_batched
!----------------------------------------------------------------------------
 subroutine cb_s_psi_batched_(npwx,npw,nvec,psi,spsi,i_batch)
!----------------------------------------------------------------------------
! for each input wfc spsi is just psi (use_overlap = .false.), or (1+exp(-ekin))**2 * psi
! The arrays are on the device (OpenACC data of the caller), so the work is done by an OpenACC kernel
! on the queue of the calling thread.

! global variables
  USE cb_module, only : DP
  USE cb_module, only : use_overlap, ekin_batched
  USE mytime, only : clock_thread
  implicit none
! input variables
  integer, intent(IN) :: npwx, npw, nvec, i_batch
  complex(DP),intent(IN) :: psi(npwx,nvec)
  complex(DP),intent(OUT) :: spsi(npwx,nvec)
! local variables
  integer :: ivec, ig

  !$acc parallel loop collapse(2) async(clock_thread) present(psi, spsi)
  do ivec=1,nvec
     do ig=1,npwx
        if (use_overlap .and. ig .le. npw) then
           spsi(ig,ivec) = (1.d0 + exp(-ekin_batched(ig,i_batch)))**2 * psi(ig,ivec)
        else
           spsi(ig,ivec) = psi(ig,ivec)
        end if
     end do
  end do

 end subroutine cb_s_psi_batched_
