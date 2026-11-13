SUBROUTINE qe_s_psi_batched( npwx, npw, nvec, psi, spsi, i_batch )
  !! S|psi> = psi + vkb qq vkb^H psi for the calling thread (batch i_batch); arrays are on the device.
  !! Only meaningful for ultrasoft/PAW (cegterg calls it when uspp is true); GEMMs on the stream of the thread.
  USE qe_batch_module, ONLY : vkb_batched, becp_s, ps_s, kb2atom, kb2ih, qe_cublas_prepare
  USE qe_ham_module, ONLY : nkb, qq_at, nh, ityp, ofsbeta
  USE laxlib_cusolver_handles, ONLY : cublas_handle
  USE mytime, ONLY : clock_thread
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npwx, npw, nvec, i_batch
  COMPLEX(8), INTENT(IN) :: psi(npwx,nvec)
  COMPLEX(8), INTENT(OUT) :: spsi(npwx,nvec)
  INTEGER :: ic, ikb, jh, na, nt, ofs
  COMPLEX(8) :: acc
  CALL qe_cublas_prepare(i_batch)
  !$acc parallel loop collapse(2) async(clock_thread) present(psi, spsi)
  DO ic = 1, nvec
     DO ikb = 1, npwx
        spsi(ikb,ic) = psi(ikb,ic)
     END DO
  END DO
  IF (nkb > 0) THEN
     !$acc host_data use_device(vkb_batched, psi, becp_s)
     CALL myzgemm_h( cublas_handle(i_batch), 'C', 'N', nkb, nvec, npw, (1.d0,0.d0), vkb_batched(1,1,i_batch), npwx, &
                     psi, npwx, (0.d0,0.d0), becp_s(1,1,i_batch), nkb )
     !$acc end host_data
     !$acc parallel loop collapse(2) async(clock_thread) present(ps_s, becp_s, qq_at, nh, ityp, ofsbeta, kb2atom, kb2ih) private(na, nt, ofs, acc)
     DO ic = 1, nvec
        DO ikb = 1, nkb
           na = kb2atom(ikb); nt = ityp(na); ofs = ofsbeta(na)
           acc = (0.d0,0.d0)
           DO jh = 1, nh(nt)
              acc = acc + qq_at(kb2ih(ikb),jh,na) * becp_s(ofs+jh,ic,i_batch)
           END DO
           ps_s(ikb,ic,i_batch) = acc
        END DO
     END DO
     !$acc host_data use_device(vkb_batched, ps_s, spsi)
     CALL myzgemm_h( cublas_handle(i_batch), 'N', 'N', npw, nvec, nkb, (1.d0,0.d0), vkb_batched(1,1,i_batch), npwx, &
                     ps_s(1,1,i_batch), nkb, (1.d0,0.d0), spsi, npwx )
     !$acc end host_data
  END IF
END SUBROUTINE qe_s_psi_batched
!
SUBROUTINE qe_g_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_comp, nb1_comp, psi, eig, nvec, i_batch, nk_batches, slot2batch)
  !! QE's preconditioner for the correction vectors of all the active slots: x = h_diag - e s_diag,
  !! denm = (1 + x + sqrt(1 + (x-1)^2)) / 2, psi = psi / denm (see cb_g_psi_batched_iter for the layout).
  USE mytime, ONLY : clock_thread
  USE qe_batch_module, ONLY : hdiag_batched, sdiag_batched
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nk_batches, npwx, nbndx, nvec, i_batch
  INTEGER, INTENT(IN) :: npw_comp(nk_batches), nbnd_comp(nk_batches), nb1_comp(nk_batches), slot2batch(nk_batches)
  COMPLEX(8), INTENT(INOUT) :: psi(npwx, nbndx, nk_batches)
  REAL(8), INTENT(IN) :: eig(nvec, nk_batches)
  INTEGER :: ik, ivec, ig, nb_max
  REAL(8) :: x, denm
  nb_max = MAXVAL(nbnd_comp(1:nk_batches))
  !$acc parallel loop collapse(3) async(clock_thread) copyin(nbnd_comp) present(psi, eig, npw_comp, nb1_comp, slot2batch, hdiag_batched, sdiag_batched) private(x, denm)
  DO ik = 1, nk_batches
     DO ivec = 1, nb_max
        DO ig = 1, npwx
           IF (ivec <= nbnd_comp(ik) .AND. ig <= npw_comp(ik)) THEN
              x = hdiag_batched(ig,slot2batch(ik)) - eig(ivec,ik) * sdiag_batched(ig,slot2batch(ik))
              denm = 0.5d0*(1.d0 + x + SQRT(1.d0 + (x-1.d0)*(x-1.d0)))
              psi(ig, nb1_comp(ik)+ivec-1, ik) = psi(ig, nb1_comp(ik)+ivec-1, ik) / denm
           END IF
        END DO
     END DO
  END DO
END SUBROUTINE qe_g_psi_batched_iter
!
SUBROUTINE qe_h_psi_batched(npwx, npw_comp, nbndx, nbnd, psi, hpsi, i_batch, nk_batches)
  !! First H|psi> of the Davidson (called once per group): the iterative routine with every slot active,
  !! all the bands to be done and the first column equal to 1.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nk_batches, npwx, nbnd, i_batch, nbndx
  INTEGER, INTENT(IN) :: npw_comp(nk_batches)
  COMPLEX(8) :: psi(npwx, nbndx, nk_batches), hpsi(npwx, nbndx, nk_batches)
  INTEGER :: nbnd_loc(nk_batches), nb1_loc(nk_batches), s2b_loc(nk_batches), ik
  DO ik = 1, nk_batches
     nbnd_loc(ik) = nbnd; nb1_loc(ik) = 1; s2b_loc(ik) = ik
  END DO
  !$acc data copyin(nb1_loc, s2b_loc)
  CALL qe_h_psi_batched_iter(npwx, npw_comp, nbndx, nbnd_loc, nb1_loc, psi, hpsi, i_batch, nk_batches, s2b_loc)
  !$acc end data
END SUBROUTINE qe_h_psi_batched
!
SUBROUTINE qe_g_psi_batched(npwx, npw, nvec, npol, psi, eig, i_batch)
  !! not used: the batched solver calls qe_g_psi_batched_iter
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npwx, npw, nvec, npol, i_batch
  COMPLEX(8), INTENT(INOUT) :: psi(npwx,nvec)
  REAL(8), INTENT(IN) :: eig(nvec)
  STOP 'qe_g_psi_batched: not used'
END SUBROUTINE qe_g_psi_batched
