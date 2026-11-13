PROGRAM qe_davidson_main
  !! Batched Davidson (cegterg of KS_Solvers) on the Hamiltonian of a QE run read from the dump files.
  !! Input (namelist &qe on stdin): dumpdir, nk_batches (k-points per group = OpenMP threads), david, ethr, fft_cols_max,
  !!   nks_max (use only the first nks_max k-points, 0 = all), out (eigenvalue file, eV)
  USE qe_ham_module
  USE qe_batch_module
  USE qe_kpoints_input
  USE cegterg_mod
  USE mytime, ONLY : t0cpu, clock_label, nclock, clock_thread, clock_cuda_stream, fftw_locker, cegterg_locker
  USE omp_lib, ONLY : omp_get_thread_num, omp_init_lock, omp_get_wtime
  USE openacc, ONLY : acc_get_cuda_stream, c_devptr
  USE laxlib_cusolver_handles, ONLY : initialize_cusolver_handles, initialize_laxlib_cuda_stream, &
                                      initialize_cublas_handles, finalize_cublas_handles
  USE fft_scalar_cuFFT
  IMPLICIT NONE
  INTEGER, PARAMETER :: npol = 1
  REAL(dp), PARAMETER :: ry2ev = 13.605693122994d0
  CHARACTER(LEN=256) :: dumpdir, out
  CHARACTER(LEN=8) :: kmode = 'dump'      ! 'dump', 'list', 'path' or 'mp' (see qe_kpoints_input)
  INTEGER :: nk1 = 1, nk2 = 1, nk3 = 1, k1 = 0, k2 = 0, k3 = 0
  INTEGER :: nk_batches, david, fft_cols_max, nks_max
  REAL(dp) :: ethr
  NAMELIST /qe/ dumpdir, nk_batches, david, ethr, fft_cols_max, nks_max, out, start_from_dump, nonlocal_batched, kmode, nk1, nk2, nk3, k1, k2, k3
  COMPLEX(dp), ALLOCATABLE :: evc_batched(:,:,:)
  REAL(dp), ALLOCATABLE :: eig_batched(:,:), eig_all(:,:), eig_new_comp(:,:), ew_old_comp(:,:), ew_c(:,:), ew_comp_itr(:,:)
  INTEGER, ALLOCATABLE :: btype(:), notcnv_batched(:), dav_iter_batched(:), nhpsi_batched(:), conv_comp(:,:)
  COMPLEX(dp), ALLOCATABLE :: hc_c(:,:,:), sc_c(:,:,:), vc_c(:,:,:), hc_comp_itr(:,:,:), sc_comp_itr(:,:,:), vc_comp_itr(:,:,:)
  COMPLEX(dp), ALLOCATABLE :: psi_fft_comp(:,:,:), hpsi_fft_comp(:,:,:)
  LOGICAL, ALLOCATABLE :: done_comp(:), final_check_arr(:)
  INTEGER, ALLOCATABLE :: nbase_comp(:), kdim_comp(:), notcnv_comp(:), npw_comp(:), nb1_comp(:), slot2batch(:), my_n_comp(:)
  INTEGER, ALLOCATABLE :: npw_batched(:)
  TYPE(c_devptr), ALLOCATABLE :: ptr_hc(:), ptr_vc(:), ptr_sc(:), ptr_psi(:), ptr_hpsi(:), ptr_spsi(:), ptr_psi_result(:), &
                                 ptr_hpsi_nb1(:), ptr_evc(:)
  INTEGER :: nbndx, ik0, n_k, i_batch, ib, ig, i_b, ind_p, ik, nks, ios, iun, ngroup
  REAL(dp) :: t_start, t_group, t_dav, rr, ri, t0, t_load, t_start_vec
  LOGICAL :: overlap, lrot
  INTEGER :: seed_n
  LOGICAL :: start_from_dump = .FALSE.   ! use the converged wavefunctions of the dump as start vectors (diagnostic)
  INTEGER, ALLOCATABLE :: igk_t(:)
  REAL(dp), ALLOCATABLE :: g_t(:), h_t(:), s_t(:), et_t(:)
  COMPLEX(dp), ALLOCATABLE :: vk_t(:,:), ev_t(:,:), hp_t(:,:), sp_t(:,:), hn_t(:,:)
  INTEGER, ALLOCATABLE :: seed(:)
  EXTERNAL qe_h_psi_batched, qe_s_psi_batched, qe_g_psi_batched, qe_h_psi_batched_iter, qe_g_psi_batched_iter

  dumpdir = 'dump'; nk_batches = 4; david = 4; ethr = 1.d-8; fft_cols_max = 0; nks_max = 0; out = 'eigenvalues_eV.txt'
  READ(5, nml=qe)
  CALL qe_read_global(TRIM(dumpdir))
  CALL qe_set_kpoints(kmode, [nk1, nk2, nk3], [k1, k2, k3], 5)
  nks = nkstot; IF (nks_max > 0) nks = MIN(nks, nks_max)
  nk_batches = MIN(nk_batches, nks)
  IF (david < 2) STOP 'david must be at least 2'
  nbndx = david * nbnd
  overlap = okvan; lrot = .FALSE.
  ALLOCATE(btype(nbnd)); btype = 1
  WRITE(*,'("qeToy: ",i5," k-points of ",i5," nbnd=",i3," npwx=",i5," nkb=",i3," grid ",3i4," okvan=",l1," nk_batches=",i3," david=",i2," ethr=",es8.1)') &
       nks, nkstot, nbnd, npwx, nkb, nr1, nr2, nr3, okvan, nk_batches, david, ethr

  CALL qe_setup_grid()
  CALL initialize_cusolver_handles(nk_batches)
  CALL initialize_cublas_handles(nk_batches)
  !$omp parallel num_threads(nk_batches) default(shared) shared(t0cpu, nclock, clock_label)
  CALL init_clocks(.TRUE.)
  CALL omp_init_lock(cegterg_locker)
  CALL omp_init_lock(fftw_locker)
  !$omp do
  DO i_batch = 1, nk_batches
     clock_thread = i_batch
     clock_cuda_stream = acc_get_cuda_stream(clock_thread)
     CALL initialize_laxlib_cuda_stream(clock_cuda_stream, clock_thread)
  END DO
  !$omp end parallel

  ! columns of one batched FFT block: power of two, memory budget 1 GiB
  IF (fft_cols_max > 0) THEN
     fft_cols_eff = MIN(fft_cols_max, nbnd*nk_batches)
  ELSE
     fft_cols_eff = nbnd*nk_batches
     DO WHILE (fft_cols_eff > 1 .AND. 16_8*nnr*fft_cols_eff > 2_8**30); fft_cols_eff = fft_cols_eff/2; END DO
  END IF
  i_b = 1; DO WHILE (2*i_b <= fft_cols_eff); i_b = 2*i_b; END DO
  fft_cols_eff = i_b
  dim_fft = nnr
  i_b = 1
  DO WHILE (i_b <= fft_cols_eff)
     CALL init_plan_cfft3d_gpu(nr1, nr2, nr3, nr1x, nr2x, nr3x, i_b, ind_p)
     i_b = 2*i_b
  END DO
  CALL qe_alloc_batch(nk_batches, fft_cols_eff, nbndx)
  IF (.NOT. kcompute) CALL qe_open_k(TRIM(dumpdir))

  ALLOCATE(evc_batched(npwx,nbnd,nk_batches), eig_batched(nbnd,nk_batches), eig_all(nbnd,nks))
  ALLOCATE(npw_batched(nk_batches), notcnv_batched(nk_batches), dav_iter_batched(nk_batches), nhpsi_batched(nk_batches))
  ALLOCATE(hc_c(nbnd,nbnd,nk_batches), sc_c(nbnd,nbnd,nk_batches), vc_c(nbnd,nbnd,nk_batches), ew_c(nbnd,nk_batches))
  ALLOCATE(done_comp(nk_batches), final_check_arr(nk_batches), nbase_comp(nk_batches), my_n_comp(nk_batches))
  ALLOCATE(kdim_comp(nk_batches), notcnv_comp(nk_batches))
  ALLOCATE(psi_fft_comp(npwx*npol,nbndx,nk_batches), hpsi_fft_comp(npwx*npol,nbndx,nk_batches))
  ALLOCATE(npw_comp(nk_batches), nb1_comp(nk_batches), slot2batch(nk_batches))
  ALLOCATE(conv_comp(nbnd,nk_batches), ew_old_comp(nbnd,nk_batches), eig_new_comp(nbnd,nk_batches))
  ALLOCATE(ptr_hc(nk_batches), ptr_vc(nk_batches), ptr_sc(nk_batches), ptr_psi(nk_batches), ptr_hpsi(nk_batches), &
           ptr_spsi(nk_batches), ptr_psi_result(nk_batches), ptr_hpsi_nb1(nk_batches), ptr_evc(nk_batches))
  !$acc enter data create(evc_batched, eig_batched)
  !$acc enter data copyin(btype)
  !$acc enter data create(hc_c, sc_c, vc_c, ew_c, psi_fft_comp, hpsi_fft_comp, npw_comp, nb1_comp, slot2batch)
  !$acc enter data create(eig_new_comp, conv_comp, ew_old_comp)

  CALL random_seed(size=seed_n); ALLOCATE(seed(seed_n)); seed = 20251113; CALL random_seed(put=seed)
  t_start = omp_get_wtime(); t_dav = 0; ngroup = 0; t_load = 0; t_start_vec = 0
  DO ik0 = 0, nks-1, nk_batches
     n_k = MIN(nk_batches, nks - ik0)
     ngroup = ngroup + 1
     done_comp = .FALSE.
     t0 = omp_get_wtime()
     CALL qe_load_group(ik0, n_k)
     t_load = t_load + omp_get_wtime() - t0
     t0 = omp_get_wtime()
     ! start wavefunctions: random, damped with the kinetic energy (zero beyond npw)
     evc_batched = (0.d0,0.d0)
     DO i_batch = 1, n_k
        npw_batched(i_batch) = ngk(ik0+i_batch)
        DO ib = 1, nbnd
           DO ig = 1, npw_batched(i_batch)
              CALL random_number(rr); CALL random_number(ri)
              evc_batched(ig,ib,i_batch) = EXP(-g2kin_batched(ig,i_batch)) * CMPLX(rr-0.5d0, ri-0.5d0, dp)
           END DO
        END DO
     END DO
     IF (start_from_dump) THEN
        IF (.NOT. ALLOCATED(igk_t)) ALLOCATE(igk_t(npwx), g_t(npwx), h_t(npwx), s_t(npwx), et_t(nbnd), vk_t(npwx,MAX(nkb,1)), &
                                             ev_t(npwx,nbnd), hp_t(npwx,nbnd), sp_t(npwx,nbnd), hn_t(npwx,nbnd))
        DO i_batch = 1, n_k
           CALL qe_read_k(ik0+i_batch, igk_t, g_t, h_t, s_t, vk_t, et_t, ev_t, hp_t, sp_t, hn_t)
           evc_batched(:,:,i_batch) = ev_t
        END DO
     END IF
     !$acc update device(evc_batched)
     t_start_vec = t_start_vec + omp_get_wtime() - t0
     t_group = omp_get_wtime()
     !$omp parallel num_threads(n_k) default(shared) private(i_batch) shared(t0cpu, nclock, clock_label)
     !$omp do
     DO i_batch = 1, n_k
        clock_thread = i_batch
        clock_cuda_stream = acc_get_cuda_stream(clock_thread)
        !$acc host_data use_device(eig_batched(1,i_batch))
        CALL cegterg( qe_h_psi_batched, qe_s_psi_batched, overlap, qe_g_psi_batched, qe_h_psi_batched_iter, qe_g_psi_batched_iter, &
                      npw_batched(i_batch), npwx, nbnd, nbndx, npol, evc_batched(1,1,i_batch), ethr, &
                      eig_batched(1,i_batch), btype, notcnv_batched(i_batch), lrot, dav_iter_batched(i_batch), &
                      nhpsi_batched(i_batch), i_batch, n_k, hc_c, sc_c, vc_c, ew_c, done_comp, final_check_arr, &
                      hc_comp_itr, sc_comp_itr, vc_comp_itr, ew_comp_itr, nbase_comp, kdim_comp, notcnv_comp, &
                      ptr_hc, ptr_vc, ptr_sc, ptr_psi, ptr_hpsi, ptr_hpsi_nb1, ptr_spsi, ptr_psi_result, ptr_evc, my_n_comp, &
                      psi_fft_comp, hpsi_fft_comp, npw_comp, nb1_comp, slot2batch, eig_new_comp, conv_comp, ew_old_comp )
        !$acc end host_data
     END DO
     !$omp end parallel
     !$acc wait
     !$acc update self(eig_batched)
     t_dav = t_dav + omp_get_wtime() - t_group
     DO i_batch = 1, n_k
        eig_all(:,ik0+i_batch) = eig_batched(:,i_batch) * ry2ev
     END DO
     WRITE(*,'("group ",i4,": k ",i5,"-",i5,"  iterations ",i3,"-",i3,"  not converged ",i3,"  time ",f8.3," s")') ngroup, ik0+1, ik0+n_k, &
          MINVAL(dav_iter_batched(1:n_k)), MAXVAL(dav_iter_batched(1:n_k)), SUM(notcnv_batched(1:n_k)), omp_get_wtime() - t_group
  END DO
  WRITE(*,'("time split: k-point setup + transfer ",f8.3," s, start vectors ",f8.3," s, Davidson ",f8.3," s")') t_load, t_start_vec, t_dav
  WRITE(*,'("Davidson time ",f9.3," s for ",i5," k-points; total ",f9.3," s")') t_dav, nks, omp_get_wtime() - t_start

  OPEN(NEWUNIT=iun, FILE=TRIM(out), STATUS='replace')
  WRITE(iun,'("# qeToy eigenvalues (eV): k (2pi/alat), weight, bands")')
  DO ik = 1, nks
     WRITE(iun,'(3es24.15e3,es24.15e3,*(f12.6))') xk(:,ik), wk(ik), eig_all(:,ik)
  END DO
  CLOSE(iun)
  CALL qe_close_k()
  CALL finalize_cublas_handles()
END PROGRAM qe_davidson_main
