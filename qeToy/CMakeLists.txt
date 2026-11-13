set(src_qe_lib
	qe_lib/qe_dump_h5.f90
	qe_lib/qe_ham_module.f90
	qe_lib/qe_batch_module.f90
	qe_operators/qe_h_psi_batched_iter.f90
	qe_operators/qe_s_psi_g_psi_batched.f90)
qe_enable_cuda_fortran("${src_qe_lib}")
qe_add_library(qe_toy_lib ${src_qe_lib})
qe_install_targets(qe_toy_lib)
target_link_libraries(qe_toy_lib
	PUBLIC
	qe_fftx
	qe_utilx
	qe_lax
	qe_device_lapack
	PRIVATE
	qe_openacc_fortran)

qe_add_executable(qe_read_test qe_read_test.f90)
set_target_properties(qe_read_test PROPERTIES OUTPUT_NAME qe_read_test.x)
target_link_libraries(qe_read_test PRIVATE qe_toy_lib)

qe_add_executable(qe_hpsi_test qe_hpsi_test.f90)
set_target_properties(qe_hpsi_test PROPERTIES OUTPUT_NAME qe_hpsi_test.x)
target_link_libraries(qe_hpsi_test PRIVATE qe_toy_lib qe_openacc_fortran)

qe_add_executable(qe_davidson qe_davidson_main.f90)
set_target_properties(qe_davidson PROPERTIES OUTPUT_NAME qe_davidson.x)
target_link_libraries(qe_davidson PRIVATE qe_toy_lib qe_kssolver_davidson qe_lax qe_openacc_fortran qe_openmp_fortran)

if(QE_ENABLE_HDF5)
    qe_add_executable(dump_to_h5 tools/dump_to_h5.f90)
    set_target_properties(dump_to_h5 PROPERTIES OUTPUT_NAME dump_to_h5.x)
    target_link_libraries(dump_to_h5 PRIVATE qe_toy_lib)
endif()
