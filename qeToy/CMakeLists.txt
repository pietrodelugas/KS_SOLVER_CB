set(src_qe_lib
	qe_lib/qe_ham_module.f90
	qe_lib/qe_batch_module.f90
	qe_operators/qe_h_psi_batched_iter.f90)
qe_enable_cuda_fortran("${src_qe_lib}")
qe_add_library(qe_toy_lib ${src_qe_lib})
qe_install_targets(qe_toy_lib)
target_link_libraries(qe_toy_lib
	PUBLIC
	qe_fftx
	qe_utilx
	qe_lax
	PRIVATE
	qe_openacc_fortran)

qe_add_executable(qe_read_test qe_read_test.f90)
set_target_properties(qe_read_test PROPERTIES OUTPUT_NAME qe_read_test.x)
target_link_libraries(qe_read_test PRIVATE qe_toy_lib)

qe_add_executable(qe_hpsi_test qe_hpsi_test.f90)
set_target_properties(qe_hpsi_test PROPERTIES OUTPUT_NAME qe_hpsi_test.x)
target_link_libraries(qe_hpsi_test PRIVATE qe_toy_lib qe_openacc_fortran)
