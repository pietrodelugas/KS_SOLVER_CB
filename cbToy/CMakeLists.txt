add_subdirectory(cb_lib) 
add_subdirectory(cb_operators)
set(src_cb_davidson
	cb_davidson_main.f90)
qe_add_executable(cb_davidson "${src_cb_davidson}")
set_target_properties(cb_davidson
	PROPERTIES
	    OUTPUT_NAME cb_davidson.x)
target_link_libraries(cb_davidson
	PRIVATE
	cb_lib
	cb_ops
        qe_kssolver_davidson
        qe_lax
	qe_openacc_fortran
	qe_openmp_fortran) 
if (QE_ENABLE_PAPI)
	target_link_libraries(cb_davidson
		PRIVATE
		papi)
endif()

set(src_checker
	my_h_psi_checker.f90)
qe_add_executable(checker  "${src_checker}")
set_target_properties(checker 
	PROPERTIES
	    OUTPUT_NAME checker.x)
target_link_libraries(checker
	PRIVATE
	cb_lib
	cb_ops)

set(src_wrapper
	my_h_psi_wrapper.f90)
qe_add_executable(wrapper  "${src_wrapper}")
set_target_properties(wrapper 
	PROPERTIES
	    OUTPUT_NAME wrapper.x)
target_link_libraries(wrapper
	PRIVATE
	cb_lib
	cb_ops
	qe_openacc_fortran)

