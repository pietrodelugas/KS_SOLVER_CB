set(src_lax_test fortuno.f90)
qe_add_executable (qe_lax_fortuno_test ${src_lax_test}) 
set_target_properties (qe_lax_fortuno_test
	PROPERTIES 
		OUTPUT_NAME qe_lax_fortuno_test.x 
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin) 
target_link_libraries (qe_lax_fortuno_test 
	PRIVATE 
	    qe_mpi_fortran 
	    qe_lax
	    qe_lapack 
	    Fortuno::Fortuno) 
add_unit_test (_qe_laxlib_fortuno_test 1 1 $<TARGET_FILE:qe_lax_fortuno_test>) 



#qe_add_executable(qe_lax_fortuno_lapack_zdotc test_lapack_zdotc.f90)
#set_target_properties(qe_lax_fortuno_lapack_zdotc
#    PROPERTIES
#        OUTPUT_NAME qe_lax_fortuno_lapack_zdotc.x
#        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
#target_link_libraries(qe_lax_fortuno_lapack_zdotc
#    PRIVATE
#        qe_openmp_fortran
#        qe_mpi_fortran
#        qe_lapack
#	Fortuno::Fortuno)
#
#add_unit_test(test_qe_lapack_zdotc 1 1 $<TARGET_FILE:test_qe_lax_lapack_zdotc>)
