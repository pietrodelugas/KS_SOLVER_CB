###########################################################
# Tests
###########################################################
set(src_fftx_test_with_fortuno test_fft_scalar_with_fortuno.f90)
qe_add_executable(qe_fftx_test_with_fortuno ${src_fftx_test_with_fortuno})
set_target_properties(qe_fftx_test_with_fortuno 
    PROPERTIES 
        OUTPUT_NAME qe_fftx_test_with_fortuno.x 
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
target_link_libraries(qe_fftx_test_with_fortuno
    PRIVATE
        qe_openmp_fortran
        qe_mpi_fortran
	    qe_fftx 
	    Fortuno::Fortuno)
add_unit_test_with_fortuno(qe_fftx_test_with_fortuno 1 1 $<TARGET_FILE:qe_fftx_test_with_fortuno>)



