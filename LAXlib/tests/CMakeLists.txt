qe_add_executable(test_qe_lax_lapack_zdotc test_lapack_zdotc.f90)
set_target_properties(test_qe_lax_lapack_zdotc
    PROPERTIES
        OUTPUT_NAME test_qe_lax_lapack_zdotc.x
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
target_link_libraries(test_qe_lax_lapack_zdotc
    PRIVATE
        qe_openmp_fortran
        qe_mpi_fortran
        qe_lapack)

add_unit_test(test_qe_lapack_zdotc 1 1 $<TARGET_FILE:test_qe_lax_lapack_zdotc>)

if(QE_ENABLE_CUDA)
    set(src_diaghg_gpu_batched_demo demo_batched/demo_diaghg_gpu_batched.f90)
    qe_enable_cuda_fortran("${src_diaghg_gpu_batched_demo}")
    qe_add_executable(qe_lax_demo_diaghg_gpu_batched ${src_diaghg_gpu_batched_demo})
    set_target_properties(qe_lax_demo_diaghg_gpu_batched
        PROPERTIES
            OUTPUT_NAME demo_diaghg_gpu_batched.x
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
    target_link_libraries(qe_lax_demo_diaghg_gpu_batched
        PRIVATE
            qe_lax
            qe_openmp_fortran)

    add_unit_test(test_qe_lax_diaghg_gpu_batched_demo 1 1 $<TARGET_FILE:qe_lax_demo_diaghg_gpu_batched>)

    if(QE_ENABLE_OPENACC)
        # Same test as above, but driving laxlib_cdiaghg_gpu_batched purely through
        # OpenACC data/host_data mapping instead of explicit CUDA-Fortran DEVICE
        # arrays (still needs the -cuda compile flag, since the "diaghg" interface
        # pulled in from laxlib.fh declares its GPU dummy arguments DEVICE).
        set(src_diaghg_gpu_batched_acc_demo demo_batched/demo_diaghg_gpu_batched_acc.f90)
        qe_enable_cuda_fortran("${src_diaghg_gpu_batched_acc_demo}")
        qe_add_executable(qe_lax_demo_diaghg_gpu_batched_acc ${src_diaghg_gpu_batched_acc_demo})
        set_target_properties(qe_lax_demo_diaghg_gpu_batched_acc
            PROPERTIES
                OUTPUT_NAME demo_diaghg_gpu_batched_acc.x
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
        target_link_libraries(qe_lax_demo_diaghg_gpu_batched_acc
            PRIVATE
                qe_lax
                qe_openacc_fortran
                qe_openmp_fortran)

        add_unit_test(test_qe_lax_diaghg_gpu_batched_acc_demo 1 1 $<TARGET_FILE:qe_lax_demo_diaghg_gpu_batched_acc>)
    endif()
endif()
