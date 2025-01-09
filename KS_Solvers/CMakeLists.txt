set(ks_headers
  ks_solver_interfaces.h)
set (ks_include_dir ${CMAKE_CURRENT_BINARY_DIR}/include)
#FIXME maybe better to move this part to an helper function 
foreach (in_h ${ks_headers})
  get_filename_component(in_h_basename ${in_h} NAME_WE)
  set (out_h "${ks_include_dir}/${in_h_basename}.fh")
  qe_preprocess_source("${CMAKE_CURRENT_SOURCE_DIR}/${in_h}" ${out_h})
  list (APPEND ks_out_headers ${out_h})
endforeach()
add_custom_target(qe_ks_headers 
  DEPENDS ${ks_out_headers}
  VERBATIM )

set(src_davidson
    Davidson/cegterg.f90
    Davidson/regterg.f90
    # GPU
    #Davidson/regterg_gpu.f90
    #Davidson/cegterg_gpu.f90
    )
qe_enable_cuda_fortran("${src_davidson}")

set(src_dense
    DENSE/rotate_HSpsi_gamma.f90
    DENSE/rotate_HSpsi_k.f90
    DENSE/rotate_wfc_gamma.f90
    DENSE/rotate_wfc_k.f90 
    DENSE/rotate_xpsi_k.f90
    DENSE/rotate_xpsi_gamma.f90
    DENSE/gram_schmidt_k.f90
    DENSE/gram_schmidt_gamma.f90
    DENSE/rotate_driver.f90
    # GPU
    DENSE/rotate_wfc_gamma_gpu.f90
    DENSE/rotate_xpsi_k_gpu.f90
    DENSE/rotate_xpsi_gamma_gpu.f90
    DENSE/gram_schmidt_k_gpu.f90
    DENSE/gram_schmidt_gamma_gpu.f90
    DENSE/rotate_wfc_k_gpu.f90
    DENSE/rotate_driver_cuf.f90)
qe_enable_cuda_fortran("${src_dense}")

qe_add_library(qe_kssolver_davidson ${src_davidson})
target_link_libraries(qe_kssolver_davidson
    PRIVATE
        qe_device_lapack
        qe_openacc_fortran
        qe_openmp_fortran
        qe_lax
        qe_utilx
        qe_mpi_fortran
        qe_devmem)

qe_add_library(qe_kssolver_dense ${src_dense})
target_link_libraries(qe_kssolver_dense
    PRIVATE
        qe_openacc_fortran
        qe_lax
        qe_utilx
        qe_mpi_fortran
        qe_devmem)
target_include_directories(qe_kssolver_dense
  PUBLIC
    $<BUILD_INTERFACE:${ks_include_dir}>
    $<INSTALL_INTERFACE:include/qe>
)
set_target_properties(qe_kssolver_dense
  PROPERTIES PUBLIC_HEADER ${ks_out_headers}
)
add_dependencies(qe_kssolver_dense qe_ks_headers)

###########################################################

qe_install_targets(
    qe_kssolver_davidson
    qe_kssolver_dense)
