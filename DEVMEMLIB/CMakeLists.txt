set(src_devmem
   device_memcpy.f90
   device_memcpy_mod.f90
   dxlib_errore.f90)

qe_enable_cuda_fortran("${src_devmem}")
qe_add_library(qe_devmem ${src_devmem})

target_link_libraries(qe_devmem
    PRIVATE
        qe_utilx
        qe_openmp_fortran)
target_include_directories(qe_devmem 
	PRIVATE
	"${CMAKE_CURRENT_SOURCE_DIR}")

qe_install_targets(qe_devmem)


