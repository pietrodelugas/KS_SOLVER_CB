#!/bin/bash
# usage: build_qe_tool.sh source.f90 output.x   (compiles and links a program against the QE CPU reference build)
export PATH=/opt/nvidia/hpc_sdk/Linux_aarch64/26.3/compilers/bin:/opt/nvidia/hpc_sdk/Linux_aarch64/26.3/comm_libs/hpcx/bin:$PATH
QB=${QE_BUILD:-$HOME/q-e/build_cpu_ref}
OMPI=/opt/nvidia/hpc_sdk/Linux_aarch64/26.3/comm_libs/13.1/hpcx/hpcx-2.25.1/ompi/lib
INC="-I$HOME/q-e/include -I/home/qespresso/q-e/build_cpu_ref/dft-d3/mod/qe_dftd3 -I/home/qespresso/q-e/build_cpu_ref/external/mbd/src/modules -I/home/qespresso/q-e/build_cpu_ref/external/mod/qe_devxlib -I/home/qespresso/q-e/build_cpu_ref/FFTXlib/src/mod/qe_fftx -I/home/qespresso/q-e/build_cpu_ref/KS_Solvers/mod/qe_kssolver_davidson -I/home/qespresso/q-e/build_cpu_ref/LAXlib/mod/qe_lax -I/home/qespresso/q-e/build_cpu_ref/Modules/mod/qe_modules -I/home/qespresso/q-e/build_cpu_ref/PW/mod/qe_pw -I/home/qespresso/q-e/build_cpu_ref/upflib/mod/qe_upflib -I/home/qespresso/q-e/build_cpu_ref/upflib/mod/qe_xml -I/home/qespresso/q-e/build_cpu_ref/UtilXlib/mod/qe_utilx -I/home/qespresso/q-e/build_cpu_ref/XClib/mod/qe_xclib -I$QB/LAXlib/include -I$HOME/q-e/LAXlib"
L=$QB/lib
mpif90 -DOMPI_SKIP_MPICXX -D__FFTW3 -D__MPI -D__MPI_MODULE -D__PGI -Mcache_align -Mlarge_arrays -Mbackslash -fast -O3 -cpp $INC "$1" \
  -Wl,-rpath -Wl,$OMPI -Wl,--enable-new-dtags -o "$2" \
  $L/libqe_pw.a $L/libqe_modules.a $L/libqe_lax.a $L/libqe_upflib.a $L/libqe_modules_c.a $L/libqe_xml.a $L/libqe_kssolver_davidson.a \
  $L/libqe_kssolver_cg.a $L/libqe_kssolver_paro.a $L/libqe_device_lapack.a $L/libqe_kssolver_rmmdiis.a $L/libqe_kssolver_dense.a \
  $L/libqe_kssolver_direct.a $L/libqe_lax.a $L/libqe_fftx.a $L/libqe_fftx_c.a /usr/lib/aarch64-linux-gnu/libfftw3.so $L/libqe_dftd3.a \
  $L/libqe_utilx.a $L/libqe_utilx_c.a $L/libqe_xclib.a $OMPI/libmpi_usempif08.so $OMPI/libmpi_usempi_ignore_tkr.so $OMPI/libmpi_mpifh.so \
  $OMPI/libmpi.so $L/libqe_libbeef.a $L/libqe_devxlib.a $L/libmbd.a /opt/nvidia/hpc_sdk/Linux_aarch64/26.3/compilers/lib/liblapack.so \
  /opt/nvidia/hpc_sdk/Linux_aarch64/26.3/compilers/lib/libblas.so
