# environment for the QE side (QE was built with the 26.3 HPC SDK: module files need the same compiler)
export PATH=/opt/nvidia/hpc_sdk/Linux_aarch64/26.3/compilers/bin:/opt/nvidia/hpc_sdk/Linux_aarch64/26.3/comm_libs/hpcx/bin:$PATH
export QE_CPU=$HOME/q-e/build_cpu_ref/bin/pw.x
export QE_GPU=$HOME/q-e/build_cuda_gb10_gpumpi/bin/pw.x
export QE_PSEUDO=$HOME/q-e/pseudo
# NOTE: the miniapp side (this repo) is built with the 25.9 compiler, see build_marco/CMakeCache.txt
