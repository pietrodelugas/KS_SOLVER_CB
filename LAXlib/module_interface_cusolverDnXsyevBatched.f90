module cusolverxsyevbatched
     use, intrinsic :: iso_c_binding
     use, intrinsic :: iso_fortran_env, only: int64
     use cusolverDn
     IMPLICIT NONE

     ! cudaDataType enumeration constants
     !integer(c_int), parameter :: CUDA_R_32F = 0
     !integer(c_int), parameter :: CUDA_R_64F = 1
     !integer(c_int), parameter :: CUDA_C_32F = 4
     !integer(c_int), parameter :: CUDA_C_64F = 5

     interface
         function cusolverDnXsyevBatched_bufferSize( &
             handle, params, jobz, uplo, n, &
             dataTypeA, dA, lda, &
             dataTypeW, dW, &
             computeType, &
             workspaceInBytesOnDevice, &
             workspaceInBytesOnHost, &
             batch_size) &
             bind(C, name="cusolverDnXsyevBatched_bufferSize")
 
             import :: c_ptr, c_int, c_size_t, int64, cusolverDnHandle, cusolverDnParams
             integer(c_int) :: cusolverDnXsyevBatched_bufferSize
             type(cusolverDnHandle), value :: handle
             type(cusolverDnParams), value :: params
             integer(c_int), value :: jobz
             integer(c_int), value :: uplo
             integer(int64), value :: n
             integer(c_int), value :: dataTypeA
             type(c_ptr), value :: dA
             integer(int64), value :: lda
             integer(c_int), value :: dataTypeW
             type(c_ptr), value :: dW
             integer(c_int), value :: computeType
             integer(c_size_t) :: workspaceInBytesOnDevice
             integer(c_size_t) :: workspaceInBytesOnHost
             integer(int64), value :: batch_size
         end function
         
         function cusolverDnXsyevBatched( &
            handle, params, jobz, uplo, n, &
            dataTypeA, dA, lda, &
            dataTypeW, dW, &
            computeType, &
            workspaceOnDevice, workspaceInBytesOnDevice, &
            workspaceOnHost, workspaceInBytesOnHost, &
            infoArray, batchSize) &
            bind(C, name="cusolverDnXsyevBatched")

            import :: c_ptr, c_int, c_size_t, int64, c_devptr, cusolverDnHandle, cusolverDnParams
            integer(c_int) :: cusolverDnXsyevBatched
            type(cusolverDnHandle), value :: handle
            type(cusolverDnParams), value :: params
            integer(c_int), value :: jobz
            integer(c_int), value :: uplo
            integer(int64), value :: n
            integer(c_int), value :: dataTypeA
            type(c_ptr), value :: dA
            integer(int64), value :: lda
            integer(c_int), value :: dataTypeW
            type(c_ptr), value :: dW
            integer(c_int), value :: computeType
            type(c_devptr), value :: workspaceOnDevice
            integer(c_size_t) :: workspaceInBytesOnDevice
            type(c_ptr), value :: workspaceOnHost
            integer(c_size_t) :: workspaceInBytesOnHost
            integer :: infoArray
            integer(int64), value :: batchSize
        end function

     end interface
end module
