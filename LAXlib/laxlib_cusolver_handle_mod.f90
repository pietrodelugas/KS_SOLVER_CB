MODULE laxlib_cusolver_handles
#if defined(__CUDA) 
  USE cudafor
  USE cusolverdn
  USE cublas !! M.Iovine - added Use cublas for Cublas initialization in the main
  IMPLICIT NONE
  TYPE(cusolverDnHandle),ALLOCATABLE      :: cusolver_handle(:)
  LOGICAL,ALLOCATABLE                     :: cusolver_initialized(:)
  TYPE(cublasHandle),ALLOCATABLE          :: cublas_handle(:) !!M.Iovine - we add the variable for the cublas handle
  LOGICAL,ALLOCATABLE                     :: cublas_initialized(:) !!M.Iovine - we add the boolean for the cublas handle
  !
  LOGICAL, SAVE                :: cusolver_initialized_host = .FALSE.
  LOGICAL, SAVE                :: cublas_initialized_host = .FALSE. !!M.Iovine - added variable for cublas
  INTEGER, SAVE                ::  cusolver_thread = 0
  INTEGER(kind=cuda_stream_kind) :: laxlib_cuda_stream = 0
  !$omp threadprivate(cusolver_thread, laxlib_cuda_stream)
  INTEGER,SAVE                     ::  cusolver_thread_max
  INTEGER,SAVE                     ::  cublas_thread_max !!M.Iovine - added variable for cublas
  PUBLIC :: initialize_cusolver_handles, get_cusolver_handle, get_cusolver_initialized, &
            set_laxlib_cuda_stream, finalize_cusolver_handles, cusolver_handle, cusolver_initialized,&
            laxlib_cuda_stream, initialize_laxlib_cuda_stream, &
            initialize_cublas_handles, get_cublas_handle, get_cublas_initialized, finalize_cublas_handles, & !!M.Iovine- cublas adds
            cublas_handle, cublas_initialized, finalize_cublas_handles, cublas_handle, cublas_initialized
  CONTAINS 
    SUBROUTINE initialize_cusolver_handles(nthreads)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nthreads
      IF ( cusolver_initialized_host ) THEN
         CALL errore( ' initialize_cusolver_handles ', 'Cusolver handles already initialized',  ABS( cusolver_thread ) )   
      END IF
      ALLOCATE(cusolver_handle(nthreads),cusolver_initialized(nthreads))
      cusolver_thread_max = nthreads
        cusolver_initialized(:) = .FALSE.
        cusolver_initialized_host = .TRUE.
    END SUBROUTINE initialize_cusolver_handles


    SUBROUTINE initialize_laxlib_cuda_stream( stream, mypippo)
       IMPLICIT NONE
       INTEGER(cuda_stream_kind), INTENT(IN) :: stream
       INTEGER, INTENT(IN) :: mypippo
       cusolver_thread = mypippo
       laxlib_cuda_stream = stream
       print '("In initialize_laxlib_cuda_stream, thread ",I5,I24)', cusolver_thread, laxlib_cuda_stream
    END SUBROUTINE initialize_laxlib_cuda_stream

    SUBROUTINE finalize_cusolver_handles()
      IMPLICIT NONE
      INTEGER :: i, info
      IF ( cusolver_initialized_host ) THEN
         DO i = 1, cusolver_thread_max
            IF ( cusolver_initialized(i) ) THEN
               info = cusolverDnDestroy(cusolver_handle(i))
               IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' finalize_cusolver_handles ', 'cusolverDnDestroy',  ABS( info ) )
               cusolver_initialized(i) = .FALSE.
            END IF
         END DO
         DEALLOCATE(cusolver_handle, cusolver_initialized)
         cusolver_initialized_host = .FALSE.
      ENDIF
    END SUBROUTINE finalize_cusolver_handles

    SUBROUTINE get_cusolver_handle( mythread, handle )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mythread
      TYPE(cusolverDnHandle), INTENT(OUT) :: handle
      IF ( mythread < 1 .OR. mythread > cusolver_thread_max ) THEN
         CALL lax_error__( ' cusolver_handle ', 'Invalid thread index', ABS( mythread ) )
      ENDIF
      handle = cusolver_handle(mythread)
    END SUBROUTINE get_cusolver_handle

    SUBROUTINE get_cusolver_initialized( mythread, initialized )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mythread
      LOGICAL, INTENT(OUT) :: initialized
      IF ( mythread < 1 .OR. mythread > cusolver_thread_max ) THEN
         CALL lax_error__( ' cusolver_initialized ', 'Invalid thread index', ABS( mythread ) )
      ENDIF
      initialized = cusolver_initialized(mythread)
    END SUBROUTINE get_cusolver_initialized

    SUBROUTINE set_laxlib_cuda_stream( stream )
      IMPLICIT NONE
      INTEGER(cuda_stream_kind), INTENT(IN) :: stream
      laxlib_cuda_stream = stream
    END SUBROUTINE set_laxlib_cuda_stream
    
    !!!! M.Iovine - created subroutines for cublas :
    SUBROUTINE initialize_cublas_handles(nthreads)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nthreads
      IF ( cublas_initialized_host ) THEN
         CALL errore( ' initialize_cublas_handles ', 'Cublas handles already initialized',  ABS( cusolver_thread ) )
      END IF
      ALLOCATE(cublas_handle(nthreads),cublas_initialized(nthreads))
      cublas_thread_max = nthreads
        cublas_initialized(:) = .FALSE.
        cublas_initialized_host = .TRUE.
    END SUBROUTINE initialize_cublas_handles

    SUBROUTINE finalize_cublas_handles()
      IMPLICIT NONE
      INTEGER :: i, info
      IF ( cublas_initialized_host ) THEN
         DO i = 1, cublas_thread_max
            IF ( cublas_initialized(i) ) THEN
               info = cublasDestroy(cublas_handle(i))
               IF ( info /= CUSOLVER_STATUS_SUCCESS ) CALL lax_error__( ' finalize_cublas_handles ', 'cublasDestroy',  ABS( info ) )
               cublas_initialized(i) = .FALSE.
            END IF
         END DO
         DEALLOCATE(cublas_handle, cublas_initialized)
         cublas_initialized_host = .FALSE.
      ENDIF
    END SUBROUTINE finalize_cublas_handles

    SUBROUTINE get_cublas_handle( mythread, handle )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mythread
      TYPE(cublasHandle), INTENT(OUT) :: handle
      IF ( mythread < 1 .OR. mythread > cublas_thread_max ) THEN
         CALL lax_error__( ' cublas_handle ', 'Invalid thread index', ABS( mythread ) )
      ENDIF
      handle = cublas_handle(mythread)
    END SUBROUTINE get_cublas_handle

    SUBROUTINE get_cublas_initialized( mythread, initialized )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mythread
      LOGICAL, INTENT(OUT) :: initialized
      IF ( mythread < 1 .OR. mythread > cublas_thread_max ) THEN
         CALL lax_error__( ' cublas_initialized ', 'Invalid thread index', ABS( mythread ) )
      ENDIF
      initialized = cublas_initialized(mythread)
    END SUBROUTINE get_cublas_initialized
    !!!!!!

#endif
END MODULE laxlib_cusolver_handles
