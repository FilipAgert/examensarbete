! example.f90
program example
    use, intrinsic :: omp_lib
    implicit none
    
    ! Set number of threads to use.
    call omp_set_num_threads(4)

    !$omp parallel

        print '("Thread: ", i0)', omp_get_thread_num()

    !$omp end parallel
end program example