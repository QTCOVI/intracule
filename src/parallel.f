      subroutine init_parallel() 
        include 'implicit.inc'
        include 'error.inc'
        include 'parallel.inc'
        ntasks=1
        task_id=0
        mp_mpi=.true.
        nthreads=1
        thread_id=0
!$omp parallel
!$omp master
!$      nthreads=omp_get_num_threads()
!$omp end master
!$      thread_id=omp_get_thread_num()
!$omp end parallel
      end subroutine init_parallel

      subroutine end_parallel()
        include 'implicit.inc'
        include 'error.inc'
        include 'parallel.inc'
      end subroutine end_parallel

      subroutine info_parallel()
        include 'implicit.inc'
        include 'io.inc'
        include 'parallel.inc'
        write(uout,'(A,I3)') " # Number of OMP threads : ", nthreads
        write(uout,'(A,I3)') " # Number of MPI process : ", ntasks
      end subroutine info_parallel

      function is_iomaster()

        include 'implicit.inc'
        include 'parallel.inc'
        logical :: is_iomaster
 
        is_iomaster = (task_id == master)
 
      end function is_iomaster
