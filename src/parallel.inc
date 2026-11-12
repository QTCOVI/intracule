!.....Scratch dir
!     character(len=132) :: scratch
!.....Parameters
      integer(kind=4), parameter :: master=0
!.....OMP variables
      include 'omp_lib.h'
      integer(kind=4) :: ichunk, nthreads, thread_id, stack
      common /gomp_parallel/ ichunk, nthreads, stack
      common /pomp_parallel/ thread_id
!.....MPI variables
!     include 'mpif.h'
      integer(kind=4) :: ntasks, task_id
      logical :: mp_mpi 
      common /gmpi_parallel/ ntasks, mp_mpi
      common /pmpi_parallel/ task_id
!
!$omp threadprivate(/pomp_parallel/)
