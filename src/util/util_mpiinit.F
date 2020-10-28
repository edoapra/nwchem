      subroutine util_mpiinit(rank)
      implicit none
#include "mpif.h"      
      integer*4 ierr
      integer*4 rank4,size4
      integer rank              ! [in/out]
#ifdef USE_MPIINITTHREAD
      CALL MPI_Init_thread(MPI_THREAD_MULTIPLE, ierr)
      IF (ierr .LT. MPI_THREAD_MULTIPLE) THEN
         WRITE(*,*) '  MPI_Init_thread_MPI_THREAD_MULTIPLE failed'
         CALL MPI_Abort(MPI_COMM_WORLD, -1)
      endif
#else
      call MPI_init(ierr)
      IF (ierr .ne. 0) THEN
         WRITE(*,*) '  MPI_Init failed'
         CALL MPI_Abort(MPI_COMM_WORLD, -1)
      endif
#endif      
      call MPI_COMM_RANK( MPI_COMM_WORLD, rank4, ierr)
      IF (ierr .ne. 0) THEN
         WRITE(*,*) '  MPI_comm_rank failed'
         CALL MPI_Abort(MPI_COMM_WORLD, -1)
      endif
      rank=rank4
      return
      end
