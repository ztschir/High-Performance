#include "mpi.h"

void my_sum_to_one( void * send_buf, void * recv_buf, int count, 
		   int root, MPI_Comm comm)
{
  MPI_Reduce( send_buf, recv_buf, count, MPI_DOUBLE, MPI_SUM, root, comm );
}


