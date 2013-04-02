#include "mpi.h"

void my_sum_to_one_util(void *, void *, int, int, MPI_Comm, int, int );

void my_sum_to_one_new( void * send_buf, void * recv_buf, int count, 
		   int root, MPI_Comm comm)
{
  int nprocs;

  MPI_Comm_size ( comm, &nprocs );
  my_sum_to_one_util(send_buf, recv_buf, count, root, comm, 0, nprocs-1);

  //  MPI_Reduce( send_buf, recv_buf, count, MPI_DOUBLE, MPI_SUM, root, comm );
}


void my_sum_to_one_util(void * send_buf, void * recv_buf, int count,
			int cur_root, MPI_Comm comm, int left, int right )
{
  int me, mid, dest;
  MPI_Status status;
  
  if ( left == right ) return;

  MPI_Comm_rank( comm, &me );

  mid = ( left + right )/2;

  if ( cur_root <= mid )
    dest = mid + 1;
  else
    dest = left;

  if ( cur_root <= mid ){
    if ( me <= mid )
      my_sum_to_one_util( send_buf, recv_buf, count, cur_root, comm, left, mid );
    else		     
      my_sum_to_one_util( send_buf, recv_buf, count, dest, comm, mid+1, right );
  }
  else{
    if ( me <= mid )
      my_sum_to_one_util( send_buf, recv_buf, count, dest, comm, left, mid );
    else		     
      my_sum_to_one_util( send_buf, recv_buf, count, cur_root, comm, mid+1, right );
  }

  if ( me == cur_root ){
    MPI_Recv( recv_buf, count, MPI_DOUBLE, dest, MPI_ANY_TAG, comm , &status );
  }
  if ( me == dest ){
    int i;
    for (i=0; i<count; i++ )
      ( ( double * ) send_buf )[ i ] += ( ( double * ) recv_buf )[ i ]; 
    MPI_Send( send_buf, count, MPI_DOUBLE, cur_root, 0, comm );
  }
  return;
}
