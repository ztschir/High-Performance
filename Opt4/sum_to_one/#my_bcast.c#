#include "mpi.h"

void my_bcast_util( void *, int, MPI_Datatype, int, MPI_Comm, int, int );


void my_bcast( void *buffer, int count, MPI_Datatype datatype,
	       int root, MPI_Comm comm )
{
  int nprocs;

  MPI_Comm_size( comm, &nprocs );

  my_bcast_util( buffer, count, datatype, root, comm, 0, nprocs-1 );

  return;
}



void my_bcast_util( void *buffer, int count, MPI_Datatype datatype,
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

  if ( me == cur_root )
    MPI_Send( buffer, count, datatype, dest, 0, comm );

  if ( me == dest )
    MPI_Recv( buffer, count, datatype, cur_root, MPI_ANY_TAG, comm,
	      &status );

  if ( cur_root <= mid ){
    if ( me <= mid )
      my_bcast_util( buffer, count, datatype, cur_root, comm, left, mid );
    else		     
      my_bcast_util( buffer, count, datatype, dest, comm, mid+1, right );
  }
  else{
    if ( me <= mid )
      my_bcast_util( buffer, count, datatype, dest, comm, left, mid );
    else		     
      my_bcast_util( buffer, count, datatype, cur_root, comm, mid+1, right );
  }

  return;
}

    
