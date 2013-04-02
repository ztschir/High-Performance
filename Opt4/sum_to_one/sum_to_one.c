#include <stdio.h>
#include "mpi.h"
#include <math.h>

#define N 10    /* Size of the vector */
#define ROOT 0   /* Root of the sum-to-one */

int my_sum_to_one( void * send_buf, void * recv_buf, int count,
                   int root, MPI_Comm comm);

int main( int argc, char **argv )
{
    int me;
    double data_in[N], data_out[N];
    int i;

    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &me );

    for (i=0; i<N; i++) {
      data_in[i] = (double) i + me * pow( 10, -me-1 );
      data_out[i] = (double) -1;
    }

    my_sum_to_one_new( data_in, data_out, N, ROOT, MPI_COMM_WORLD );
    
    if (me == ROOT) {
      printf("Process %d printing result vector:\n", me);
      for (i=0; i<N; i++)
	printf("result[%2d] = %7.5f\n", i, data_out[i]);
    }

    MPI_Finalize( );

    return 0;
}

