#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

// Prototypes
#include "prototypes.h"

int main(int argc, char **argv )
{
  int 
    me,              /* holds the index of "this" process */
    nprocs,          /* holds the number of processes involved */
    nprows, npcols,  /* mesh sizes */
    myrow, mycol,    /* this node's mesh coordinates */
    n,               /* global matrix size */
    local_m,         /* local row size of A, B, C */
    local_n;         /* local column size of A, B, C */

  double
    *global_A,       /* array in which to hold matrix A */
    *global_B,       /* array in which to hold matrix B */
    *global_C,       /* array in which to hold matrix C */
    *local_A,        /* array in which to hold local part of matrix A */
    *local_B,        /* array in which to hold local part of matrix B */
    *local_C,        /* array in which to hold local part of matrix C */
    *local_C_ref,    /* array in which to hold local part of matrix C */
    local_diff,      /* hold difference between sequential and */
    diff,            /*			parallel result */
    d_one = 1.0;     /* double precision one, to pass by address */

  MPI_Comm
    comm_row, comm_col;  /* communicators for the row and col in which this 
			    node exists */

  /* Initialize MPI, passing in the command-line parameters.  MPI_Init
     strips out the command-line parameters for MPI (e.g., the -np 5
     in our example) and then returns argc and argv with without those
     parameters. */
  MPI_Init( &argc, &argv );

  /* Inquire how many processes were started up by mpiexec */
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

  /* Inquire the index (rank) of "this" process within the proceses
     started up by mpiexec */
  MPI_Comm_rank( MPI_COMM_WORLD, &me );

  /* Process 0 accepts input */
  if ( me == 0 ){
    printf("enter matrix size n:");
    scanf( "%d", &n );
    printf("enter nprows, npcols:");
    scanf( "%d%d", &nprows, &npcols );
  }

  /* share parameters  with all nodes */
  MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &nprows, 1, MPI_INT, 0, MPI_COMM_WORLD );
  MPI_Bcast( &npcols, 1, MPI_INT, 0, MPI_COMM_WORLD );

  if ( nprows * npcols != nprocs ){
    printf( "mesh not of right size\n" );
    exit( 0 );
  }
  
  /* Figure out what my index is */
  mycol = me / nprows;
  myrow = me % nprows;

  /* create a communicator for the row of which I am part */
  MPI_Comm_split( MPI_COMM_WORLD, myrow, mycol, &comm_row );

  /* create a communicator for the column of which I am part */
  MPI_Comm_split( MPI_COMM_WORLD, mycol, myrow, &comm_col );

  /* create buffers into which to hold the global A, B, C (everyone will have a copy) */
  global_A = ( double * ) malloc ( sizeof( double ) * n * n );
  global_B = ( double * ) malloc ( sizeof( double ) * n * n );
  global_C = ( double * ) malloc ( sizeof( double ) * n * n );

  /* create random matrices on node zero and share with all nodes */
  if ( me == 0 ){
    random_matrix( n, n, global_A, n );
    random_matrix( n, n, global_B, n );
    random_matrix( n, n, global_C, n );
  }
  MPI_Bcast( global_A, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( global_B, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD );
  MPI_Bcast( global_C, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  /* compute local matrix sizes */
  local_m   = n / nprows + ( myrow < n % nprows ? 1 : 0 );
  local_n   = n / npcols + ( mycol < n % npcols ? 1 : 0 );

  /* create buffer into which to hold the ocal A, B, C */
  local_A = ( double * ) malloc ( sizeof( double ) * local_m * local_n );
  local_B = ( double * ) malloc ( sizeof( double ) * local_m * local_n );
  local_C = ( double * ) malloc ( sizeof( double ) * local_m * local_n );

  /* copy the local parts */
  CopyMatrixGlobalToLocal( n, n, 
			   global_A, n, 
			   local_A, local_m, 
			   comm_row, comm_col );
  CopyMatrixGlobalToLocal( n, n, 
			   global_B, n, 
			   local_B, local_m, 
			   comm_row, comm_col );
  CopyMatrixGlobalToLocal( n, n, 
			   global_C, n, 
			   local_C, local_m, 
			   comm_row, comm_col );

  /* Compute parallel matrix-matrix multiply */
  ParallelMMult( n, n, n, 
		 local_A, local_m, 
		 local_B, local_m, 
		 local_C, local_m, 
		 comm_row, comm_col );

  /* Compute sequential matrix-matrix multiply on all nodes */
  dgemm_( "N", "N", &n, &n, &n,
  	  &d_one, global_A, &n, global_B, &n,
  	  &d_one, global_C, &n );

  local_C_ref = ( double * ) malloc ( sizeof( double ) * local_m * local_n );

  CopyMatrixGlobalToLocal( n, n, 
			   global_C, n, 
			   local_C_ref, local_m, 
			   comm_row, comm_col );

  local_diff = compare_matrices( local_m, local_n, 
				 local_C, local_m, 
				 local_C_ref, local_m );

  MPI_Allreduce( &local_diff, &diff, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

  if ( me == 0 )
    printf("\ndiff = %le\n", diff );

  free( global_A );
  free( global_B );
  free( global_C );
  free( local_A );
  free( local_B );
  free( local_C );
  free( local_C_ref);

  MPI_Comm_free( &comm_row );
  MPI_Comm_free( &comm_col );

  /* Cleanup up the MPI environment */
  MPI_Finalize();

  exit( 0 );
}

