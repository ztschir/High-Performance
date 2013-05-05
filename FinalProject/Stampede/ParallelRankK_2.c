#include <malloc.h>
#include "mpi.h"
#include "prototypes.h"

#define LOCAL_A( i, j ) local_A[ (j)*local_ldimA + (i) ]
#define LOCAL_B( i, j ) local_B[ (j)*local_ldimB + (i) ]

void ParallelRankK( 
		   // global matrix dimensions m, n, k
		   int global_m, int global_n, int global_k,  
		   // the global index of the current column of A and current row of B 
		   // to be used for the rank-1 update
		   int p,
		   int bk, //the block size for this rank-k update
		   // address of where local A, B, and C are stored, as well as their 
		   // local "leading dimensions"
		   double *local_A, int local_ldimA, 
		   double *local_B, int local_ldimB, 
		   double *local_C, int local_ldimC,
		   // communicator for the row in which this node is a member
		   MPI_Comm comm_row, 
		   // communicator for the column in which this node is a member
		   MPI_Comm comm_col )
{
  /* This is a utility routine for implementing the matrix-matrix multiplication
     C = A B + C.  The idea is that C is updated by a rank-1 update with the p-th column of A 
     and p-th row of B.  A loop around a call to this routine will then implement the 
     matrix-matrix multiplication.
  */

  int 
    nrows, ncols,     // number of rows and columns in the mesh of nodes
    myrow, mycol,     // this node's row and column in the mesh of nodes
    local_m, local_n, // the number of rows and columns of C assigned to this node
    currow, curcol,   // index of the node row and node column in which the p-th row of B
                      // and p-th column of A exist
    local_p_A,        // on node that owns p-th column of A, this is the index of that local column
    local_p_B,        // on node that owns p-th row of B, this is the index of that local row
    i, j,             // indices over rows and columns
    i_one = 1;        // integer one, so we can pass by address
  double 
    *work_A,          // address of work buffer in which the copy of current column of A will be put
    *work_B,          // address of work buffer in which the copy of current row    of B will be put
    d_one = 1.0;      // double precision one, so we can pass by address
  // extract information about the mesh of nodes
  MPI_Comm_size( comm_col, &nrows );
  MPI_Comm_size( comm_row, &ncols );
  MPI_Comm_rank( comm_col, &myrow );
  MPI_Comm_rank( comm_row, &mycol );




  // Compute the local row and column dimension of C, A, and B
  local_m   = global_m / nrows + ( myrow < global_m % nrows ? 1 : 0 );
  local_n   = global_n / ncols + ( mycol < global_n % ncols ? 1 : 0 );

  //  local_k_A = global_k / ncols + ( mycol < global_k % ncols ? 1 : 0 );
  //  local_k_B = global_k / nrows + ( myrow < global_k % nrows ? 1 : 0 );

  // Create a work array into which to receive the current column of A when
  // broadcast within rows
  work_A = ( double *) malloc( sizeof( double ) * local_m *bk);

  // Create a work array into which to receive the current row of B when
  // broadcast within columns
  work_B = ( double * ) malloc( sizeof( double ) * local_n *bk);

int x, offset_a, offset_b;
for(x=p; p<x+bk && p<global_k ; p++)
{
  // Compute which column of nodes owns the p-th column of A
  curcol = p % ncols;
  offset_a =  local_m * (p-x);


  // On the column of nodes that owns the p-th column of A, compute which local
  // column of A we need and pack that column into the work array.
  if ( curcol == mycol ){
    local_p_A = p / ncols;
    for ( i=0; i<local_m; i++ ){
      work_A[ i + offset_a] = LOCAL_A( i, local_p_A );
//	printf("P = %d    (Row, Col) =  (%d, %d) has %d rows and %d cols of C value = %d \n", p, myrow, mycol, local_m,
//		local_n, work_A[i+offset]);

     }
  }




  // Broadcast the copy of the local column of A within this node's row of nodes
  MPI_Bcast( &(work_A[offset_a]), local_m, MPI_DOUBLE, curcol, comm_row );
}
p=x;
for(x=p; p<x+bk && p<global_k ; p++)
{
  // Compute which row of nodes owns the p-th row of B
  currow = p % nrows;
  offset_b = local_n*(p-x);

  // On the row of nodes that owns the p-th row of B, compute which local
  // row of B we need and pack that row into the work array.
  if ( currow == myrow ){
    local_p_B = p / nrows;
    for ( j=0; j<local_n; j++ ){
      work_B[ j + offset_b] = LOCAL_B( local_p_B, j );
//	   printf("P = %d    (Row, Col) =  (%d, %d) has %d rows and %d cols of C value = %d \n", p, myrow, mycol, local_m,
//               local_n, work_B[j + offset]);
     }
  }

  // Broadcast the copy of the local row of B within this node's column of nodes

  MPI_Bcast( &(work_B[offset_b]), local_n, MPI_DOUBLE, currow, comm_col );

}
p=x;
for(x=p; p<x+bk && p<global_k ; p++)
{
	offset_a = local_m*(p-x);
	offset_b = local_n*(p-x);
  // Perform local rank-1 update
  dger_( &local_m, &local_n, 
	 &d_one, &(work_A[offset_a]), &i_one, &(work_B[offset_b]), &i_one, 
	 local_C, &local_ldimC );
}
  free( work_A );
  free( work_B );


}
