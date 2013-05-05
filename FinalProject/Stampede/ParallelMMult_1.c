#include <malloc.h>
#include "mpi.h"
#include "prototypes.h"

void ParallelMMult( 
		   // global matrix dimensions m, n, k
		   int global_m, int global_n, int global_k,  
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
  /* This routine implements the matrix-matrix multiplication
     C = A B + C.  The idea is that C is updated by  rank-1 updates with the p-th column of A 
     and p-th row of B, for all p.  
  */

  int 
    p;                // the global index of the current column of A and current row of B 
		      // to be used for the rank-1 update

  for ( p = 0; p<global_k; p++ ){
     ParallelRank1( global_m, global_n, global_k, 
		    p,
		    local_A, local_ldimA,
		    local_B, local_ldimB,
		    local_C, local_ldimC,
		    comm_row, comm_col );
  }
}
