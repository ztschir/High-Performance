#include "mpi.h"
#include "prototypes.h"

#define GLOBAL_A( i, j ) global_A[ (j)*global_ldimA + (i) ]
#define LOCAL_A( i, j )  local_A[  (j)*local_ldimA  + (i) ]

void CopyMatrixGlobalToLocal(
		   // global matrix dimensions m, n
		   int global_m, int global_n, 
		   // address of where global, A and local A are stored, as well as their 
		   // local "leading dimensions"
		   double *global_A, int global_ldimA, 
		   double *local_A, int local_ldimA, 
		   // communicator for the row in which this node is a member
		   MPI_Comm comm_row, 
		   // communicator for the column in which this node is a member
		   MPI_Comm comm_col )
{
  /* This is a utility routine for copying a matrix A that is a complete copy on each node 
     into the array that will hold the local part of matrix A.
  */

  int 
    nprows, npcols,   // number of rows and columns in the mesh of nodes
    myrow, mycol,     // this node's row and column in the mesh of nodes
    i, j,             // indices for striding through the global A
    local_i, local_j; // indices for striding through the local A

  // extract information about the mesh of nodes
  MPI_Comm_size( comm_col, &nprows );
  MPI_Comm_size( comm_row, &npcols );
  MPI_Comm_rank( comm_col, &myrow );
  MPI_Comm_rank( comm_row, &mycol );


  // Copy the part that this node locally owns into local_A from global_A
  // Note: the first element (top-left element) of global_A that this node owns is
  // GLOBAL_A( myrow, mycol ).  From there, in the column (n) direction, every npcols-th column 
  // is assigned to this node, and in the row (m) direction, every nprows-th row is assigned.
  local_j = 0;
  for ( j=mycol; j<global_n; j+=npcols ){
    local_i = 0;
    for ( i=myrow; i<global_m; i+=nprows ){
      LOCAL_A( local_i, local_j ) = GLOBAL_A( i, j );
      local_i++;
    }
    local_j++;
  }
}
	

