#include <math.h>
#include "FLAME.h"
#include "LU_prototypes.h"

#define AA(i,j) buff_A[ (j)*ldim_A + (i) ]

FLA_Error REF_LU( int time_lapack, FLA_Obj A, int nb_alg )
{
  int n, ldim_A, info;
  
  double *buff_A, sqrt();


  n = FLA_Obj_length( A );
  ldim_A = FLA_Obj_col_stride( A );

  buff_A = (double *) FLA_Obj_buffer_at_view( A );

  {
    int i, j, k;

    for ( j=0; j<n; j++ ){

      /* a21 = a21 / alpha11 */
      for ( i=j+1; i<n; i++ )
	AA( i,j ) = AA( i,j )  / AA( j,j );

      /* A22 = A22 - a21 * a12t */
      for ( k=j+1; k<n; k++ )
	for ( i=j+1; i<n; i++ )
	  AA( i,k ) = AA( i,k ) - AA( i,j ) * AA( j,k );
    }
  }

  return FLA_SUCCESS;
}

