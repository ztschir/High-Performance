/* Create macros so that the matrices are stored in column-major order */

#include "parameters.h"

#define A(i,j) a[ (j)*lda + (i) ]
#define AT(i,j) at[ (j)*(k) + (i) ]
#define B(i,j) b[ (j)*ldb + (i) ]
#define C(i,j) c[ (j)*ldc + (i) ]

#define min( x, y ) ( x < y ? x : y )

#define MB 144
#define KB 144

/* MY_MMult_inner prototype */
void MY_MMult_inner( int m, int n, int k, double *a, int lda, 
                                    double *b, int ldb,
                                    double *c, int ldc );
                                                                        

/* Routine for computing C = A * B + C */

void MY_MMult( int m, int n, int k, double *a, int lda, 
                                    double *b, int ldb,
                                    double *c, int ldc )
{
  int i, p, ib, pb;

  for ( p=0; p<k; p+=KB ){
    pb = min( k-p, KB );
    for ( i=0; i<m; i+=MB ){
      ib = min( m-i, MB );

      MY_MMult_inner( ib, n, pb, 
                      &A( i, p ), lda, 
                      &B( p, 0 ), ldb, 
                      &C( i, 0 ), ldc );
    }
  }
}


  


/* Inner kernel for computing C = A * B + C */

void MY_MMult_inner( int m, int n, int k, double *a, int lda, 
                                    double *b, int ldb,
                                    double *c, int ldc )
{
  int i, j, p;
  double at[ MB * KB ];

  transpose_matrix( m, k, a, lda, at, k );

  for ( j=0; j<n; j++ ){
    double *cp = &C( 0, j );

    for ( i=0; i<m; i+=4 ){
      register double c0=0.0, c1=0.0, c2=0.0, c3=0.0;
      double *a0p = &AT(0,i  ), *a1p = &AT(0,i+1), 
             *a2p = &AT(0,i+2), *a3p = &AT(0,i+3), 
             *bp=&B(0,j);

      for ( p=0; p<k; p++ ){
        register double bpj = *bp;

        c0 = c0 +  *a0p++ * bpj;
        c1 = c1 +  *a1p++ * bpj;
        c2 = c2 +  *a2p++ * bpj;
        c3 = c3 +  *a3p++ * bpj;

        bp++;
      }
      
      cp[ 0 ] += c0;
      cp[ 1 ] += c1;
      cp[ 2 ] += c2;
      cp[ 3 ] += c3;
      
      cp += 4;
    }
  }
}
