/* Create macros so that the matrices are stored in column-major order */

#define A(i,j) a[ (j)*lda + (i) ]
#define B(i,j) b[ (j)*ldb + (i) ]
#define C(i,j) c[ (j)*ldc + (i) ]

/* Routine for computing C = A * B + C */

void MY_MMult( int m, int n, int k, double *a, int lda, 
                                    double *b, int ldb,
                                    double *c, int ldc )
{
  int i, j, p;

  for ( j=0; j<n; j++ ){
    for ( i=0; i<m; i+=4 ){
      register double c0=0.0, c1=0.0, c2=0.0, c3=0.0;

      for ( p=0; p<k; p++ ){
        c0 += A( i,p ) * B( p,j );
        c1 += A( i+1,p ) * B( p,j );
        c2 += A( i+2,p ) * B( p,j );
        c3 += A( i+3,p ) * B( p,j );
      }
      
      C( i,j )   += c0;
      C( i+1,j ) += c1;
      C( i+2,j ) += c2;
      C( i+3,j ) += c3;
    }
  }
}
