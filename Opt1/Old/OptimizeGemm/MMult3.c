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
      double *ap = &A(i,0), *bp=&B(0,j);

      for ( p=0; p<k; p++ ){
        c0 = c0 + ap[0] * *bp;
        c1 = c1 + ap[1] * *bp;
        c2 = c2 + ap[2] * *bp;
        c3 = c3 + ap[3] * *bp;

        bp++;
        ap += lda;
      }
      
      C( i,j )   += c0;
      C( i+1,j ) += c1;
      C( i+2,j ) += c2;
      C( i+3,j ) += c3;
    }
  }
}


  
