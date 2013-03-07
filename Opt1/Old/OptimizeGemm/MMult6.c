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
    double *cp = &C( 0, j );

    for ( i=0; i<m; i+=4 ){
      register double c0=0.0, c1=0.0, c2=0.0, c3=0.0;
      double *ap = &A(i,0), *bp=&B(0,j);

      for ( p=0; p<k; p++ ){
        register double bpj = *bp;
        register double r0, r1, r2;

        r0 = ap[ 0 ];
        r1 = ap[ 1 ];
        r2 = ap[ 2 ];
        c0 = c0 +  r0 * bpj;

        r0 = ap[ 3 ];
        c1 = c1 +  r1 * bpj;
        c2 = c2 +  r2 * bpj;
        c3 = c3 +  r0 * bpj;

        bp++;
        ap+=lda;
      }
      
      cp[ 0 ] += c0;
      cp[ 1 ] += c1;
      cp[ 2 ] += c2;
      cp[ 3 ] += c3;
      
      cp += 4;
    }
  }
}
