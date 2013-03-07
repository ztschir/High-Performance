/* Create macros so that the matrices are stored in column-major order */

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

  for ( j=0; j<n; j+=4 ){
    double *c0p = &C(0,j), *c1p = &C(0,j+1), *c2p = &C(0,j+2), *c3p = &C(0,j+3);

    for ( i=0; i<m; i++ ){
      register double c0=0.0, c1=0.0, c2=0.0, c3=0.0;
      double *ap = &AT(0,i), *b0p=&B(0,j  ), *b1p=&B(0,j+1), 
                            *b2p=&B(0,j+2), *b3p=&B(0,j+3);

      for ( p=0; p<k; p+=4 ){
        register double aip;

        aip = *ap++;

        c0 = c0 + aip * b0p[0];
        c1 = c1 + aip * b1p[0];
        c2 = c2 + aip * b2p[0];
        c3 = c3 + aip * b3p[0];

        aip = *ap++;

        c0 = c0 + aip * b0p[1];
        c1 = c1 + aip * b1p[1];
        c2 = c2 + aip * b2p[1];
        c3 = c3 + aip * b3p[1];

        aip = *ap++;

        c0 = c0 + aip * b0p[2];
        c1 = c1 + aip * b1p[2];
        c2 = c2 + aip * b2p[2];
        c3 = c3 + aip * b3p[2];

        aip = *ap++;

        c0 = c0 + aip * b0p[3];
        c1 = c1 + aip * b1p[3];
        c2 = c2 + aip * b2p[3];
        c3 = c3 + aip * b3p[3];

        b0p+=4;
        b1p+=4;
        b2p+=4;
        b3p+=4;
      }
      
      *c0p++ += c0;
      *c1p++ += c1;
      *c2p++ += c2;
      *c3p++ += c3;
    }
  }
}
