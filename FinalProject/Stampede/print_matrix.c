

#define A( i,j ) a[ (j)*lda + (i) ]

void print_matrix( int m, int n, double *a, int lda )
{
  int i,j;

  for ( i=0; i<m; i++ ){
    for ( j=0; j<n; j++ )
      printf("%5.1le ", A( i,j ) );
    printf("\n");
  }
}
