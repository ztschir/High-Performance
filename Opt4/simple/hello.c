#include <stdio.h>
#include "mpi.h"

void main(int argc, char **argv )
{
  MPI_Init( &argc, &argv );

  printf("hello world\n" );

  MPI_Finalize();

  exit( 0 );
}
