void CopyMatrixGlobalToLocal( int, int, 
			      double *, int, 
			      double *, int, 
			      MPI_Comm,
			      MPI_Comm );

void ParallelRank1( int, int, int,
		    int, 
		    double *, int,
		    double *, int,
		    double *, int,
		    MPI_Comm, MPI_Comm );

void ParallelMMult( int, int, int,
		    double *, int,
		    double *, int,
		    double *, int,
		    MPI_Comm, MPI_Comm );

double compare_matrices( int, int, double *, int, double *, int );

void random_matrix( int, int, double *, int );
void print_matrix( int, int, double *, int );

void dger_( int *, int *, double *, double *, int *, double *, int *, double *, int * );

void dgemm_( char *, char *, int *, int *, int *,
	     double *, double *, int *, double *, int *, 
	     double *, double *, int * );



