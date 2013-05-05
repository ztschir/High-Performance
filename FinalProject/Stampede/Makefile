CC = mpicc

OBJ = driver.o ParallelMMult_2.o ParallelRankK_4.o CopyMatrixGlobalToLocal.o random_matrix.o compare_matrices.o print_matrix.o

BLAS = -mkl 

mmmult.x: $(OBJ)
	 $(CC) -o mmmult.x $(OBJ) $(BLAS)

run:	mmmult.x
	sbatch batch.tacc

clean:
	rm -f *.o *.x *~ core *.out *.err

