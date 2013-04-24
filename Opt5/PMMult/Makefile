CC = mpicc

NPROCS = 6

RUN = mpiexec -np $(NPROCS)

OBJ = driver.o ParallelMMult_2.o ParallelRankK_2.o CopyMatrixGlobalToLocal.o random_matrix.o compare_matrices.o print_matrix.o

BLAS = ~rvdg/flame/OpenBLAS/libopenblas_core2-r0.1.1.a

mmmult.x: $(OBJ)
	 $(CC) -o mmmult.x $(OBJ) $(BLAS)

run:	mmmult.x
	echo "100 3 2" | $(RUN) mmmult.x

clean:
	rm -f *.o *.x *~ core
