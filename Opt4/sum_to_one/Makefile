CC = mpicc

NPROCS = 5

RUN = mpiexec  -machinefile hostfile -n $(NPROCS)

compile: sum_to_one.c my_sum_to_one_new.c
	$(CC) -o sum_to_one.x sum_to_one.c my_sum_to_one_new.c -lm

run:	compile
	$(RUN) ./sum_to_one.x

clean:
	rm *.o *.x *~


