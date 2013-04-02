CC = mpicc

NPROCS = 5

RUN = mpirun -np $(NPROCS)

SIMPLE = hello.o

hello.x: $(SIMPLE)
	 $(CC) -o hello.x $(SIMPLE)

run:	hello.x
	$(RUN) hello.x

clean:
	rm -f *.o *.x *~ core
