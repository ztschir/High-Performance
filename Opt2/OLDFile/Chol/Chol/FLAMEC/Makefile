# BLAS_LIB   := $(HOME)/flame/GotoBLAS2/libgoto2.a
BLAS_LIB   := $(HOME)/flame/OpenBLAS/libopenblas_penryn-r0.1.1.a
LAPACK_LIB := # $(HOME)/flame/lapack-3.3.0/lapack_LINUX.a
FLAME_LIB  := $(HOME)/flame/lib/libflame.a
FLAME_INC  := $(HOME)/flame/include

# indicate where the object files are to be created
CC         := gcc
FC         := gfortran
LINKER     := $(FC)
CFLAGS     := -O2 -Wall -I$(FLAME_INC) -m64 -msse3
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.
MAX_GFLOPS := 40.0#     GFLOPS that defines the top of the graph.  If smaller 
                    #  than observed, the top of the graph is the best observed
                    #  performance.
NB_ALG     := 128#     algorithmic block size to be used.
NFIRST     := 200#     smallest size to be timed.
NLAST      := 3000#    largest size to be timed.
NINC       := 200#     increment between sizes.
NLAST_UNB  := 2000#    largest problem size to be timed for unblocked variants

LDFLAGS    := -lm -lpthread -m64

# indicate where the FLAME include files reside

TEST_OBJS  := driver.o REF_Chol.o \
	Chol_unb_var3.o \
	Chol_blk_var3.o \
	Chol_unb_var1.o \
	Chol_unb_var2.o \
	Chol_blk_var1.o \
	Chol_blk_var2.o \

LAPACK_OBJS = # dpotrfx.o

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

driver.x: $(TEST_OBJS)  $(LAPACK_OBJS)
	$(LINKER) $(TEST_OBJS) $(LAPACK_OBJS) $(LDFLAGS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) -o driver.x

test:   driver.x
	echo "$(NREPEATS) $(MAX_GFLOPS) $(NB_ALG) $(NFIRST) $(NLAST) $(NINC) $(NLAST_UNB)" | ./driver.x > output.m

clean:
	rm -f *.o *~ core *.x output.m


