CC         =  mpicc
CCFLAGS    =  -Wall -O2 -std=c99
LIBS       =  -lmpi -lm

all:
	$(CC) $(CCFLAGS) -o qsort qsort.c $(LIBS)

qsort_test: qsort_test.c
	$(CC) $(CCFLAGS) -o qsort_test qsort_test.c $(LIBS)

clean:
	rm -f qsort