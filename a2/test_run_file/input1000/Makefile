CC         =  mpicc
CCFLAGS    =  -Wall -O3 -std=c99
LIBS       =  -lmpi -lm

all:
	$(CC) $(CCFLAGS) -o qsort qsort.c $(LIBS)

clean:
	rm -f qsort