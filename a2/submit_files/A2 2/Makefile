CC         =  mpicc
CCFLAGS    =  -O3 -std=c99
LIBS       =  -lmpi -lm

all:
	$(CC) $(CCFLAGS) -o quicksort quicksort.c $(LIBS)

clean:
	rm -f quicksort