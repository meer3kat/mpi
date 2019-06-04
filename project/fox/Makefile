CC         =  mpicc
CCFLAGS    =  -Wall -O3 -std=c99
LIBS       =  -lmpi -lm

all:
	$(CC) $(CCFLAGS) -o mm_fox mm_fox.c $(LIBS)

clean:
	rm -f mm_fox
