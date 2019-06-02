CC         =  mpicc
CCFLAGS    =  -Wall -O3 -std=c99
LIBS       =  -lmpi -lm

all:
	$(CC) $(CCFLAGS) -o cannon cannon.c $(LIBS)

clean:
	rm -f cannon