CC         =  mpicc
CCFLAGS    =  -O3 -std=c99
LIBS       =  -lmpi -lm

mm_cannon: mm_cannon.c
	$(CC) $(CCFLAGS) -o mm_cannon mm_cannon.c $(LIBS)

mm_fox: mm_fox.c
	$(CC) $(CCFLAGS) -o mmfox mm_fox.c $(LIBS)


clean:
	rm -f mm_cannon
	rm -f mm_fox
	rm -f output.txt
