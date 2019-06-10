There are two version of matrix matrix multipication implemented: Fox's and Cannon's algorithm.
To compile fox's algorithm use: make mm_fox
To compile Cannon's algorithm use: make mm_cannon

the code takes 2 arguments, an input file and an output file. all input file are located in the course project folder on uppmax. 

to run, for example Cannon's algorithm use:
mpirun -np 16 ./mm_cannon /proj/g2019005/nobackup/matmul_indata/input3600.txt output.txt
