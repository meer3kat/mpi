#!/bin/bash -l
#SBATCH -A g2019005
#SBATCH -t 1000:00
#SBATCH -p core -n 9


module load gcc openmpi 
mpirun -np 9 ./mm_fox /proj/g2019005/nobackup/matmul_indata/input3600.txt output.txt
mpirun -np 9 ./mm_fox /proj/g2019005/nobackup/matmul_indata/input7488.txt output.txt
mpirun -np 9 ./mm_fox /proj/g2019005/nobackup/matmul_indata/input9072.txt output.txt
