#!/bin/bash -l
#SBATCH -A g2019005
#SBATCH -t 100:00
#SBATCH -p node -n 64


module load gcc openmpi 
mpirun -np 64 ./mm_fox /proj/g2019005/nobackup/matmul_indata/input3600.txt output.txt
mpirun -np 64 ./mm_fox /proj/g2019005/nobackup/matmul_indata/input7488.txt output.txt
mpirun -np 64 ./mm_fox /proj/g2019005/nobackup/matmul_indata/input9072.txt output.txt


