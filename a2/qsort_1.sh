#!/bin/bash -l
#SBATCH -A g2019005
#SBATCH -t 25:00
#SBATCH -p core -n 1

module load gcc openmpi 
mpirun -np 1 ./qsort /proj/g2019005/nobackup/qsort_indata/backwards125000000.txt outputa.txt 0
mpirun -np 1 ./qsort /proj/g2019005/nobackup/qsort_indata/backwards125000000.txt outputa.txt 1
mpirun -np 1 ./qsort /proj/g2019005/nobackup/qsort_indata/backwards125000000.txt outputa.txt 2

