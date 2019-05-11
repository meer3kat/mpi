#!/bin/bash -l
#SBATCH -A g2019005
#SBATCH -t 10:00
#SBATCH -p core -n 4

module load gcc openmpi 
mpirun -np 4 ./qsort /proj/g2019005/nobackup/qsort_indata/backwards93.txt outputa.txt 0
mpirun -np 4 ./qsort /proj/g2019005/nobackup/qsort_indata/backwards93.txt outputa.txt 1

