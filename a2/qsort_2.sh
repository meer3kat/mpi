#!/bin/bash -l
#SBATCH -A g2019005
#SBATCH -t 5:00
#SBATCH -p core -n 8

module load gcc openmpi 
mpirun -np 8 ./qsort /proj/g2019005/nobackup/qsort_indata/input125000000.txt outputa.txt 0
