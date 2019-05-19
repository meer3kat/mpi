#!/bin/bash -l
#SBATCH -A g2019005
#SBATCH -t 10:00
#SBATCH -p core -n 8

module load gcc openmpi 
mpirun -np 8 ./quicksort input10.txt outputb.txt 1
