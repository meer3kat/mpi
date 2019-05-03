#!/bin/bash -l

#SBATCH -A g2019005
#SBATCH -p core -n 4
#SBATCH -t 5:00

module load gcc openmpi
mpirun -np 4 alltoall
