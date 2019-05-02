#!/bin/bash -l

#SBATCH -A g2019005
#SBATCH --reservation g2019005_1
#SBATCH -p core -n 2
#SBATCH -t 5:00

module load gcc openmpi
mpirun -np 2 alltoall
