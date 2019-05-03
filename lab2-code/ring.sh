#!/bin/bash -l

#SBATCH -A g2019005
#SBATCH -p core -n 10
#SBATCH -t 5:00

module load gcc openmpi
mpirun -np 10 ./ring
