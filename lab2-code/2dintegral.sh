#!/bin/bash -l

#SBATCH -A g2019005
#SBATCH -p core -n 2
#SBATCH -t 3:00

module load gcc openmpi
mpirun -np 2 ./integral2d