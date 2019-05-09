#!/bin/bash -l
#SBATCH -A g2018006
#SBATCH  g2018006_1
#SBATCH -t 5:00
module load gcc openmpi

make A2


mpirun ./A2 40000
mpirun ./A2 200000
mpirun ./A2 1000000
mpirun ./A2 10000000
mpirun ./A2 100000000

