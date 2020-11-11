#!/bin/bash 

#SBATCH --job-name=test
#SBATCH --partition=v100
#SBATCH --output=test24.out
#SBATCH --error=test24.err
#SBATCH --array=0-0
#SBATCH --nodes=1
#SBATCH -n 4
#SBATCH --exclude=""
#SBATCH --exclusive
#SBATCH --gres=gpu:4


export OMP_NUM_THREADS=8

mpirun -np 4 ./chroma -i overlap_24.xml -geom 1 1 1 4 > output-test.24
