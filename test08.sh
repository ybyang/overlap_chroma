#!/bin/bash 

#SBATCH --job-name=test08
#SBATCH --partition=fat3t
#SBATCH --output=test08.out
#SBATCH --error=test08.err
#SBATCH --array=0-0
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --exclude=""
#SBATCH --exclusive
#SBATCH --gres=gpu:1


export OMP_NUM_THREADS=8

mpirun -np 1 ./chroma -i overlap_08.xml -geom 1 1 1 1 >output.08
