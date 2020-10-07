#!/bin/bash 

#SBATCH --job-name=test
#SBATCH --partition=fat3t
#SBATCH --output=test.sh.out
#SBATCH --error=test.sh.err
#SBATCH --array=0-0
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --exclude=""
##SBATCH --ntasks-per-node=16
##SBATCH --ntasks-per-socket=16
#SBATCH --exclusive
#SBATCH --gres=gpu:1


export OMP_NUM_THREADS=16

mpirun -np 1 ./chroma -i overlap.xml -geom 1 1 1 1 >output.08
#QUDA_ENABLE_TUNING=0 QUDA_ENABLE_P2P=0 mpirun -np 1 ./bind_gpu.sh_quda  ./chroma_cuda -i ini.1005.xml.test -geom 1 1 1 1 >out.test 2>&1
