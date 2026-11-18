#!/bin/bash -l
#SBATCH --partition=demo
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=%x-%j.out

set -e
module load openmpi
N=$1
BIN=build-$N/src

echo "=== sequential cached, 1 process ==="
mpirun -np 1 $BIN/algebra_determinants
echo "=== mpi efficient, 8 ranks ==="
mpirun $BIN/mpi_efficient
echo "=== mpi inefficient, 8 ranks ==="
mpirun $BIN/mpi_inefficient