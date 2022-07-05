#!/bin/sh
#SBATCH --partition=batch
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --job-name="water2000"
#SBATCH --output=eam.out

module load gcc/6.2.0


echo "$NPROC"
#mpiexec -n $NPROC lmp_mpi <../in.eam> log.lammps 
time poreblazer input.dat >poreblazer.log
