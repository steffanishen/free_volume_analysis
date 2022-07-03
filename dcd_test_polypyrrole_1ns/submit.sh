#!/bin/sh
#SBATCH --partition=batch
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name="water2000"
#SBATCH --output=eam.out

module load gcc/6.2.0


echo "$NPROC"
#mpiexec -n $NPROC lmp_mpi <../in.eam> log.lammps 
time /projects/project_interface/mshen/post_processing/Fortran/poreblazer_v3.0.2_full_neighbors_no_z_PBC_n_only_read_dcd/poreblazer.exe input.dat >poreblazer.log
