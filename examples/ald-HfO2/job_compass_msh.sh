#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=phys.plexicom.q
#SBATCH --error=slurm-%j.err
#SBATCH --output=slurm-%j.out
#SBATCH --time=120:00:00

export exe="./spk_mpi"

source /home/phys/mshirazi/intel/compilers_and_libraries_2020.0.166/linux/bin/compilervars.sh  intel64

mpirun $exe < in.ald
