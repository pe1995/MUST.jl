#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=10
##SBATCH --mem=200gb                     # Job memory request
#SBATCH --mem=92GB                     # Job memory request
#SBATCH --time=5:00:00                # Time limit hrs:min:sec

julia convert.jl $1
