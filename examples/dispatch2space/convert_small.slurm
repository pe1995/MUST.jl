#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb                     # Job memory request
#SBATCH --time=01:00:00                # Time limit hrs:min:sec

julia convert.jl $1
