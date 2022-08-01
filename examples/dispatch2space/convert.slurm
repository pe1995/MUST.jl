#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=3
#SBATCH --ntasks=30
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=2
#SBATCH --mem=96gb                     # Job memory request
#SBATCH --time=01:00:00                # Time limit hrs:min:sec

julia convert.jl $1
