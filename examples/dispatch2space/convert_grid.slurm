#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=5
#SBATCH --ntasks=50
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=4
#SBATCH --mem=96gb                     # Job memory request
#SBATCH --time=05:00:00                # Time limit hrs:min:sec

julia convert_grid.jl $1
