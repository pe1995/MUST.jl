#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb                     # Job memory request
#SBATCH --time=05:00:00                # Time limit hrs:min:sec

julia convert_fast.jl $1
