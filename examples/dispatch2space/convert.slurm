#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=10
#SBATCH --mem=96gb                     # Job memory request
#SBATCH --time=16:00:00                # Time limit hrs:min:sec

julia convert.jl $1
