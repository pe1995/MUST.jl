#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32gb                     # Job memory request
#SBATCH --time=05:00:00                # Time limit hrs:min:sec

julia convert_TSO.jl $1
