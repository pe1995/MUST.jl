#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=2
#SBATCH --ntasks=42
##SBATCH --partition=medium
##SBATCH --nodes=1
##SBATCH --ntasks=20
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=00:40:00                # Time limit hrs:min:sec

julia convert.jl $1
