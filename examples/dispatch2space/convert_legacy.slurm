#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=2
#SBATCH --ntasks=62
##SBATCH --partition=medium
##SBATCH --nodes=1
##SBATCH --ntasks=20
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=03:00:00                # Time limit hrs:min:sec

julia convert_legacy.jl $1
