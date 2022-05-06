#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=2
#SBATCH --ntasks=42
#SBATCH --partition=medium
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=02:00:00                # Time limit hrs:min:sec

julia convert.jl $1
