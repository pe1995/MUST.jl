#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --partition=medium
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --time=01:00:00                # Time limit hrs:min:sec

julia convert.jl data/stellar_i80/