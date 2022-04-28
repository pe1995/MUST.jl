#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --partition=debug
#SBATCH --time=7-00:00:00

julia convert.jl data/stellar_i120/