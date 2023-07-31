#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8

julia convert.jl $1
