#!/bin/sh
#SBATCH --job-name='TSOGRID'
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
##SBATCH --mem=90G
#SBATCH --partition=debug

julia "$@"