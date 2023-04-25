#!/bin/sh
#SBATCH --job-name='TSOGRID'
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
##SBATCH --mem=90G
#SBATCH --partition=debug

julia $1