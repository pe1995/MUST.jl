#!/bin/sh
#SBATCH --job-name='TSOGRID'
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=30
#SBATCH --mem=200G
#SBATCH --time=05:00:00

julia $1