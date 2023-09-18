#!/bin/sh
#SBATCH --job-name='TSOGRID'
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=100G
#SBATCH --partition=debug
#SBATCH --time=5:00:00                # Time limit hrs:min:sec

julia "$@"