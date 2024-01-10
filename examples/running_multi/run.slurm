#!/bin/sh
#SBATCH --job-name=convert
#SBATCH --nodes=1
#SBATCH --ntasks=2
##SBATCH --ntasks-per-node=7
#SBATCH --cpus-per-task=35
#SBATCH --mem=200gb                     # Job memory request
#SBATCH --time=12:00:00                # Time limit hrs:min:sec

julia $1