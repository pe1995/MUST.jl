#!/bin/sh
#SBATCH --job-name='TSOGRID'
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
##SBATCH --mem=200G
#SBATCH --mem=92G
#SBATCH --time=05:00:00

julia $1