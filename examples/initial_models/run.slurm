#!/bin/sh
#SBATCH --job-name='MUST'
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
##SBATCH --mem=200G
#SBATCH --mem=90G
#SBATCH --time=05:00:00

julia $@