#!/bin/bash
#SBATCH --job-name=MUST                # Job name
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1
#SBATCH --mem=32gb                     # Job memory request
#SBATCH --time=1:00:00                # Time limit hrs:min:sec

julia $1 $2 $3