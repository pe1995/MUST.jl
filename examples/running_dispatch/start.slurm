#!/bin/bash
#SBATCH --job-name=RUN_DISPATCH           # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=eitner@mpia.de    # Where to send mail
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --time=24:00:00                # Time limit hrs:min:sec


#module load intel/21.3.0
#module load impi/2021.3

module load julia/1.6

julia running_namelists.jl 150
