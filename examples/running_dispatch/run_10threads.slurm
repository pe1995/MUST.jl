#!/bin/bash
#SBATCH --job-name=MUST           # Job name
##SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=eitner@mpia.de    # Where to send mail
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
##SBATCH --partition=medium
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --time=01:00:00                # Time limit hrs:min:sec
##SBATCH -o MUST_12threads.log        # Standard output and error log

#module load intel/21.3.0
#module load impi/2021.3

export OMP_STACKSIZE=516M
export OMP_NUM_THREADS=10

srun ./dispatch.x $1
