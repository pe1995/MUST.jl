#!/bin/bash
#SBATCH --ntasks=3                    
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --mem=92gb                     # Job memory request
#SBATCH --time=24:00:00                # Time limit hrs:min:sec

srun --exclusive --mem=30G -n1 julia --project monitor.jl grid_t50g40m00_h10_e1
srun --exclusive --mem=30G -n1 julia --project monitor.jl grid_t50g40m00_h13_e1
srun --exclusive --mem=30G -n1 julia --project monitor.jl grid_t50g40m00_h18_e1