#!/bin/bash
#SBATCH --job-name=MUST           
#SBATCH --ntasks=5   
#SBATCH --ntasks-per-node=6                
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=96GB
##SBATCH --mem-per-cpu=1333mb    
#SBATCH --ntasks-per-core=1      # No hyperthreading            
#SBATCH --time=01:30:00             

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=516M

julia create_MUSTGrid.jl 2