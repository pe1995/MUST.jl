#!/bin/bash
#SBATCH --job-name=MUST           
#SBATCH --ntasks=12              # Number of parallel workers
#SBATCH --ntasks-per-node=4      # Maximum tasks per node
#SBATCH --ntasks-per-core=1      # No hyperthreading               
#SBATCH --cpus-per-task=10       # Number of cores per task
#SBATCH --nodes=3                # distributed over how many nodes
#SBATCH --mem=64GB               # Memory per node, each step gets constant (!) mem / tasks_per_node  
#SBATCH --time=12:00:00             

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=516M

julia create_MUSTGrid.jl -i [1,100] -p summary_j1