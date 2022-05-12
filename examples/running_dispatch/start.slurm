#!/bin/bash
#SBATCH --job-name=MUST           
#SBATCH --ntasks=30              # Number of parallel workers
#SBATCH --ntasks-per-node=6      # Maximum tasks per node
#SBATCH --ntasks-per-core=1      # No hyperthreading               
#SBATCH --cpus-per-task=12       # Number of cores per task
#SBATCH --nodes=5                # distributed over how many nodes
#SBATCH --mem=96GB               # Memory per node, each step gets constant (!) mem / tasks_per_node  
#SBATCH --time=24:00:00             

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=516M

julia create_MUSTGrid.jl 200