#!/bin/bash
#SBATCH --job-name=MUST1   
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=eitner@mpia.de      
#SBATCH --ntasks=5               # Number of parallel workers
#SBATCH --ntasks-per-node=1      # Maximum tasks per node
#SBATCH --ntasks-per-core=1      # No hyperthreading               
#SBATCH --cpus-per-task=40       # Number of cores per task
#SBATCH --nodes=5                # distributed over how many nodes
#SBATCH --mem=96GB               # Memory per node, each step gets constant (!) mem / tasks_per_node  
#SBATCH --time=24:00:00             

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=516M

julia create_MUSTGrid_onePhase.jl -i [1,30] -p summary_highRho_1