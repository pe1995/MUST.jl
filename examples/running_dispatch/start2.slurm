#!/bin/bash
#SBATCH --job-name=MUST2     
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=eitner@mpia.de        
#SBATCH --ntasks=12              # Number of parallel workers
#SBATCH --ntasks-per-node=4      # Maximum tasks per node
#SBATCH --ntasks-per-core=1      # No hyperthreading               
#SBATCH --cpus-per-task=10       # Number of cores per task
#SBATCH --nodes=3                # distributed over how many nodes
#SBATCH --mem=128GB              # Memory per node, each step gets constant (!) mem / tasks_per_node  
#SBATCH --time=12:00:00             

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=516M

julia create_MUSTGrid.jl -i [51,100] -p summary_solar2