#!/bin/bash
#SBATCH --job-name=MUST           
#SBATCH --ntasks=12                    
#SBATCH --cpus-per-task=12
#SBATCH --nodes=2
#SBATCH --mem-per-cpu=1333mb    
#SBATCH --ntasks-per-core=1      # No hyperthreading            
#SBATCH --time=01:30:00             

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_STACKSIZE=516M

#srun -N 1 -n 1 -c $SLURM_CPUS_PER_TASK --mem=16gb --exclusive -D ../../../dispatch2/experiments/stellar_atmospheres/ ./dispatch.x grid1_phase1.nml > grid1_phase1.log &
#srun -N 1 -n 1 -c $SLURM_CPUS_PER_TASK --mem=16gb --exclusive -D ../../../dispatch2/experiments/stellar_atmospheres/ ./dispatch.x grid2_phase1.nml > grid2_phase1.log &
#srun -N 1 -n 1 -c $SLURM_CPUS_PER_TASK --mem=16gb --exclusive -D ../../../dispatch2/experiments/stellar_atmospheres/ ./dispatch.x grid3_phase1.nml > grid3_phase1.log &
#srun -N 1 -n 1 -c $SLURM_CPUS_PER_TASK --mem=16gb --exclusive -D ../../../dispatch2/experiments/stellar_atmospheres/ ./dispatch.x grid4_phase1.nml > grid4_phase1.log &

#wait

julia test_dispatch.jl