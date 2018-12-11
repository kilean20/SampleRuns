#!/bin/bash -l

#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --job-name=my_job
#SBATCH --output my_job.o%j
#SBATCH --license=SCRATCH
#SBATCH -C haswell

#Cori Haswell has 32 physical cores per node
srun -n 1 time xmain

