#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:30:00
#SBATCH -p eqa-mse_4270_6270-18sp
#SBATCH -A mse_4270_6270
#SBATCH --output=output.txt

# Run program
module load intel
cd $SLURM_SUBMIT_DIR
./runmd
