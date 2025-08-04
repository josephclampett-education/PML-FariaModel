#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks=36
#SBATCH --cpus-per-task=1
#SBATCH --mem=100g
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jclampett@unc.edu
#SBATCH -t 07-00:00:00

module load matlab/2024b
matlab -nodesktop -nosplash -r main -logfile output_main.out