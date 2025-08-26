#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jclampett@unc.edu
#SBATCH -t 07-00:00:00

module load matlab/2024b
matlab -nodesktop -nosplash -singleCompThread -r global_postprocess -logfile output_global_postprocess.out
