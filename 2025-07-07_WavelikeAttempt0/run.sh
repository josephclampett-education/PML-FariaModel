#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jclampett@unc.edu
#SBATCH -t 07-00:00:00

module load matlab/2023b
matlab -nodesktop -nosplash -r main -logfile output.out


