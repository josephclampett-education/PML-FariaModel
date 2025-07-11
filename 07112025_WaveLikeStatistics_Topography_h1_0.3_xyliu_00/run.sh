#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=xyliu@unc.edu
#SBATCH -t 07-00:00:00

module load matlab/2023b

matlab -nodesktop -nosplash -r main_bath -logfile output_bath.out
matlab -nodesktop -nosplash -r main -logfile output.out