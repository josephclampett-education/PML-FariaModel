#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=64g
#SBATCH -t 07-00:00:00

module load matlab/2023b
matlab -nodesktop -nosplash -singleCompThread -r main_2025_07_03_WavelikeAttempt0 -logfile output.out


