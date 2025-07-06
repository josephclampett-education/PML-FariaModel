#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64g
#SBATCH -t 07-00:00:00

module load matlab/2023b
matlab -nodesktop -nosplash -r main_2025_07_05_WavelikeAttempt1 -logfile output.out


