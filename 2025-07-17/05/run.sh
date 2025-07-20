#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mem=32g
#SBATCH --array=1-10
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=jclampett@unc.edu
#SBATCH -t 07-00:00:00

eval "scontrol update jobid=${SLURM_JOB_ID} name=ID[${SLURM_ARRAY_TASK_ID}]"
echo "ID[${SLURM_ARRAY_TASK_ID}]"
module load matlab/2024b
matlab -nodesktop -nosplash -singleCompThread -r "batch_main($SLURM_ARRAY_TASK_ID)" -logfile output_main.out
# matlab -nodesktop -nosplash -singleCompThread -r postprocess -logfile output_postprocess.out