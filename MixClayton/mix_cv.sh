#!/bin/bash
#SBATCH --account=def-junpark
#SBATCH --time=00:40:00
#SBATCH --job-name="Real"
#SBATCH --array=105-208
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=ruyi.pan@mail.utoronto.ca
#SBATCH --mail-type=END
#SBATCH -o testoutputReal.txt
#SBATCH -e testerrorReal.txt

module load StdEnv gcc/9.3.0 r/4.3.1
cd /home/panruyi/scratch/MixCopula/Codes


r_output_path=//home/panruyi/scratch/MixCopula/Results

export R_LIBS=/home/panruyi/scratch/R_libs/


  srun Rscript PM_O3.R $SLURM_JOB_NAME $SLURM_ARRAY_TASK_ID $r_output_path 

