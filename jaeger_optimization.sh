#! /bin/bash -l
#SBATCH -D /home/jkta/projects/jaeger
#SBATCH -o /home/jkta/projects/jaeger/logs/out-%j.txt
#SBATCH -e /home/jkta/projects/jaeger/logs/error-%j.txt
#SBATCH -J jaeger_opt

module load R
R

R CMD BATCH --no-save --no-restore "--args ${SLURM_ARRAY_TASK_ID}" Optimization_Cluster.R logs/Rout-${SLURM_ARRAY_TASK_ID}.txt
echo "Done"