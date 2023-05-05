#$ -S /bin/bash

#$ -q serial
#$ -N misspecification
#$ -t 1:8

source /etc/profile

module add R

echo Job task $SGE_TASK_ID 

R CMD BATCH --no-save misspecsims.R 
