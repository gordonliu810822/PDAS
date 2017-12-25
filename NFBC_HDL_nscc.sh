#!/bin/bash
#PBS -N NFBC_HDL
#PBS -l select=1:ncpus=12:mem=48GB
#PBS -l walltime=02:00:00
#PBS -l place=pack
#PBS -j oe
#PBS -P 11000368
#PBS -e /home/projects/11000369/Work/NFBC/NFBC_HDL.error
#PBS -o /home/projects/11000369/Work/NFBC/NFBC_HDL.log

module load gcc/5.1.0
module load R/3.3.1
R CMD BATCH /home/projects/11000369/Work/NFBC/qsub_sh/NFBC_HDL_nscc.R
