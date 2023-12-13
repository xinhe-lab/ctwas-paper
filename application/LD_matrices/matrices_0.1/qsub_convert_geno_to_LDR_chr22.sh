#!/bin/bash

#PBS -N ukb_LDR
#PBS -S /bin/bash
#PBS -l mem=16gb
#PBS -l walltime=24:00:00
#PBS -l nodes=1

module load gcc/4.9.4 R/3.6.2 gsl; cd /gpfs/data/xhe-lab/ukb_LDR/matrices_0.1; Rscript convert_geno_to_LDR_chr.R 22 1
