#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=14
#SBATCH --mem=56G
#SBATCH --time=36:00:00

module load R/4.1.0

Rscript run_SBP.R $TRAIT_FILE $LD_DIR $WEIGHT_FILE $CONFIG_FILE $OUTNAME_E $OUTNAME $OUTDIR
