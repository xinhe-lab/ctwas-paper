
#!/bin/bash

#PBS -N bgen2plink
#PBS -S /bin/bash
#PBS -l mem=32gb
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=20

cd /gpfs/data/xhe-lab/ukb_LDR/genotype_data_0.1; bash prep_ukbiobank_chr14_genotype.sh
