cd /gpfs/data/pierce-lab/uk-biobank-genotypes/; /gpfs/data/xhe-lab/uk-biobank/tools/plink2 --bgen ukb_imp_chr16_v3.bgen --chr 16 --sample /gpfs/data/xhe-lab/uk-biobank/data/genotypes/ukb27386_imp_v3_s487324.sample --keep /gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_samples_0.1.txt --extract /gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_variants.txt --threads 19 --memory 30000 --make-pgen --out /gpfs/data/xhe-lab/ukb_LDR/genotype_data/ukb_chr16
