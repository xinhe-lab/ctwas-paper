Download the list of samples and variants used by the Neale Lab: 'samples.both_sexes.txt' and 'variants.txt' as documented in their manifest: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679
Run 'neale_process_samples_variants.R' in the 'neale_lab' folder
Run 'liftover_ukbiobank.sh' in the 'genotype_data_0.1' folder to convert variant positions from b37 to b38
Run the 'prep_ukbiobank_chr*_genotype.sh' files in the 'genotype_data_0.1' folder to subset to the specified variants and samples. The 'qsub_' files specify the resources needed.
Run 'convert_geno_to_LDR_chr.R' in the 'matrices_0.1' folder for each chromosome. The 'qsub_' files specify the resources needed.
