#!/bin/bash

N_JOBS=$(cat trait_weight_pairs_nolnc.csv|wc -l)

for JOB_INDEX in `seq 1 $N_JOBS`
do
	TRAIT=$(cut -d"," -f"2" trait_weight_pairs_nolnc.csv | sed -n "$JOB_INDEX"p)
	TRAIT_FILE="/project2/compbio/gwas_summary_statistics/ukbb_neale_v3/"$TRAIT".vcf"
	WEIGHT=$(cut -d"," -f"3" trait_weight_pairs_nolnc.csv | sed -n "$JOB_INDEX"p)
	LD_DIR="/project2/mstephens/wcrouse/UKB_LDR_0.1"
	WEIGHT_FILE="/project2/mstephens/wcrouse/predictdb_nolnc/mashr_"$WEIGHT"_nolnc.db"
	CONFIG_FILE="/project2/mstephens/wcrouse/UKB_analysis_allweights_corrected/ctwas_config.R"
	OUTNAME_E=$TRAIT"_"$WEIGHT"_expr"
	OUTNAME=$TRAIT"_"$WEIGHT"_ctwas"
	OUTDIR="/project2/mstephens/wcrouse/UKB_analysis_allweights_corrected/"$TRAIT"/"$WEIGHT"_nolnc"
	job_name=$TRAIT"_"$WEIGHT
	job_out="out/"$job_name".out"
	job_err="out/"$job_name".err"
  
	FINAL_FILE=$OUTDIR"/"$OUTNAME".susieIrss.txt"
	if [ ! -f $FINAL_FILE ]
	then 
		echo sbatch --export=TRAIT_FILE=$TRAIT_FILE,LD_DIR=$LD_DIR,WEIGHT_FILE=$WEIGHT_FILE,CONFIG_FILE=$CONFIG_FILE,OUTNAME_E=$OUTNAME_E,OUTNAME=$OUTNAME,OUTDIR=$OUTDIR --job-name=$job_name --output=$job_out --error=$job_err run_IBD.sbatch
	fi
done

