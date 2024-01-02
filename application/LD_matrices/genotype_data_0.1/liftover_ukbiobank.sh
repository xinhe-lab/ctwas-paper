#!/bin/bash

module load gcc/6.2.0 LiftOver

for chr in $(seq 1 1 2)
do
	echo $i
	/gpfs/data/xhe-lab/uk-biobank/tools/plink2 --make-bed -pfile ukb_chr$chr --out ukb_chr$chr
	awk '{ print "chr"$1, $4-1, $4, $4, $2 }' ukb_chr$chr.bim > ukb_chr$chr.ucsc_bed
	liftOver ukb_chr$chr.ucsc_bed hg19ToHg38.over.chain.gz ukb_chr$chr.ucsc_bed_b38 ukb_chr$chr.unmapped_temp
	grep -v "#" ukb_chr$chr.unmapped_temp > ukb_chr$chr.unmapped; rm ukb_chr$chr.unmapped_temp
	for i in $(seq 1 1 22); do if [ $i == $chr ]; then continue; fi; cat ukb_chr$chr.ucsc_bed_b38 |grep chr$i$'\t'  >> ukb_chr$chr.unmapped; done
	awk '{print $5,$3}' ukb_chr$chr.ucsc_bed_b38 > ukb_chr$chr.b38_positions
	/gpfs/data/xhe-lab/uk-biobank/tools/plink2 --bed ukb_chr$chr.bed --bim ukb_chr$chr.bim --fam ukb_chr$chr.fam --exclude ukb_chr$chr.unmapped --update-map ukb_chr$chr.b38_positions --make-pgen --sort-vars  --out ukb_chr$chr.b38
	/gpfs/data/xhe-lab/uk-biobank/tools/plink2 --make-bed -pfile ukb_chr$chr.b38 --out ukb_chr$chr.b38
	rm ukb_chr$chr.b38.pgen ukb_chr$chr.b38.psam ukb_chr$chr.b38.pvar ukb_chr$chr.b38.log ukb_chr$chr.bed ukb_chr$chr.bim ukb_chr$chr.fam ukb_chr$chr.log ukb_chr$chr.ucsc_bed ukb_chr$chr.ucsc_bed_b38 ukb_chr$chr.unmapped ukb_chr$chr.b38_positions
done

