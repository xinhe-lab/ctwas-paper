setwd("/gpfs/data/xhe-lab/ukb_LDR/neale_lab")

library(data.table)
library(dplyr)

####################
#load phenotype data from UK Biobank

# input.file <- file.path("/gpfs/data/xhe-lab/uk-biobank/data/phenotypes",
#                         "12-feb-2019","ukb26140.csv.gz")
# 
# cols      <- c("eid","22007-0.0","22008-0.0")
# col_names <- c("id","plate_name","well")
# 
# out <- system.time({
#   dat <- fread(input.file,sep = ",",header = TRUE,verbose = FALSE,
#                showProgress = T,colClasses = "character");
# })
# 
# class(dat) <- "data.frame"
# cat(sprintf("Data loading step took %d seconds.\n",round(out["elapsed"])))
# cat(sprintf("Table contains %d rows.\n",nrow(dat)))
# 
# cat("Preparing data.\n")
# dat        <- dat[,cols]
# names(dat) <- col_names
# 
# save(dat, file = "/gpfs/data/xhe-lab/ukb_LDR/neale_lab/population_ukb26140.Rd")
load("/gpfs/data/xhe-lab/ukb_LDR/neale_lab/population_ukb26140.Rd")

####################

#load sample file from Neale lab
samples <- fread("data/samples.both_sexes.txt")

#create identifier from plate_name and well
samples$plate_name_well <- paste(samples$plate_name, samples$well, sep="_")

####################

#drop individuals without plate_name or well
dat <- dat[dat$plate_name!="" & dat$well!="",]

#create identifier from plate_name and well
dat$plate_name_well <- paste(dat$plate_name, dat$well, sep="_")

#subset to individuals in Neale lab samples
dat <- dat[dat$plate_name_well %in% samples$plate_name_well,]

#drop individuals not uniquely identified by plate and well
dat <- dat[!(dat$plate_name_well %in% dat$plate_name_well[duplicated(dat$plate_name_well)]),]

cat(sprintf("Table contains %d rows.\n",nrow(dat)))

#output sample file
output.file <- file.path( "/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_samples.txt")
out <- dat[,c("id","id")]
colnames(out) <- c("id","id")
write.table(out, file=output.file, quote=F, row.names=F, col.names=T, sep="\t")

#output a subset of 10% of samples for testing
set.seed(12375548)

subset_index <- sample(1:nrow(out), ceiling(nrow(out)*0.1))

out_subset <- out[subset_index,]

output.file <- file.path( "/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_samples_0.1.txt")
write.table(out_subset, file=output.file, quote=F, row.names=F, col.names=T, sep="\t")

####################

#load variants file from Neale lab
variants <- fread("data/variants.txt")

#subset to MAF > 0.01 in Neale analysis
variants <- variants[variants$minor_AF > 0.01,]

####################

for (i in 1:22){
  print(i)
  
  #load variant lists by chromosome from UK Biobank
  ukb_bim_chr <- fread(paste0("/gpfs/data/xhe-lab/ukb_LDR/neale_lab/data/ukb_imp_chr", i, "_v3.bim"))
  colnames(ukb_bim_chr) <- c("chr", "id", "cm", "pos", "alt", "ref")
  
  #subset to variants included in Neale analysis after MAF filter
  ukb_bim_chr <- ukb_bim_chr[ukb_bim_chr$id %in% variants$rsid,]
  
  #store variants
  if (i==1){
    ukb_bim <- ukb_bim_chr
  } else {
    ukb_bim <- rbind(ukb_bim, ukb_bim_chr)
  }
}

#drop variants not uniquely identified by IDs (e.g. multiallelic)
ukb_bim <- ukb_bim[!(ukb_bim$id %in% ukb_bim$id[duplicated(ukb_bim$id)]),]

#output full bim file
output.file <- file.path( "/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_variants.bim")
write.table(ukb_bim, file=output.file, quote=F, row.names=F, col.names=F, sep="\t")

#output variant list
output.file <- file.path( "/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_variants.txt")
out <- ukb_bim$id
write.table(out, file=output.file, quote=F, row.names=F, col.names=F, sep="\t")

#print number of variants per chromosome
table(ukb_bim$chr)

####################
#liftover to hg38

neale_bim <- fread("/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_variants.bim")
colnames(neale_bim) <- c("chr", "id", "cm", "pos", "alt", "ref")

neale_bim <- as.data.frame(neale_bim)

rownames(neale_bim) <- neale_bim$id

library(liftOver)

ch = import.chain("hg19ToHg38.over.chain")

neale_hg38 <- GRanges(seqnames = paste0("chr",neale_bim$chr), 
                      ranges = IRanges(neale_bim$pos, neale_bim$pos),
                      id = neale_bim$id,
                      cm = neale_bim$cm,
                      alt = neale_bim$alt,
                      ref = neale_bim$ref)

neale_hg38 <- unlist(liftOver(neale_hg38, ch))
neale_hg38 <- as.data.table(neale_hg38)
neale_hg38$chr <- as.numeric(substring(neale_hg38$seqnames, 4))
neale_hg38$pos <- neale_hg38$start
neale_hg38 <- neale_hg38[!is.na(neale_hg38$chr),
                         c("chr", "id", "cm", "pos", "alt", "ref")]

#drop variants that changed chromosomes
neale_hg38 <- neale_hg38[neale_bim[neale_hg38$id,"chr"]==neale_hg38$chr,]

#sort output
neale_hg38 <- neale_hg38[order(neale_hg38$chr, neale_hg38$pos),]

#output full bim file
output.file <- file.path( "/gpfs/data/xhe-lab/ukb_LDR/neale_lab/neale_variants_hg38.bim")
write.table(neale_hg38, file=output.file, quote=F, row.names=F, col.names=F, sep="\t")
