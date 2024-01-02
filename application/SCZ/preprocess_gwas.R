#!/software/R-4.0.4-el7-x86_64/bin/R
args=commandArgs(trailingOnly=TRUE)


library(data.table)
library(bigsnpr)
library(tibble)
library(dplyr)

source("~/scripts/R/utility_func.R")
inDIR="/home/jinggu/cluster/data/gwas/1_Processed"
outDIR="/home/jinggu/cluster/data/gwas/2_Filtered"
#functions:
# merge SNP index information from reference panels to gwas
# flip strand if gwas effect alleles mismatch with reference allele
# snp_match func. automatically puts an opposite sign if strand is flipped
merge.bigsnp.gwas <- function(gwas, bigSNP){
  
  map <- bigSNP$map
  snp_info <- map[,c('chromosome','physical.pos','allele1','allele2')]
  colnames(snp_info) <- c('chr','pos','a0','a1')
  
  matched.gwas <- as_tibble(bigsnpr::snp_match(gwas, 
                                               snp_info, 
                                               strand_flip = T, 
                                               match.min.prop = 0.1)) %>% dplyr::rename(og_index = `_NUM_ID_.ss`) %>% dplyr::rename(bigSNP_index = `_NUM_ID_`) %>% mutate(zscore = beta/se)
  
  matched.gwas["rsID"]= map[matched.gwas$bigSNP_index,"marker.ID"]
  return(matched.gwas)
}

##load LD blocks
bigSNP <- bigsnpr::snp_attach(rdsfile = '/project2/xinhe/1kg/bigsnpr_QC/1000G_EUR_Phase3_QC.rds')


#load formatted GWAS summary stats.
gwas<-fread(args[1])
#the order of a0, a1 doesn't matter as long as they are consistent between zscore and annotation files
colnames(gwas)<-c('chr', 'pos', 'beta','se','a1', 'a0','rawSNP', 'pval')


#drop XY or mitochondrial chromosomes
gwas<-gwas[!(gwas$chr %in% c("X", "Y", "M"))]

#compute Zscores
gwas[,"zscore"]<-gwas$beta/gwas$se
gwas<-gwas[!is.na(gwas$zscore),]

#convert alleles to upper case
gwas$a0<-toupper(gwas$a0)
gwas$a1<-toupper(gwas$a1)
#remove indels
snplist<-c('A', 'T', 'C', 'G')
a0IsSNP<-(gwas$a0 %in% snplist)
a1IsSNP<-(gwas$a1 %in% snplist)
gwas<-gwas[a0IsSNP & a1IsSNP,]

# drop duplicate SNPs
chrpos <- paste0(gwas$chr, '_', gwas$pos)
gwas <- gwas[!duplicated(chrpos), ]
print('duplicate snps removed')

gwas$chr<-as.integer(gwas$chr)

gwas<-merge.bigsnp.gwas(gwas, bigSNP)
print('Assining SNPs to LD blocks...')
#load European independent LD blocks
load("~/scripts/R/Euro_LD_Chunks.RData")
#assign SNPs to LD blocks
gwas <- assign.locus.snp(cleaned.sumstats = gwas, ld = LD_Blocks)
print('Complete cleaning')
print(paste("Total number of GWAS SNPs after filtering", dim(gwas)[1], sep=":"))

#change snp format to be chr:pos:a1:a0:rsID
gwas["snp"]=paste(as.character(gwas$chr), as.character(gwas$pos), gwas$a1, gwas$a0, gwas$rsID, sep=":")
#save cleaned summary stats(gwas, args[2])
write.table(gwas, gzfile(args[2]), sep="\t", quote = FALSE, row.names=FALSE)
write.table(gwas[, c("snp", "locus", "zscore")], gzfile(args[3]), sep='\t', quote = FALSE, row.names=FALSE)
