
library(tools)
library(ctwas)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4 ) {
  stop(" 6 arguments (last two optional):
       * expr file (txt file)
       * phenotype file (Rd file)
       * out file name
       * outputdir
       * ncore
       * nsplits", call.= FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "gwas.R"))

exprfs <- read.table(args[1], header = F, stringsAsFactors = F)[,1]
exprvarfs <- sapply(exprfs, prep_exprvar)

phenofile <- args[2]
load(phenofile)
pheno <- phenores$Y

outputdir <- args[4]
outnames <- paste0( args[3], "-B", 1:length(exprfs), ".exprgwas.txt")

ncore <- as.numeric(args[5])
if (is.na(args[5])) ncore <- 1

nsplits <- as.numeric(args[6])
if (is.na(args[6])) nsplits <- 10

for (b in 1:length(exprfs)){

  outname <- outnames[b]

  # load expr data

  geneinfo <- read_exprvar(exprvarfs[b])
  anno <- geneinfo[, - "id"] # chrom p0 p1

  GWAA(exprfs[b], mode = "expr", pheno, snpname = geneinfo$id, anno = anno, outname, outputdir, family = gaussian, ncore = ncore, nSplits = nsplits, compress = T)

}

# combine all results
outdflist <- list()

for (b in 1:length(exprfs)){
  outdflist[[b]] <- read.table(paste0(outnames[b], ".gz"),
                               header = T, stringsAsFactors = F, comment.char = "")
}
outdf <- do.call(rbind, outdflist)

outname <-  paste0(outputdir, "/", args[3], ".exprgwas.txt")
write.table(outdf, file = outname , row.names=F, col.names=T, sep="\t", quote = F)
system(paste0("bgzip ", outname))

# clean
for (b in 1:length(exprfs)){
  system(paste0("rm ", outnames[b], ".gz"))
}




