library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop(" 6 arguments:
       * prune ld score file
       * gwas file
       * eqtl summary statistics
       * weight dir
       * TWAS results files
       * out file name", call.=FALSE)
}

mrjti <- "/home/simingz/causalTWAS/MR-JTI/mr/MR-JTI.r"

ldfile <- as.character(args[1])
gwas <- as.character(args[2])
eqtl <- as.character(args[3])
weightf <- as.character(args[4])
twas <- as.character(args[5])
outname <- as.character(args[6])

####################

harmonize_eqtl_gwas <- function(eqtl_sumstats, gwas_sumstats){
  snpnames <- intersect(eqtl_sumstats$id, gwas_sumstats$id)

  if (length(snpnames) != 0){
    ss.idx <- which(eqtl_sumstats$id %in% snpnames)
    ld.idx <- match(eqtl_sumstats$id[ss.idx], gwas_sumstats$id)
    qc <- ctwas:::allele.qc(eqtl_sumstats$alt[ss.idx], eqtl_sumstats$ref[ss.idx],
                            gwas_sumstats$alt[ld.idx], gwas_sumstats$ref[ld.idx])
    ifflip <- qc[["flip"]]
    ifremove <- !qc[["keep"]]

    flip.idx <- ss.idx[ifflip]
    eqtl_sumstats[flip.idx, c("alt", "ref")] <- eqtl_sumstats[flip.idx, c("ref", "alt")]
    eqtl_sumstats[flip.idx, "slope"] <- -eqtl_sumstats[flip.idx, "slope"]

    remove.idx <- ss.idx[ifremove]
    if (length(remove.idx) != 0) {
      eqtl_sumstats <- eqtl_sumstats[-remove.idx,,drop = F]
    }
  }
  return(eqtl_sumstats)
}

####################
#load weight names
weights <- as.data.frame(fread(weightf, header = T))
weights$ENSEMBL_ID <- sapply(weights$WGT, function(x){unlist(strsplit(unlist(strsplit(x,"/"))[2], "[.]"))[2]})

# get genes that pass TWAS FDR cut off
a <- fread(twas, header = T)
a$TWAS.P.FDR <- p.adjust(a$TWAS.P, method = "fdr")
genes <- a[a$TWAS.P.FDR <0.05, ID]
rm(a)

genes <- weights$ENSEMBL_ID[match(genes, weights$ID)]
rm(weights)

####################
# GWAS
gwas <- as.data.frame(fread(gwas, header = T))

#drop variants that are non uniquely identified by ID
gwas <- gwas[!(gwas$id %in% gwas$id[duplicated(gwas$id)]),]

####################
#load eQTL; in GTEx, the eQTL effect allele is the ALT allele
eqtl <- as.data.frame(fread(eqtl, header = T))

#trim version number from ensembl IDs
eqtl_gene_id_crosswalk <- unique(eqtl$gene_id)
eqtl_gene_id_crosswalk <- data.frame(original=eqtl_gene_id_crosswalk,
                                     trimmed=sapply(eqtl_gene_id_crosswalk, function(x){unlist(strsplit(x, "[.]"))[1]}))
eqtl$gene_id <- eqtl_gene_id_crosswalk[eqtl$gene_id, "trimmed"]
eqtl <- eqtl[, c("rs_id_dbSNP147_GRCh37p13", "gene_id", "gene_name", "ref", "alt", "slope", "slope_se", "pval_nominal")]
eqtl <- dplyr::rename(eqtl, id="rs_id_dbSNP147_GRCh37p13")

#keep genes detected by TWAS FDR
eqtl <- eqtl[eqtl$gene_id %in% genes,]

#drop entries not uniquely identified by gene_id and variant id (variant not biallelic)
eqtl_unique_id_gene <- paste0(eqtl$id, eqtl$gene_id)
eqtl <- eqtl[!(eqtl_unique_id_gene %in% eqtl_unique_id_gene[duplicated(eqtl_unique_id_gene)]),]

# harmonize gwas and eqtl
eqtl <- harmonize_eqtl_gwas(eqtl, gwas)

# merge gwas, eqtl and ldinfo
ldfiles <- as.data.frame(fread(ldfile, header = F))

ldinfolist <- list()
for (i in 1:22){
  prune <- as.data.frame(fread(ldfiles[[i,1]], header = F))
  colnames(prune) <- "SNP"
  ldsc <- as.data.frame(fread(ldfiles[[i,2]], header = T))
  ldinfolist[[i]] <- merge(prune, ldsc, by = "SNP")
}
ldinfo <- do.call(rbind, ldinfolist)

snpall <- merge(merge(gwas, eqtl, by.x = "id", by.y = "id"),
      ldinfo, by.x = "id", by.y = "SNP")

rm(eqtl)

####################

# count number of genes with >=20 SNPs for multiple testing correction
ngene <- sum(table(snpall$gene_id[snpall$gene_id %in% genes])>=20)

# run MR-JTI
gnlist <- list()
gene.test <- NULL
reslist <- list()
for (gene in genes){
  gout <- snpall[snpall$gene_id == gene, c("id", "ldscore", "slope", "slope_se", "pval_nominal", "Estimate", "PVALUE")]
  colnames(gout) <- c("rsid", "ldscore", "eqtl_beta", "eqtl_se", "eqtl_p", "gwas_beta", "gwas_p")
  gnlist[[gene]] <- nrow(gout)
  if (nrow(gout) < 20) next
  gfile <- paste0(outname, "_", gene,"_input.txt")
  resfile <- paste0(outname, "_", gene, "_res.txt")
  write.table(gout, file= gfile, row.names=F, col.names=T, sep="\t", quote = F)
  gene.test<- c(gene.test, gene)
  system(paste("Rscript", mrjti,  "--df_path", gfile, "--n_genes", ngene, "--out_path", resfile))

  if (file.exists(resfile)){
    res <- read.table(resfile, header = T, stringsAsFactors = F, sep = ",")
    res[1, 1] <- gene
    reslist[[gene]] <- res[1, ]
    file.remove(resfile)
  }

  file.remove(gfile)
}

outdf <- do.call(rbind, reslist)
write.table(outdf, file=outname, row.names=F, col.names=T, sep="\t", quote = F)
