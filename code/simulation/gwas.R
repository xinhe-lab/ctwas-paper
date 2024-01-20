library(foreach)
library(doParallel)
library(ctwas)

#' @description this is a function to run GWAS analysis for small genotype matrix (using glm),for faster results use plink (start from plink files). Credit to http://www.stat-gen.org/tut/tut_analysis.html, requires tabix and bgzip if compress is true.
#'
#' @param geno genotype matrix, columns: snps, rows: individuals.
#' @param pheno phenotype, matrix, ncol = 1
#' @param snpname a vector of N, correpond to geno type matrix, will be written to output.
#' @param anno a data frame/matrix, anno column with N rows.
#' @param compress T/F, if T will be compressed and indexed by tabix. (TODO)
#'
GWAA <- function(genof, mode = c("snp", "expr"), pheno, snpname, anno = NULL, outname, outputdir = getwd(), family = gaussian, ncore = 1, nSplits = 10, compress = F){

  cl <- makeCluster(ncore)
  show(cl)
  registerDoParallel(cl)

  colnames(pheno) <- "pheno"

  mode <- match.arg(mode)

  columns<- c(colnames(anno), "id", c("Estimate", "Std.Error", "t-value", "PVALUE"))
  write.table(t(columns), paste0(outputdir, "/",outname), row.names=FALSE, col.names=FALSE, quote=FALSE, sep = "\t")
  nSNPs <- length(snpname)
  if (nSNPs > 0){

    genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset

    snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
    snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group
    res <- foreach(part = 1:length(snp.start),.combine='rbind', .packages = "ctwas")  %dopar% {

      cat(sprintf("GWAS %s%% finished\n", part/length(snp.start)))

      get_geno <- function(genof, mode = mode,
                           variantidx = NULL, outputdir = NULL){
        if (mode == "expr"){
          geno <- ctwas:::read_expr(genof, variantidx = variantidx)
        } else if (mode == "snp"){
          pvarf <- ctwas:::prep_pvar(genof, outputdir = outputdir)
          pgen <- ctwas:::prep_pgen(genof, pvarf)
          geno <- ctwas:::read_pgen(pgen, variantidx = variantidx)
        } else {
          stop("undefined mode")
        }
        geno
      }
      
      geno.i <- get_geno(genof, mode = mode,
                         variantidx = snp.start[part]:snp.stop[part],
                         outputdir = outputdir)

      snpname.i <- snpname[snp.start[part]:snp.stop[part]]
      if (!is.null(anno)) {
        anno.i <- anno[snp.start[part]:snp.stop[part], ]
      } else {
        anno.i <- NULL
      }

      res.i <- matrix(, nrow = length(snpname.i), ncol =4)
      for (snp.idx in 1: length(snpname.i)) {
        snp <- geno.i[, snp.idx, drop = F]
        a <- summary(glm(pheno ~ snp,
                         family=family,
                         data=data.frame("pheno" = pheno, "snp" = snp)))
        res.i[snp.idx, ] <- a$coefficients['snp',]
      }
      
      colnames(res.i) <- c("Estimate", "Std.Error", "t-value", "PVALUE")

      res.out <- cbind(anno.i, snpname.i, res.i)
      res.out
    }
      # write results so far to a file
    write.table(res, paste0(outputdir, "/",outname),
                  append=TRUE,
                  quote=FALSE, col.names=FALSE, row.names=FALSE,
                  sep = "\t")


  }

  stopCluster(cl)

  # compress
  if (isTRUE(compress)){
    system(paste0("gzip -f ", outname))
  }

  return(print("Done."))
}

