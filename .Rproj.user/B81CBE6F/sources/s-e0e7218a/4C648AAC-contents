# Combine FUSION assoc test results across chromosomes

library(optparse)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,"--trait",type="character",default=NULL)
parser <- add_option(parser,"--weights",type="character",default="WEIGHTS")
parser <- add_option(parser,"--outdir",type="character",default="FUSION_assoc_test")
parser <- add_option(parser,"--coloc",action = "store_true",default = FALSE)
out    <- parse_args(parser)
trait   <- out$trait
weights <- out$weights
outdir  <- out$outdir
coloc   <- out$coloc
rm(parser,out)


combine_TWAS_table <- function(trait, weights, outdir){

  cat(sprintf("Combine TWAS results for %s using %s weights.\n", trait, weights))

  result_TWAS_noMHC <- data.frame()
  for(CHR in 1:22){

    result_TWAS_chr <- read.table(paste0(outdir, "/", trait, ".", weights, ".", CHR, ".dat"), header = T, stringsAsFactors = F)
    result_TWAS_chr <- result_TWAS_chr[, -c(1,2)]
    result_TWAS_noMHC <- rbind(result_TWAS_noMHC, result_TWAS_chr)
  }

  result_TWAS_MHC <- read.table(paste0(outdir, "/",trait, ".", weights, ".6.dat.MHC"), header = T, stringsAsFactors = F)
  result_TWAS_MHC <- result_TWAS_MHC[, -c(1,2)]
  result_TWAS <- rbind(result_TWAS_noMHC, result_TWAS_MHC)

  result_TWAS_noMHC$TWAS.P.Bonferroni <- p.adjust(result_TWAS_noMHC$TWAS.P, method = "bonferroni")
  result_TWAS$TWAS.P.Bonferroni <- p.adjust(result_TWAS$TWAS.P, method = "bonferroni")

  filename_output <- paste0(outdir, "/", trait, ".", weights, ".noMHC.result")
  write.table(result_TWAS_noMHC, filename_output, col.names = T, row.names = F, quote = F, sep = "\t")
  filename_output <- paste0(outdir, "/", trait, ".", weights, ".result")
  write.table(result_TWAS, filename_output, col.names = T, row.names = F, quote = F, sep = "\t")

}


combine_TWAS_coloc_table <- function(trait, weights, outdir){

  cat(sprintf("Combine TWAS COLOC results for %s using %s weights.\n", trait, weights))

  result_TWAS_noMHC <- data.frame()
  for(CHR in 1:22){
    fn <- paste0(outdir, "/", trait, ".", weights, ".", CHR, ".coloc.dat")
    if (file.exists(fn)){
        result_TWAS_chr <- read.table(fn, header = T, stringsAsFactors = F)
        result_TWAS_chr <- result_TWAS_chr[, -c(1,2)]
        result_TWAS_noMHC <- rbind(result_TWAS_noMHC, result_TWAS_chr)
    }
  }
    
  fn <- paste0(outdir, "/",trait, ".", weights, ".6.coloc.dat.MHC")
  if (file.exists(fn)){
      result_TWAS_MHC <- read.table(fn, header = T, stringsAsFactors = F)
      result_TWAS_MHC <- result_TWAS_MHC[, -c(1,2)]
      result_TWAS <- rbind(result_TWAS_noMHC, result_TWAS_MHC)
  } else{
      result_TWAS <- result_TWAS_noMHC
  }

  result_TWAS_noMHC$TWAS.P.Bonferroni <- p.adjust(result_TWAS_noMHC$TWAS.P, method = "bonferroni")
  result_TWAS$TWAS.P.Bonferroni <- p.adjust(result_TWAS$TWAS.P, method = "bonferroni")

  filename_output <- paste0(outdir, "/", trait, ".", weights, ".coloc.noMHC.result")
  write.table(result_TWAS_noMHC, filename_output, col.names = T, row.names = F, quote = F, sep = "\t")
  filename_output <- paste0(outdir, "/", trait, ".", weights, ".coloc.result")
  write.table(result_TWAS, filename_output, col.names = T, row.names = F, quote = F, sep = "\t")

}

if (coloc) {
  combine_TWAS_coloc_table(trait, weights, outdir)
} else {
  combine_TWAS_table(trait, weights, outdir)
}
