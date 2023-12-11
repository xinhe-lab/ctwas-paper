library(ctwas)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("4 arguments must be supplied:
       * genotype file name (multiple file names in a .txt file)
       * weight dir
       * out file name
       * outputdir", call.=FALSE)
}


pgenfs <- read.table(args[1], header = F, stringsAsFactors = F)[,1]
weight <- args[2]

outname <- args[3]
outputdir <- args[4]

exprfs <- impute_expr(pgenfs = pgenfs, weight = weight,
                           method = "lasso", outputdir = outputdir,
                           outname = outname)

out <- data.frame("file" = exprfs, stringsAsFactors = F)
write.table(out, file = paste0(outputdir, "/", outname, ".expr.txt"),
             row.names=F, col.names=F, sep="\t", quote = F)

