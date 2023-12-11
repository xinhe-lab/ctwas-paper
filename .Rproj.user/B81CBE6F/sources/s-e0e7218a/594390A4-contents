library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
  stop(" 1 argument:
       * outname", call.=FALSE)
}

print(args)

outname <- as.character(args[1])

#outname <- "/project2/mstephens/causalTWAS/simulations/simulation_ctwas_rss_20210416_compare/ukb-s80.45-adi_simu4-1.Adipose_Subcutaneous.pmr.result_pi_080"

####################

outname_base <- rev(unlist(strsplit(outname, "/")))[1]
outputdir <- paste(rev(rev(unlist(strsplit(outname, "/")))[-1]), collapse="/")

results_files <- list.files(outputdir, outname_base)
results_files <- results_files[grep(".batch", results_files)]
results_files <- paste(outputdir, results_files, sep="/")

results <- lapply(results_files, fread, header=T)

results <- as.data.frame(do.call(rbind, results))

write.table(results, file=outname, row.names=F, col.names=T, sep="\t", quote = F)
