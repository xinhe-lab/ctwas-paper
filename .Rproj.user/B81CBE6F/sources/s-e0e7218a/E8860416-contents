library(ctwas)

args = commandArgs(trailingOnly=TRUE)
if (length(args) <7) {
  stop(" 7 arguments:
       * zscore file name
       * ld R directory
       * weight
       * config file name (.R)
       * out expr z file name
       * out file name
       * outputdir", call.=FALSE)
}

z_snp <- data.table::fread(args[1], header = T)
data.table::setnames(z_snp, old = c("alt", "ref", "t.value"), new = c("A1", "A2", "z"))
z_snp <- z_snp[, c("id", "A1", "A2", "z")]
ld_R_dir <- args[2]
weight <- args[3]
outname.e <- args[5]
outname <- args[6]
outputdir <- args[7]


ld_regions_custom <- "/home/simingz/ctwas/inst/extdata/example_regions.bed"
thin <- 1
ncore <- 1
ncore.rerun <- 1
prob_single <- 0.8
thin <- 1
max_snp_region <- 1e4

source(args[4]) # config
# get gene z score
if (file.exists(paste0(outputdir, "/", outname.e, "_z_gene.Rd"))){
  ld_exprfs <- paste0(outputdir, "/", outname.e, "_chr", 1:22, ".expr.gz")
  load(file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
} else {
    res <- impute_expr_z(z_snp, weight = weight, ld_R_dir = ld_R_dir,
                         method = "lasso", outputdir = outputdir,
                         outname = outname.e)
    z_gene <- res$z_gene
    ld_exprfs <- res$ld_exprfs

    save(z_gene, file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
}

# save.image("temp.Rd")
# run ctwas_rss
ctwas_rss(z_gene, z_snp, ld_exprfs, ld_pgenfs = NULL, ld_R_dir = ld_R_dir, ld_regions = "EUR", ld_regions_custom = ld_regions_custom, thin = thin, max_snp_region = max_snp_region, outputdir = outputdir, outname = outname, ncore = ncore, ncore.rerun = ncore.rerun, prob_single = prob_single, harmonize = T)

