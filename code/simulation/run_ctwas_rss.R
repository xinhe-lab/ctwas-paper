library(ctwas)

args = commandArgs(trailingOnly=TRUE)
if (length(args) <7) {
  stop(" 7 arguments:
       * zscore file name
       * ld genotype file (multiple file names in a .txt file)
       * weight
       * config file name (.R)
       * out expr z file name
       * out file name
       * outputdir", call.=FALSE)
}

z_snp <- data.table::fread(args[1], header = T)
data.table::setnames(z_snp, old = c("alt", "ref", "t.value"), new = c("A1", "A2", "z"))
z_snp <- z_snp[, c("id", "A1", "A2", "z")]
ld_pgenfs <- read.table(args[2], header = F, stringsAsFactors = F)[,1]
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
    res <- impute_expr_z(z_snp, weight = weight, ld_pgenfs = ld_pgenfs,
                         method = "lasso", outputdir = outputdir,
                         outname = outname.e)
    z_gene <- res$z_gene
    ld_exprfs <- res$ld_exprfs

    save(z_gene, file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
}

# save.image("~/temp.Rd")
# run ctwas_rss
ctwas_rss(z_gene, z_snp, ld_exprfs, ld_pgenfs = ld_pgenfs, ld_R_dir = NULL, ld_regions = "EUR", ld_regions_custom = ld_regions_custom, thin = thin, max_snp_region = max_snp_region, outputdir = outputdir, outname = outname, ncore = ncore, ncore.rerun = ncore.rerun, prob_single = prob_single, harmonize_z = F)

