library(ctwas)

args = commandArgs(trailingOnly=TRUE)

# args <- c("/project2/compbio/gwas_summary_statistics/ukbb_neale_v3/ukb-a-360.vcf",
#           "/project2/mstephens/wcrouse/UKB_LDR_0.1",
#           "/project2/mstephens/wcrouse/predictdb_nolnc/mashr_Spleen_nolnc.db",
#           "/project2/mstephens/wcrouse/UKB_analysis_allweights_corrected/ctwas_config.R",
#           "ukb-a-360_Spleen_expr",
#           "ukb-a-360_Spleen_ctwas",
#           "/project2/mstephens/wcrouse/UKB_analysis_allweights_corrected/ukb-a-360/Spleen_nolnc")

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

print(args)

outputdir <- args[7]

dir.create(outputdir, showWarnings=F, recursive=T)

z_snp_stem <- unlist(strsplit(rev(unlist(strsplit(args[1], "/")))[1],"[.]"))[1]
z_snp_outfile <- paste0(outputdir, "/", z_snp_stem, ".RDS")

if (file.exists(z_snp_outfile)){
  z_snp <- readRDS(z_snp_outfile)
} else {
  z_snp <- VariantAnnotation::readVcf(args[1])
  z_snp <- as.data.frame(gwasvcf::vcf_to_tibble(z_snp))
  z_snp$Z <- z_snp$ES/z_snp$SE
  z_snp <- z_snp[,c("rsid", "ALT", "REF", "Z", "SS")]
  colnames(z_snp) <- c("id", "A1", "A2", "z", "ss")
  z_snp <- z_snp[!(z_snp$id %in% z_snp$id[duplicated(z_snp$id)]),] #drop multiallelic variants (id not unique)
  saveRDS(z_snp, file=z_snp_outfile)
}

ld_R_dir <- args[2]
weight <- args[3]
outname.e <- args[5]
outname <- args[6]

#args[4] has been replaced with these hard-coded parameters
thin <- 0.1
ncore <- 10
ncore.rerun <- 1
prob_single <- 0.8
max_snp_region <- 20000
ld_regions <- "EUR"
ld_regions_version <- "b38"

# get gene z score
if (file.exists(paste0(outputdir, "/", outname.e, "_z_gene.Rd"))){
  ld_exprfs <- paste0(outputdir, "/", outname.e, "_chr", 1:22, ".expr.gz")
  load(file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
  load(file = paste0(outputdir, "/", outname.e, "_z_snp.Rd"))
} else {
  res <- impute_expr_z(z_snp, weight = weight, ld_R_dir = ld_R_dir,
                       method = NULL, outputdir = outputdir, outname = outname.e,
                       harmonize_z = T, harmonize_wgt = T,
                       strand_ambig_action_z = "none", recover_strand_ambig_wgt = T)
  z_gene <- res$z_gene
  ld_exprfs <- res$ld_exprfs
  z_snp <- res$z_snp
  
  save(z_gene, file = paste0(outputdir, "/", outname.e, "_z_gene.Rd"))
  save(z_snp, file = paste0(outputdir, "/", outname.e, "_z_snp.Rd"))
}

# run ctwas_rss
ctwas_rss(z_gene, z_snp, ld_exprfs, ld_pgenfs = NULL, ld_R_dir = ld_R_dir, ld_regions = ld_regions, ld_regions_version = ld_regions_version, thin = thin, max_snp_region = max_snp_region, outputdir = outputdir, outname = outname, ncore = ncore, ncore.rerun = ncore.rerun, prob_single = prob_single)
