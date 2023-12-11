library(ctwas)
library(logging)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("5 arguments must be supplied:
       * genotype file name (multiple file names in a .txt file)
       * imputed expr file (multiple file names in a .txt file)
       * param R file
       * out file name
       * outputdir", call.=FALSE)
}

codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir, "stats_func.R"))       
source("~/causalTWAS/causal-TWAS/code/simulate_phenotype.R")

pgenfs <- read.table(args[1], header = F, stringsAsFactors = F)[,1]
exprfs <- read.table(args[2], header = F, stringsAsFactors = F)[,1]

prior_dist_causal <- "normal"
source(args[3])
outname <- args[4]
outputdir <- args[5]

loginfo("simulating phenotype for %s", outname)
set.seed(SED)
phenores <- simulate_phenotype(pgenfs, exprfs,
                               pve.expr = pve.expr,
                               pve.snp = pve.snp,
                               pi_beta =  pi_beta,
                               pi_theta = pi_theta,
                               prior_dist_causal = prior_dist_causal)

save(phenores, file = paste0(outputdir,"/", outname, "-pheno.Rd"))

