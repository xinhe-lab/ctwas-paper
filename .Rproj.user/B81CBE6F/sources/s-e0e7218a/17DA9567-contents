library(ctwas)

pgenfs <- paste0("/project2/mstephens/causalTWAS/ukbiobank/ukb_pgen_s40.22/ukb_chr", 1:22, "_s40.22_thinned.pgen")
weight <- "/project2/mstephens/causalTWAS/fusion_weights/Adipose_Subcutaneous"
load("/project2/mstephens/causalTWAS/simulations/simulation_ashtest_20201001/simu_20201001-1-1-pheno.Rd")
Y <- phenores$Y
s1 <- data.table::fread("/project2/mstephens/causalTWAS/ukbiobank/ukb_pgen_s40.22/ukb_chr10_s40.22_thinned.psam")
s2 <- data.table::fread("/project2/mstephens/causalTWAS/ukbiobank/ukbiobank_samples_s40.22.txt")
Y <- Y[match(s1$IID, s2$V1),]

outputdir <- "/project2/mstephens/causalTWAS/simulations/test_package_temp"

exprfs <- vector()
for (i in 1:22){
  pgenf <- pgenfs[i]
  exprfs[i] <- impute_expr(pgenf = pgenf, weight = weight,
                           method = "lasso", outputdir = outputdir,
                           outname = "test")
}

ctwas(pgenfs, exprfs, Y, ld_regions = "EUR", thin = 1, outputdir = outputdir, outname = "test", ncore = 14)
