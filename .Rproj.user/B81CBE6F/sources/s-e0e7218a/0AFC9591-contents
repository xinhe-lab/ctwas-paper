codedir <- "/project2/mstephens/causalTWAS/causal-TWAS/code/"
source(paste0(codedir,"input_reformat.R"))

# Fusion weights
weight <- "/home/simingz/causalTWAS/fusion_weights/Adipose_Subcutaneous"
method <- "lasso"
wgtcut <- 1e-8 # abs(wight) > wgtcut

wgtfs <- read.table(paste0(weight,".pos"), header = T)
wgtdir <-dirname(weight)
wgtfs <- paste0(wgtdir, "/", wgtfs$WGT)
snpall <- NULL
snpannoall <- NULL
for (wgtf in wgtfs){
  load(wgtf)
  snp <- rownames(wgt.matrix[abs(wgt.matrix[, method]) > wgtcut,])
  snpanno <- snps[snps[,2] %in% snp, ]
  snpanno <- cbind(snpanno, wgt.matrix[snpanno[,2],])
  snpall <- c(snpall, snp)
  snpannoall <- rbind(snpannoall, snpanno)
}

write.table(snpannoall, file= paste(weight, method, "eqtl.txt", sep= ".") , row.names=F, col.names=T, sep="\t", quote = F)

#----------------------------chr17-22 combined-------------------
eqtl <- "/home/simingz/causalTWAS/fusion_weights/Adipose_Subcutaneous.lasso.eqtl.txt"
prune <- c(T, rep(F,9)) # T will be kept
pfiles <- paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", 17:22, "_s20000.pgen")

snpall <- read.table(eqtl, header = T, stringsAsFactors = F)[, 2]
mtxall <- NULL
mtxall.unscaled <- NULL
varall <- NULL
for (pfile in pfiles){
  pfile <- drop_ext(pfile)
  pvar <- paste0(pfile, ".pvar")
  var <- read.table(pvar, header = F)
  eqtl <- var[,3] %in% snpall
  prunetag <- rep_len(prune, nrow(var))
  select <- (1:nrow(var))["|"(eqtl, prunetag)]
  varall <- rbind(varall, var[select,])
  mtx <- pgen2mtx(pfile, nb_parts = 1, select=select, progress = T)
  mtxall.unscaled <- cbind(mtxall.unscaled, mtx)
  # mtx.scaled <- scaleRcpp(mtx)
  # mtxall <- cbind(mtxall, mtx.scaled)
}

# dat <- list("G"       = mtxall,
#             "chr"     = as.matrix(varall[,1]),
#             "pos"     = as.matrix(varall[,2]),
#             "snp"     = as.matrix(varall[,3]),
#             "counted" = as.matrix(varall[,4]),
#             "alt"     = as.matrix(varall[,5]))
#
#
# m <- as_FBM(mtxall, backingfile = "/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000", is_read_only = T)$save()
# dat$G <- m
# save(dat, file = "/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000.FBM.Rd")

dat <- list("G"       = mtxall.unscaled,
            "chr"     = as.matrix(varall[,1]),
            "pos"     = as.matrix(varall[,2]),
            "snp"     = as.matrix(varall[,3]),
            "counted" = as.matrix(varall[,4]),
            "alt"     = as.matrix(varall[,5]))

m <- as_FBM(mtxall.unscaled, backingfile = "/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000.unscaled", is_read_only = T)$save()
dat$G <- m
save(dat, file = "/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000.unscaled.FBM.Rd")

#----------------------------chr22 only-------------------------
pfile <- "/home/simingz/causalTWAS/ukbiobank/ukb_chr22_s20000"
pgen2fbm(pfile, select = NULL, scale = T, type = "double")


#----------------------------whole genome ----------------------
pgen2dat <- function(eqtl, prune, pfiles, outf){
  snpall <- read.table(eqtl, header = T, stringsAsFactors = F)[, 2]
  mtxall <- NULL
  mtxall.unscaled <- NULL
  varall <- NULL
  for (pfile in pfiles){
    pfile <- drop_ext(pfile)
    pvar <- paste0(pfile, ".pvar")
    var <- read.table(pvar, header = F)
    eqtl <- var[,3] %in% snpall
    prunetag <- rep_len(prune, nrow(var))
    select <- (1:nrow(var))["|"(eqtl, prunetag)]
    varall <- rbind(varall, var[select,])
    mtx <- pgen2mtx(pfile, nb_parts = 50, select=select, progress = T)
    mtxall.unscaled <- cbind(mtxall.unscaled, mtx)
    #mtx.scaled <- scaleRcpp(mtx)
    #mtxall <- cbind(mtxall, mtx.scaled)
  }

  dat <- list(
              #"G"       = mtxall,
              "G" = mtxall.unscaled,
              "chr"     = as.matrix(varall[,1]),
              "pos"     = as.matrix(varall[,2]),
              "snp"     = as.matrix(varall[,3]),
              "counted" = as.matrix(varall[,4]),
              "alt"     = as.matrix(varall[,5]))


  m <- as_FBM(mtxall.unscaled, backingfile = outf, is_read_only = T)$save()
  dat$G <- m
  save(dat, file = paste0(outf, ".FBM.Rd"))
}

eqtl <- "/home/simingz/causalTWAS/fusion_weights/Adipose_Subcutaneous.lasso.eqtl.txt"
prune <- c(T, rep(F,9)) # T will be kept


for (chrom in 1:22){
  pfiles <- paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", chrom, "_s40000.pgen")
  outf <- paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", chrom, "_s40000")
  pgen2dat(eqtl, prune, pfiles, outf)
}

for (chrom in 1:22){
  pfiles <- paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", chrom, "_s40000.pgen")
  outf <- paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", chrom, "_s40000.unscaled")
  pgen2dat(eqtl, prune, pfiles, outf)
}

# regions file
chunksize <- 5e5
regions <- NULL
for (chrom in 1:22){
  load(paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", chrom, "_s40000.FBM.Rd"))
  stpos <- min(dat$pos) - chunksize/2
  edpos <- max(dat$pos) + chunksize/2
  p0 <- stpos + 0: floor((edpos - stpos)/chunksize) * chunksize
  p1 <- stpos + 1: ceiling((edpos - stpos)/chunksize) * chunksize
  itv <- cbind(chrom, p0, p1)
  regions <- rbind(regions, itv)
  regfile <- paste0("/home/simingz/causalTWAS/simulations/shared_files/chr", chrom, "-500kb_bins.txt")
  write.table(itv, file= regfile , row.names=F, col.names=T, sep="\t", quote = F)
}
write.table(regions, file= "/home/simingz/causalTWAS/simulations/shared_files/chr1-22-500kb_bins.txt",
            row.names=F, col.names=T, sep="\t", quote = F)

#-------------------------filter s40000 samples---------------------------------------
s.s <- read.table("/home/simingz/causalTWAS/ukbiobank/ukbiobank_samples_s40.22.txt", header =F, stringsAsFactors = F)
s.ori <- read.table("/home/simingz/causalTWAS/ukbiobank/ukbiobank_samples40000.txt", header =F, stringsAsFactors = F)

s.s.idx <- match(s.s[,1], s.ori[,1])
for (chrom in 1:22){
  outf <- paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", chrom, "_s40000.unscaled.FBM.Rd")
  load(outf)
  G.new <- dat$G[]
  G.new <- G.new[s.s.idx, ]
  m <- as_FBM(G.new, backingfile = paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", chrom, "_s40.22.unscaled"), is_read_only = T)$save()
  dat$G <- m
  save(dat, file = paste0("/home/simingz/causalTWAS/ukbiobank/ukb_chr", chrom, "_s40.22.unscaled.FBM.Rd"))
}

#---------chr17-22 multiple copies with different names---
setwd("/home/simingz/causalTWAS/ukbiobank")
for (i in 1:6) {
  load("/home/simingz/causalTWAS/ukbiobank/ukb_chr17to22_s20000.FBM.Rd")
  tag <- paste0("-B", i)
  dat$chr <- as.matrix(paste0(dat$chr, tag), ncol = 1)
  dat$snp <- as.matrix(paste0(dat$snp, tag), ncol = 1)
  save(dat, file = paste0("ukb_chr17to22_s20000", tag, ".FBM.Rd"))
}

# weight files
setwd("/home/simingz/causalTWAS/fusion_weights")
posf0 <- "Adipose_Subcutaneous.pos"
posf = "Adipose_Subcutaneous_B.pos"
pos0 = read.table(posf0, header = T, stringsAsFactors = F)
pos <- read.table(posf, header = T, stringsAsFactors = F)

pos0 <- pos0[pos0$CHR <=22 & pos0$CHR >=17,]
pos <- pos[pos$CHR <=22 & pos$CHR >=17,]

renamesnps <- function(row, tag) {
  f0 <- row[1]
  f <- row[2]
  load(f0)
  rownames(wgt.matrix) <- paste0(rownames(wgt.matrix), tag)
  snps[,1] <- paste0(snps[,1], tag)
  snps[,2] <- paste0(snps[,2], tag)
  save(wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot, file = f)
}


outdfall <- NULL
for (i in 1:6) {
  tag <- paste0("-B", i)
  outdf <- pos
  outdf$WGT <- gsub(".wgt", paste0(tag, ".wgt"), outdf$WGT)
  outdf$ID <- paste0(outdf$ID, tag)
  outdf$CHR <- paste0(outdf$CHR, tag)
  rd <- cbind(pos0$WGT, outdf$WGT)
  apply(rd, 1, renamesnps, tag = tag)
  outdfall <- rbind(outdfall, outdf)
}

write.table(outdfall , file= posf, row.names=F, col.names=T, sep="\t", quote = F)

# region files
setwd("/home/simingz/causalTWAS/simulations/shared_files")
regfile0 <- "chr17-22-500kb_bins.txt"
for (i in 1:6) {
  regfile <- paste0("chr17-22-500kb_B", i, "_bins.txt")
  tag <- paste0("-B", i)
  reg <- read.table(regfile0, header =T)
  reg[,1] <- paste0(reg[,1], tag)
  write.table(reg, file= regfile , row.names=F, col.names=T, sep="\t", quote = F)
}





