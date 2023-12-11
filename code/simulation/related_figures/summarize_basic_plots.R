library(ggplot2)
library(cowplot)
library(plotrix)
library(ggpubr)

# Calculate power in a GWAS study
pow <- function(total, n, beta, cutp){
  rec <- rep(0, total)
  for (i in 1:total){
    x <- rnorm(n)
    y <- x * rnorm(1, sd = beta) + rnorm(n, sd = sqrt(2.5))
    lm.s <- lm(y~x)
    pv <- summary(lm.s)$coefficients[2,4]
    rec[i] <- pv
  }
  length(rec[rec < cutp])/length(rec)
}

# p value distribution
pdist_plot <- function(gwasf, chr, cau, thin = 1){

  a <- data.table::fread(gwasf, header = T)
  a <- a[a$chrom == chr]

  a$ifcausal <- ifelse(a$id %in% cau, 1, 0)
  a.causal <- a[a$ifcausal == 1, ]
  a <- a[sample(nrow(a), nrow(a)* thin), ]
  a.nonc <- a[a$ifcausal == 0, ]
  a <- rbind(a.causal, a.nonc)

  ax <- pretty(0:max(-log10(a$PVALUE)), n = 30)

  if (!("p0" %in% colnames(a))) a$p0 <- a$pos

  par(mfrow=c(3,1))

  h1 <- hist(-log10(a$PVALUE), breaks = 100, xlab = "-log10(p)", main = "P value distribution-all", col = "grey", xlim= c(3,20), ylim =c(0,50)); grid()

  h2 <- hist(-log10(a[a$ifcausal == 1, ]$PVALUE), breaks = h1$breaks, xlab = "-log10(p)", main = "P value distribution-causal", col = "salmon", xlim= c(3,20), ylim =c(0,50));grid()

  plot(a$p0, -log10(a$PVALUE), col = a$ifcausal + 1, xlab = paste0("chr", chrom), ylab = "-log10(pvalue)")
  points(a[a$ifcausal ==1, ]$p0, -log10(a[a$ifcausal ==1, ]$PVALUE), col = "red", pch =19)
  grid()

  # cat("number of p < ", p.cut, ": ",nrow(a[a$PVALUE < p.cut, ]))
  # cat("number of causal p < ", p.cut, ": ", nrow(a[a$PVALUE < p.cut & a$ifcausal ==1, ]))

}

# used by show_param, caliPIP_plot
get_causal_id <- function(phenores){
  gene.causal <- unlist(sapply(1:22, function(x) phenores$batch[[x]]$id.cgene))
  snp.causal <- unlist(sapply(1:22, function(x) phenores$batch[[x]]$id.cSNP))
  return(c(gene.causal, snp.causal))
}

# used by nca_plot
.obn <- function(pips, ifcausal, mode = c("PIP", "FDR")){
  if (mode == "PIP"){
    a_bin <- cut(pips, breaks= seq(0, 1, by=0.1))
  } else if (mode == "FDR"){
    a_bin <- cut(pips, breaks= c(0, 0.01, 0.05, 0.1, seq(0.2, 1, by=0.1)))
  }

  obca = c(by(ifcausal, a_bin, FUN = sum))
  obnon = c(by((1-ifcausal), a_bin, FUN = sum))
  return(list("ncausal" = obca,
              "nnoncausal" = obnon))
}

# used by ncausal plot
nca_plot <- function(pips, ifcausal, runtag = NULL, mode = c("PIP", "FDR"), xmin =0, main = mode[1], ...){
  # ifcausal:0,1, runtag: for adding std.
  if (is.null(runtag)){
    runtag <- rep(1, length(pips))
  }

  if (mode == "PIP"){
    bins <- seq(0, 1, by=0.1)[1:10]
  } else if (mode == "FDR"){
    bins <- c(0, 0.01, 0.05, 0.1, seq(0.2, 1, by=0.1))[1:12]
  }

  calist <- list()
  nonlist <- list()
  for (rt in unique(runtag)){
    pips.rt <- pips[runtag == rt]
    ifcausal.rt <- ifcausal[runtag == rt]
    res <- .obn(pips.rt, ifcausal.rt, mode = mode)
    calist[[rt]] <- cbind(res[["ncausal"]], "causal", bins)
    nonlist[[rt]] <- cbind(res[["nnoncausal"]], "noncausal", bins)
  }

  df <- rbind(do.call(rbind, calist),
          do.call(rbind, nonlist))
  df <- data.frame("count"= as.numeric(df[,1]),
                   "ifcausal" = factor(df[,2], levels = c("noncausal", "causal")),
                   "bins" = as.numeric(df[,3]))

  if (mode == "PIP"){
    ymax <- 1.1* (max(df[df$bins > xmin & df$ifcausal == "causal", "count"], na.rm = T) + max(df[df$bins > xmin & df$ifcausal == "noncausal", "count"], na.rm = T))
  } else {
    ymax <- 1.1* (max(df[df$bins < xmin & df$ifcausal == "causal", "count"], na.rm = T) + max(df[df$bins < xmin & df$ifcausal == "noncausal", "count"], na.rm = T))
  }


  fig <- ggbarplot(df, x = "bins", y = "count", add = "mean_se", fill = "ifcausal", palette = "jco", ylim=c(0,ymax), main =main)

  return(fig)
}

# used by cp_plot
.exob <- function(pips, ifcausal, mode = c("PIP", "FDR")){
  a_bin <- cut(pips, breaks= seq(0, 1, by=0.1))
  if (mode == "PIP") {
    ex = c(by(pips, a_bin, FUN = mean))
    ob = c(by(ifcausal, a_bin, FUN = mean))
  } else if (mode == "FDR"){
    ex = c(by(pips, a_bin, FUN = mean))
    ob = 1 - c(by(ifcausal, a_bin, FUN = mean))
  }
  return(list("expected" = ex, "observed" = ob))
}

# used by caliPIP_plot, caliFDP_plot
cp_plot <- function(pips, ifcausal, runtag = NULL, mode = c("PIP", "FDR"), main = mode[1]){

  # ifcausal:0,1, runtag: for adding std.
  if (is.null(runtag)){
    se = 0
  } else{
    dflist <- list()
    for (rt in unique(runtag)){
      pips.rt <- pips[runtag == rt]
      ifcausal.rt <- ifcausal[runtag == rt]
      dflist[[rt]] <- .exob(pips.rt, ifcausal.rt, mode = mode)
    }
    mean_pip <- colMeans(do.call(rbind, lapply(dflist, '[[', "expected")),na.rm = T)
    observed_freq <- colMeans(do.call(rbind, lapply(dflist, '[[', "observed")), na.rm =T)
    se <- apply(do.call(rbind, lapply(dflist, '[[', "observed")), 2, plotrix::std.error, na.rm =T)
  }

  df <- data.frame("mean_pip" = mean_pip, "observed_freq"= observed_freq, "se" = se)

  ggplot(df, aes(x=mean_pip, y=observed_freq)) +
    geom_errorbar(aes(ymin=observed_freq-se, ymax=observed_freq+se), colour="black", size = 0.5, width=.01) +
    geom_point(size=1.5, shape=21, fill="#002b36") + # 21 is filled circle
    xlab("Expected") +
    ylab("Observed") +
    coord_cartesian(ylim=c(0,1), xlim=c(0,1)) +
    geom_abline(slope=1,intercept=0,colour='red', size=0.2) +
    ggtitle(main) +
    expand_limits(y=0) +                        # Expand y range
    theme_cowplot() +
    theme(panel.grid.major = element_line(colour = "grey",size=0.2,linetype="dashed"), plot.title = element_text(size=20))

  #plot(Expected, Observed, xlim= c(0,1), ylim=c(0,1), pch =19, main = main, ...)
  #lines(x = c(0,1), y = c(0,1), col ="grey", lty = 2)
}

# used by show_param
get_pi1 <- function(phenores){
  gene.causal <- sum(sapply(1:22, function(x) length(phenores$batch[[x]]$idx.cgene)))
  gene.total <- sum(sapply(1:22, function(x) phenores$batch[[x]]$J))

  snp.causal <- sum(sapply(1:22, function(x) length(phenores$batch[[x]]$idx.cSNP)))
  snp.total <- sum(sapply(1:22, function(x) phenores$batch[[x]]$M))

  return(c(gene.causal/gene.total, snp.causal/snp.total))
}

# used by show_param, get mean effect size (abs of effect size and take mean)
get_effect_size <- function(phenores){
  gene_effect <- mean(unlist(sapply(1:22, function(x) abs(phenores$batch[[x]]$e.beta))))
  snp_effect <- mean(unlist(sapply(1:22, function(x) abs(phenores$batch[[x]]$s.theta))))
  return(c(gene_effect, snp_effect))
  }

# return data.frame, rows are each simu run, columns are different parameters.
show_param <- function(phenofs, susieIfs, susieIfs2, thin = 1){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  # overall truth
  truth <- do.call(rbind, lapply(phenofs, function(x) {load(x);
    c(phenores$param$pve.gene.truth, phenores$param$pve.snp.truth,
      get_pi1(phenores), get_effect_size(phenores))}))
  colnames(truth) <- c("PVE.gene_truth", "PVE.SNP_truth", "pi1.gene_truth", "pi1.SNP_truth",
                       "sigma.gene_truth", "sigma.SNP_truth")
  truth <- cbind(truth, truth[, "pi1.gene_truth"]/truth[, "pi1.SNP_truth"])
  colnames(truth)[7] <- "enrich_truth"

  # truth in selected regions for param estimation
  truth.se <- do.call(rbind, lapply(1: length(susieIfs2), function(x) {
    a <- fread(susieIfs2[x], header = T)
    a$ifcausal <- ifelse(a$id %in% cau[[x]], 1, 0)
    c(nrow(a[a$ifcausal == 1 & a$type == "gene" ])/ nrow(a[a$type == "gene"]),
      nrow(a[a$ifcausal == 1 & a$type == "SNP"])/ nrow(a[a$type == "SNP"]))
  }))
  colnames(truth.se) <- c("pi1.gene_se_truth", "pi1.SNP_se_truth")

  # estimated parameters
  est <- do.call(rbind, lapply(susieIfs, function(x) {load(x); group_prior_rec[, ncol(group_prior_rec)]}))
  colnames(est) <- c("pi1.gene_est", "pi1.SNP_est")
  est[, "pi1.SNP_est"] <- est[, "pi1.SNP_est"] * thin
  est <- cbind(est, est[, "pi1.gene_est"]/est[, "pi1.SNP_est"])
  colnames(est)[3] <- "enrich_est"

  est2 <- do.call(rbind, lapply(susieIfs, function(x) {load(x); sqrt(group_prior_var_rec[, ncol(group_prior_var_rec)]/n)}))
  colnames(est2) <- c("sigma.gene_est", "sigma.SNP_est")

  est3 <- cbind( est2[,1]**2 * J * est[,1], est2[,2]**2 * p * est[,2])
  colnames(est3) <- c("PVE.gene_est", "PVE.SNP_est")

  mtx <- cbind(truth, truth.se, est, est2, est3)
  return(mtx)

  # df %>%
  # kable("html", escape = F) %>%
  # kable_styling("striped", full_width = F) %>%
  #  row_spec(c(1:5, 11:15), background = "#FEF3B9") %>%
  # scroll_box(width = "100%", height = "600px", fixed_thead = T)
}

plot_param <- function(mtx){
  n <- nrow(mtx)
  # plot gene pi1
  y1 <- mtx[,"pi1.gene_truth"]
  y2 <- mtx[,"pi1.gene_se_truth"]
  y3 <- mtx[,"pi1.gene_est"]
  plot(jitter(rep(1, n)), y1, col = 1:n, pch = 19, ylab = "gene pi1", xaxt = "n", xlab="", xlim = c(0.8, 3.5), ylim = c(0, max(c(y1, y2, y3)) * 1.05) ,frame.plot=FALSE)
  points(jitter(rep(2, n)), y2, col = 1:n, pch =19)
  points(jitter(rep(3, n)), y3, col = 1:length(y1), pch =19)
  axis(side=1, at=1:2, labels = FALSE, tick = F)
  text(x=1:3, 0, labels = c( "truth", "truth(selected)", "ctwas"), xpd = T, pos =1)
  grid()

  # plot SNP pi1
  y1 <- mtx[,"pi1.SNP_truth"]
  y2 <- mtx[,"pi1.SNP_se_truth"]
  y3 <- mtx[,"pi1.SNP_est"]
  plot(jitter(rep(1, n)), y1, col = 1:n, pch = 19, ylab = "SNP pi1", xaxt = "n", xlab="", xlim = c(0.8, 3.5), ylim = c(0, max(c(y1, y2, y3)) * 1.05) ,frame.plot=FALSE)
  points(jitter(rep(2, n)), y2, col = 1:n, pch =19)
  points(jitter(rep(3, n)), y3, col = 1:length(y1), pch =19)
  axis(side=1, at=1:2, labels = FALSE, tick = F)
  text(x=1:3, 0, labels = c( "truth", "truth(selected)", "ctwas"), xpd = T, pos =1)
  grid()

  #plot enrichment
  y1 <- mtx[,"pi1.gene_truth"]/mtx[,"pi1.SNP_truth"]
  y2 <- mtx[,"pi1.gene_se_truth"]/mtx[,"pi1.SNP_se_truth"]
  y3 <- mtx[,"pi1.gene_est"]/mtx[,"pi1.SNP_est"]
  plot(jitter(rep(1, n)), y1, col = 1:n, pch = 19, ylab = "Enrich", xaxt = "n", xlab="", xlim = c(0.8, 3.5), ylim = c(0, max(c(y1, y2, y3)) * 1.05) ,frame.plot=FALSE)
  points(jitter(rep(2, n)), y2, col = 1:n, pch =19)
  points(jitter(rep(3, n)), y3, col = 1:length(y1), pch =19)
  axis(side=1, at=1:2, labels = FALSE, tick = F)
  text(x=1:3, 0, labels = c( "truth", "truth(selected)", "ctwas"), xpd = T, pos =1)
  grid()
}
