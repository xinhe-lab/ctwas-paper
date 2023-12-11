
library(plyr)

# PIP calibration plot
caliPIP_plot <- function(phenofs, pipfs, main = "PIP Calibration"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(pipfs)) {
    res <- fread(pipfs[i], header = T)
    res <- data.frame(res[res$type  =="gene", ])
    res$runtag <- i
    res <- res[complete.cases(res),]
    res$ifcausal <- ifelse(res$id %in% cau[[i]], 1, 0)
    df <- rbind(df, res)
  }

  fig <- cp_plot(df$susie_pip, df$ifcausal, df$runtag, mode ="PIP", main = main)
  return(fig)
}

# Power plot, plot the number of causal genes in PIP bins.
ncausal_plot <- function(phenofs, pipfs, main = "PIP"){

  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})

  df <- NULL
  for (i in 1:length(pipfs)) {
    res <- fread(pipfs[i], header = T)
    res <- data.frame(res[res$type  =="gene", ])
    res$ifcausal <- ifelse(res$id %in% cau[[i]], 1, 0)
    res$runtag <- i
    res <- res[complete.cases(res),]
    df <- rbind(df, res)
  }

  fig <- nca_plot(df$susie_pip, df$ifcausal, df$runtag, mode ="PIP", xmin = 0.5, main = main)
  return(fig)
}

scatter_plot_PIP_p <- function(phenofs, pipfs, gwasfs, main ="PIP-p"){
  cau <- lapply(phenofs, function(x) {load(x);get_causal_id(phenores)})
  df <- NULL
  for (i in 1:length(pipfs)) {
    pipres <- fread(pipfs[i], header = T)
    pipres <- data.frame(pipres[pipres$type == "gene", ])
    pipres$runtag <- i
    pipres$ifcausal <- ifelse(pipres$id %in% cau[[i]], 1, 0)
    gwasres <- read.table(gwasfs[i], header = T)
    res <- merge(gwasres, pipres, by = "id", all = T)
    res <- res[complete.cases(res),]
    df <- rbind(df, res)
  }

  colnames(df)[colnames(df)== "PVALUE"] <- "TWAS.p"
  df[,"TWAS.p"] <- -log10(df[, "TWAS.p"])

  df$ifcausal <- mapvalues(df$ifcausal, from=c(0,1),to=c("darkgreen", "salmon"))
  plot(df$TWAS.p, df$susie_pip, bg = df$ifcausal, main = main, xlab = "-log10(TWAS p value)", ylab = "PIP", pch =21, bty='n')
  grid()

  # df$ifcausal <- mapvalues(df$ifcausal, from=c(0,1), to=c("Non causal", "Causal"))
  # fig <- plot_ly(data = df, x = ~ TWAS.p, y = ~ susie_pip, color = ~ ifcausal,
  #                colors = c( "salmon", "darkgreen"), type ="scatter", text = ~ paste("Name: ", paste0(runtag,":", id),
  #                                                                  "\nChr: ", chrom.x,  "\nPos:", pos))
  # fig
  return(df)
}
