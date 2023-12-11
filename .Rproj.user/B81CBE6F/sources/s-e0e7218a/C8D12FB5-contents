#library(mr.ash.alpha)
library(data.table)
suppressMessages({library(plotly)})

#source("~/causalTWAS/causal-TWAS/code/before_package/fit_mr.ash.R")

get_files <- function(tag, tag2){
  par <- paste0(outputdir, tag, "-mr.ash2s.", tag2, ".param.txt")
  gmrash <- paste0(outputdir, tag, "-mr.ash2s.", tag2, ".expr.txt")
  ggwas <- paste0(outputdir, tag, ".exprgwas.txt.gz")
  smrash <- paste0(outputdir, tag, "-mr.ash2s.", tag2, ".snp.txt")
  sgwas <- paste0(outputdir, tag, ".snpgwas.txt.gz")
  gsusie <- paste0(susiedir, tag, ".", tag2, ".L3.susieres.expr.txt")
  ssusie <- paste0(susiedir, tag, ".", tag2, ".L3.susieres.snp.txt")
  rsq <- paste0(outputdir, tag, ".LD.txt")
  return(tibble::lst(par, gmrash, ggwas, smrash, sgwas, gsusie, ssusie, rsq))
}

plot_pip <- function(dt){

  pal <- c("salmon", "darkgreen")
  pal <- setNames(pal, c("Causal", "Non causal"))

  if ( "max_r2" %in% colnames(dt)) {
    fig <- plot_ly(dt, x = ~p0, y = ~ PIP, type = 'scatter',
                   mode = 'markers', color = ~ ifcausal, colors = pal,
                   text = ~ paste("Name: ", name,
                                  "\nmax R^2: ", round(max_r2, digits = 2)))
    #,"\nMax pair: ", max_pair))
  } else {
    fig <- plot_ly(dt, x = ~p0, y = ~ PIP, type = 'scatter', mode = 'markers', color = ~ ifcausal, colors = pal, text = ~paste("Name: ", name))
  }

  fig <- fig %>% layout(yaxis = list(title = "PIP"))
  return(fig)
}

plot_beta <- function(dt){
  fig <- plot_ly(dt, x = ~p0, y = ~ TrueBeta, type = 'scatter', mode = 'markers', name ="Truth", text = ~paste("Name: ", name))
  fig <- fig %>% add_trace(y = ~ EstBeta, name = 'Estimated',mode = 'markers')

  fig <- fig %>% layout(yaxis = list(title = "beta"))
  return(fig)
}

plot_lbf <- function(dt){
  pal <- c("salmon", "darkgreen")
  pal <- setNames(pal, c("Causal", "Non causal"))

  fig <- plot_ly(dt, x = ~p0, y = ~ lbf , type = 'scatter', mode = 'markers', color = ~ ifcausal, colors = pal, text = ~paste("Name: ", name))

  fig <- fig %>% layout(yaxis = list(title = "log10(BF)"))
  return(fig)
}


plot_p <- function(dt){
  pal <- c("salmon", "darkgreen")
  pal <- setNames(pal, c("Causal", "Non causal"))

  fig <- plot_ly(dt, x = ~p0, y = ~ PVALUE , type = 'scatter', mode = 'markers', color = ~ ifcausal, colors = pal, text = ~paste("Name: ", name))

  fig <- fig %>% layout(yaxis = list(title = "p value"))
  return(fig)
}

plot_susie_pip <- function(dt) {
  pal <- c("salmon", "darkgreen")
  pal <- setNames(pal, c("Causal", "Non causal"))

  fig <- plot_ly(dt, x = ~ p0, y = ~ susiePIP , type = 'scatter', mode = 'markers', color = ~ ifcausal, colors = pal, text = ~paste("Name: ", name))

  fig <- fig %>% layout(yaxis = list(title = "SUSIE PIP"))
}

plot_coloc_pp4 <- function(dt) {
  pal <- c("salmon", "darkgreen")
  pal <- setNames(pal, c("Causal", "Non causal"))

  fig <- plot_ly(dt, x = ~ p0, y = ~ COLOC.PP4 , type = 'scatter', mode = 'markers', color = ~ ifcausal, colors = pal, text = ~paste("Name: ", name))

  fig <- fig %>% layout(yaxis = list(title = "coloc PP4"))
}


plot_all <- function(mrashres.gene, gwasres.gene, mrashres.snp, gwasres.snp, susieres.gene, susieres.snp, ldres = NULL, chrom=NULL){
  a <- read.table(mrashres.gene , header =T)
  a[a$ifcausal == 1, "ifcausal"] = "Causal"
  a[a$ifcausal == 0, "ifcausal"] = "Non causal"

  if (!is.null(ldres)){
    ld <- read.table(ldres, header =T)
    a <- merge(a, ld, by = "name", all.x = T)}


  su <- read.table(susieres.gene, header =T)
  colnames(su)[colnames(su) == "pip"] <-  "susiePIP"
  a <- merge(a, su, by = "name", all.x = T)

  p <- fread(gwasres.gene, header =T)
  p[,"PVALUE"] <- -log10(p[,"PVALUE"])
  colnames(p)[colnames(p) == "MARKER_ID"] <-  "name"
  dt.gene <- merge(a, p, by = "name", all.x = T)

  a <- read.table(mrashres.snp , header =T)
  a[a$ifcausal == 1, "ifcausal"] = "Causal"
  a[a$ifcausal == 0, "ifcausal"] = "Non causal"

  su <- read.table(susieres.snp, header =T)
  colnames(su)[colnames(su) == "pip"] <-  "susiePIP"
  a <- merge(a, su, by = "name", all.x = T)

  p <- fread(gwasres.snp, header =T)
  p[,"PVALUE"] <- -log10(p[,"PVALUE"])
  colnames(p)[colnames(p) == "MARKER_ID"] <-  "name"
  dt.snp <- merge(a, p, by = "name", all.x = T)

  if (!is.null(chrom)){
    dt.gene <- dt.gene[dt.gene$chr.x == chrom,]
    dt.snp <- dt.snp[dt.snp$chr.x == chrom,]
  }


  fig1 <- plot_pip(dt.gene)
  fig2 <- plot_beta(dt.gene)
  fig3 <- plot_p(dt.gene)
  fig4 <- plot_susie_pip(dt.gene)
  fig5 <- plot_pip(dt.snp)
  fig6 <- plot_beta(dt.snp)
  fig7 <- plot_p(dt.snp)
  fig8 <- plot_susie_pip(dt.snp) %>% layout( plot_bgcolor='#EFF6FC')

  fig <- subplot(fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, nrows=8, shareX = T, titleY = T)

  fig <- fig %>% layout(autosize = F, width = 600, height = 1100, margin = 0.02,
                        xaxis = list(title = "Genomic position"))

  fig <- fig %>% toWebGL()
  fig
}
