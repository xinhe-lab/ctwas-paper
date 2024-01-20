
#' Simulate phenotype
#' @description simulate Y under the sparse model (spike and slab prior).
#' will standardize gene expression and SNP genotype when simulate.
#' @param pve.expr scalar [0,1], not used when if `sigma_beta` is given
#' @param pve.snp scalar [0,1], not used if `sigma_theta` is given
#' @param pi_beta scalar [0,1]
#' @param pi_theta scalar [0,1]
#' @param sigma_beta default NULL, if not given, will calculate from pve.
#' @param sigma_theta default NULL, if not given, will calculate from pve.
#' @param prior_dist the prior distribution for effect size of causal variants, `mixnormal` prior is a mixture of normal distribution
#' @param pi_k the mixture proportion for `mixnormal` prior. sum(pi_k) = 1
#' @param sa2_kr relative sa2 for `mixnormal` prior
#'       sa2_k/sigma_beta^2 or sa2_k/sigma_theta^2. sum(pi_k * sa2_kr) = 1
#' @return a list, Y.g is phenotype determined by genetic, will need to add
#'  error term (pve calculation is based on N(0,1)).
#' requires tibble, pgenlibr, scaleRcpp
simulate_phenotype<- function(pgenfs,
                              exprfs,
                              pve.expr = 0.1,
                              pve.snp = 0.1,
                              pi_beta = 0.01,
                              pi_theta = 0.0001,
                              sigma_beta = NULL,
                              sigma_theta = NULL,
                              prior_dist_causal = c("normal", "mixnormal"),
                              pi_k = NULL,
                              sa2_kr = NULL
) {
  pvarfs <- sapply(pgenfs, ctwas:::prep_pvar, outputdir = outputdir)
  exprvarfs <- sapply(exprfs, ctwas:::prep_exprvar)

  pgen1 <- ctwas:::prep_pgen(pgenf = pgenfs[1], pvarfs[1])
  N <- pgenlibr::GetRawSampleCt(pgen1)

  M.b <- sapply(pvarfs, function(x) nrow(ctwas:::read_pvar(x)))
  J.b <- sapply(exprvarfs, function(x) nrow(ctwas:::read_exprvar(x)))

  M <- sum(M.b)
  J <- sum(J.b)

  if (is.null(sigma_beta)){
    expr.meanvar <- 1 # always scale before simulate
    J.c <- round(J * pi_beta)
    sigma_beta <- sqrt(pve.expr / (J.c * expr.meanvar * (1 - pve.snp - pve.expr)))
    if (is.infinite(sigma_beta)) sigma_beta <- 0
  }

  if (is.null(sigma_theta)){
    M.c <- round(M * pi_theta)
    sigma_theta <- sqrt(pve.snp / (M.c * (1 - pve.snp - pve.expr)))
  }

  phenores <- list("batch" = list())
  for ( b in 1:length(pgenfs)){

    idx.cgene <- which(rbinom(J.b[b], 1, pi_beta) == 1)
    idx.cSNP <- which(rbinom(M.b[b], 1, pi_theta) == 1)

    X.g <- ctwas:::read_expr(exprfs[b], variantidx = idx.cgene)
    if (is.null(X.g)){
      X.g <- matrix(, nrow= N, ncol =0)
    }

    X.g <- scaleRcpp(X.g) # always scale expr

    pgen <- ctwas:::prep_pgen(pgenf = pgenfs[b], pvarfs[b])

    X.s <- ctwas:::read_pgen(pgen, variantidx = idx.cSNP)
    X.s <- scaleRcpp(X.s) # always scale genotype

    prior_dist_causal <- match.arg(prior_dist_causal)

    if (prior_dist_causal == "normal"){
      e.beta <- rnorm(length(idx.cgene), mean = 0, sd = sigma_beta)
      s.theta <- rnorm(length(idx.cSNP), mean = 0, sd = sigma_theta)
    } else if (prior_dist_causal == "mixnormal"){

      if (is.null(pi_k)){
        # mixture of 4 normal, equal proportion
        nmix <- 4
        pi_k <- rep(1/nmix, nmix)
      }

      if (sum(pi_k) != 1){
        stop("Improper prior distribution. Stopped.")
      }

      if (is.null(sa2_kr)){
        # sa2 difference is 2 fold
        ra <- 2**(1: length(pi_k))
        sa2_kr <-  ra / (pi_k * sum(ra))
      }

      if (sum(pi_k * sa2_kr) != 1){
        stop("Improper prior distribution. Stopped.")
      }

      e.beta <- mixtools::rnormmix(length(idx.cgene), lambda = pi_k,
                          mu = rep(0, length(pi_k)),
                          sigma = sigma_beta * sqrt(sa2_kr))

      s.theta <- mixtools::rnormmix(length(idx.cSNP), lambda = pi_k,
                           mu = rep(0, length(pi_k)),
                           sigma = sigma_theta * sqrt(sa2_kr))
    }

    Y.g <- X.g %*% e.beta + X.s %*% s.theta

    var.gene <- var(X.g %*% e.beta)
    var.snp <- var(X.s %*% s.theta)

    id.cgene <- ctwas:::read_exprvar(exprvarfs[b])[idx.cgene,][["id"]]

    id.cSNP <-  ctwas:::read_pvar(pvarfs[b])[idx.cSNP, ][["id"]]

    phenores[["batch"]][[b]] <- tibble::lst(Y.g, s.theta, e.beta,
                                             sigma_theta, sigma_beta,
                                             idx.cSNP, idx.cgene,
                                             J = J.b[b], M = M.b[b],
                                             N, var.gene, var.snp,
                                             id.cSNP, id.cgene, pi_k,
                                             sa2_kr)
  }

  Y <- matrix(rowSums(do.call(cbind, lapply(phenores[["batch"]], '[[', "Y.g"))) + rnorm(N), ncol = 1)
  phenores$Y <- Y

  var.y <- var(Y)

  phenores$param$pve.snp.truth <-
    sum(unlist(lapply(phenores[["batch"]], '[[', "var.snp")))/var.y
  phenores$param$pve.gene.truth <-
    sum(unlist(lapply(phenores[["batch"]], '[[', "var.gene")))/var.y

  return(phenores)

}

simulate_phenotype_multi <- function(pgenfs,
                                     exprfs,
                                     weight_1,
                                     pve_expr_1 = 0.1,
                                     pi_beta_1 = 0.01,
                                     weight_2,
                                     pve_expr_2 = 0.1,
                                     pi_beta_2 = 0.01,
                                     pve_snp = 0.1,
                                     pi_theta = 0.0001,
                                     sigma_beta = NULL,
                                     sigma_theta = NULL){

  weight_1 <- rev(unlist(strsplit(tools::file_path_sans_ext(weight_1), "/")))[1]
  weight_2 <- rev(unlist(strsplit(tools::file_path_sans_ext(weight_2), "/")))[1]
  weight <- c(weight_1, weight_2)

  pi_beta <- c(pi_beta_1, pi_beta_2)
  pve_expr <- c(pve_expr_1, pve_expr_2)

  pvarfs <- sapply(pgenfs, ctwas:::prep_pvar, outputdir = outputdir)
  exprvarfs <- sapply(exprfs, ctwas:::prep_exprvar)

  pgen1 <- ctwas:::prep_pgen(pgenf = pgenfs[1], pvarfs[1])
  N <- pgenlibr::GetRawSampleCt(pgen1)

  M.b <- sapply(pvarfs, function(x) nrow(ctwas:::read_pvar(x)))
  J.b <- t(sapply(exprvarfs, function(x){df <- ctwas:::read_exprvar(x); df$weight <- sapply(df$id, function(y){unlist(strsplit(y, "[|]"))[2]}); table(df$weight)[weight]}))

  M <- sum(M.b)
  J <- apply(J.b, 2, sum)

  if (is.null(sigma_beta)){
    expr_meanvar <- 1 # always scale before simulate
    J.c <- round(J * pi_beta)
    sigma_beta <- sqrt(pve_expr / (J.c * expr_meanvar * (1 - pve_snp - sum(pve_expr))))
    sigma_beta[is.infinite(sigma_beta) | is.nan(sigma_beta)] <- 0
  }

  if (is.null(sigma_theta)){
    M.c <- round(M * pi_theta)
    sigma_theta <- sqrt(pve_snp / (M.c * (1 - pve_snp - sum(pve_expr))))
  }

  phenores <- list("batch" = list())

  for (b in 1:length(pgenfs)){
    idx.cSNP <- which(rbinom(M.b[b], 1, pi_theta) == 1)

    offset <- c(0,rev(rev(J.b[b,])[-1]))
    idx.cgene <- lapply(1:ncol(J.b), function(x){which(rbinom(J.b[b,x], 1, pi_beta[x])==1)+offset[x]})
    n_cgene <- sapply(idx.cgene, length)
    idx.cgene <- unlist(idx.cgene)
    names(idx.cgene) <- NULL

    X.g <- ctwas:::read_expr(exprfs[b], variantidx = idx.cgene)
    if (is.null(X.g)){
      X.g <- matrix(NA, nrow=N, ncol=0)
    }

    X.g <- scale(X.g) # always scale expr

    # change invariant SNPs from NaN to 0
    X.g[,apply(X.g, 2, function(x){any(is.nan(x))})] <- 0

    pgen <- ctwas:::prep_pgen(pgenf = pgenfs[b], pvarfs[b])

    X.s <- ctwas:::read_pgen(pgen, variantidx = idx.cSNP)
    X.s <- scale(X.s) # always scale genotype

    s.theta <- rnorm(length(idx.cSNP), mean = 0, sd = sigma_theta)
    e.beta <- unlist(lapply(1:length(n_cgene), function(x){rnorm(n_cgene[x], mean=0, sd=sigma_beta[x])}))

    Y.g <- X.g %*% e.beta + X.s %*% s.theta

    var.gene <- var(X.g %*% e.beta)
    var.snp <- var(X.s %*% s.theta)

    id.cgene <- ctwas:::read_exprvar(exprvarfs[b])[idx.cgene,][["id"]]
    id.cSNP <-  ctwas:::read_pvar(pvarfs[b])[idx.cSNP, ][["id"]]

    phenores[["batch"]][[b]] <- tibble::lst(Y.g, s.theta, e.beta,
                                            sigma_theta, sigma_beta,
                                            idx.cSNP, idx.cgene,
                                            J = J.b[b,], M = M.b[b],
                                            N, var.gene, var.snp,
                                            id.cSNP, id.cgene)
  }

  Y <- matrix(rowSums(do.call(cbind, lapply(phenores[["batch"]], '[[', "Y.g"))) + rnorm(N), ncol = 1)
  phenores$Y <- Y

  var.y <- var(Y)

  phenores$param$pve.snp.truth <-
    sum(unlist(lapply(phenores[["batch"]], '[[', "var.snp")))/var.y
  phenores$param$pve.gene.truth <-
    sum(unlist(lapply(phenores[["batch"]], '[[', "var.gene")))/var.y

  return(phenores)
}
