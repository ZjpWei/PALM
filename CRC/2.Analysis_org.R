# =============================================== #
#  CRC 2. Meta-analyses on original data (K401)   #
# =============================================== #

  # Packages ----
  library('tidyverse')
  library('coin')
  library('glmnet')
  library("ANCOMBC")
  library('MicrobiomeStat')
  library('Maaslin2')
  library("PALM")

  # General ----
  rm(list = ls())
  set.seed(2024)

  load("./CRC/Processed_data/data.org.K401.Rdata")
  L <- length(data.rel)
  K <- ncol(data.rel[[1]]$Y)
  rel.abd <- list()
  covariate.interest <- list()
  outcome <- NULL
  feature.table <- NULL
  study <- NULL
  for(l in 1:L){
    study <- c(study, rep(l, length(data.rel[[l]]$X)) )
    feature.table <- rbind(feature.table, data.rel[[l]]$Y)
    outcome <- c(outcome, data.rel[[l]]$X)
    rel.abd[[paste0("S",l)]] <- data.rel[[l]]$Y
    covariate.interest[[paste0("S",l)]] <- matrix(data.rel[[l]]$X, ncol = 1,
                                                  dimnames = list(rownames(data.rel[[l]]$Y), "disease"))
  }

  ## Combine data
  names(outcome) <- rownames(feature.table)
  feature.table = data.frame(t(feature.table))
  meta.data = data.frame(labels = factor(outcome), study = factor(study))
  colnames(feature.table) <- rownames(meta.data)
  feature.ID <- rownames(feature.table)

  ## ANCOM-BC2
  ## function use for ANCOM-BC2
  source("./utility/ancombc.R")

  ## ANCOM-BC2 with sensitivity score filter
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    ANCOMBC2.model <- ancombc.fun(feature.table = feature.table[,meta.data$study == l],
                                  meta = meta.data %>% dplyr::filter(study == l),
                                  formula = "labels",
                                  adjust.method = "fdr",
                                  group = NULL,
                                  subject = NULL,
                                  method = "ancombc2")

    AA.est[ANCOMBC2.model$res$taxon[ANCOMBC2.model$res$passed_ss_labels1],l] <-
      ANCOMBC2.model$res$lfc_labels1[ANCOMBC2.model$res$passed_ss_labels1]
    AA.var[ANCOMBC2.model$res$taxon[ANCOMBC2.model$res$passed_ss_labels1],l] <-
      (ANCOMBC2.model$res$se_labels1[ANCOMBC2.model$res$passed_ss_labels1])^2
  }

  ## Calculate FDR
  tmp.pval <- NULL
  for(k in 1:nrow(AA.est)){
    nonna.id <- !is.na(AA.est[k,])
    if(sum(nonna.id) > 1){
      m <- metafor::rma(yi = AA.est[k,nonna.id], vi = AA.var[k,nonna.id], method = "EE")
      tmp.pval <- c(tmp.pval, m$QEp)
    }else{
      tmp.pval <- c(tmp.pval, NA)
    }
  }

  ## Combined statistics
  est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
  var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
  meta.coef <- est.statics / var.statics
  meta.var <- 1 / var.statics
  q.coef <- (meta.coef)^2 / meta.var
  pval.sin <- 1 - pchisq(q.coef, df = 1)
  qval.sin <- p.adjust(pval.sin, method = "fdr")

  ANCOMBC2.res <- data.frame(features = rownames(AA.est),
                             est = meta.coef,
                             var = meta.var,
                             pval = pval.sin,
                             qval = qval.sin,
                             het.pval = tmp.pval,
                             het.qval = p.adjust(tmp.pval, method = "fdr"))

  ANCOMBC2.model <- list(est = AA.est, var = AA.var)

  ## MaAsLin2
  source("./utility/Maaslin2.R")

  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  ## MaAsLin2 model
  Maaslin2.model <- Maaslin2_meta(feature.table = feature.table,
                                  meta.data = meta.data,
                                  fixed_effects = "labels",
                                  random_effects = NULL,
                                  study = "study",
                                  rma.method = "EE")

  for(l in 1:L){
    AA.est[Maaslin2.model$Maaslin2.lst[[l]]$feature,l] <- Maaslin2.model$Maaslin2.lst[[l]]$coef
    AA.var[Maaslin2.model$Maaslin2.lst[[l]]$feature,l] <- (Maaslin2.model$Maaslin2.lst[[l]]$stderr)^2
  }

  ## Calculate FDR
  tmp.pval <- NULL
  for(k in 1:nrow(AA.est)){
    nonna.id <- !is.na(AA.est[k,])
    if(sum(nonna.id) > 1){
      m <- metafor::rma(yi = AA.est[k,nonna.id], vi = AA.var[k,nonna.id], method = "EE")
      tmp.pval <- c(tmp.pval, m$QEp)
    }else{
      tmp.pval <- c(tmp.pval, NA)
    }
  }

  ## Combined statistics
  est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
  var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
  meta.coef <- est.statics / var.statics
  meta.var <- 1 / var.statics
  q.coef <- (meta.coef)^2 / meta.var
  pval.sin <- 1 - pchisq(q.coef, df = 1)
  qval.sin <- p.adjust(pval.sin, method = "fdr")

  Maaslin2.res <- data.frame(features = rownames(AA.est),
                             est = meta.coef,
                             var = meta.var,
                             pval = pval.sin,
                             qval = qval.sin,
                             het.pval = tmp.pval,
                             het.qval = p.adjust(tmp.pval, method = "fdr"))

  Maaslin2.model <- list(est = AA.est, var = AA.var)

  ## LM-CLR
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    tmp.data <- (data.rel[[l]]$Y)[,colSums(data.rel[[l]]$Y)!=0]
    Z.pool <- log(tmp.data + 0.5)
    clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

    inormal <- apply(clrx_Z, 2, function(d){
      data.d <- data.frame(X = data.rel[[l]]$X, d = d)
      lm.fit <- lm(X ~ d, data = data.d)
      sum.fit <- summary(lm.fit)
      if("d" %in% rownames(sum.fit$coefficients)){
        return(sum.fit$coefficients["d",c("Estimate", "Std. Error")])
      }else{
        return(c(NA, NA))
      }
    })

    AA.est[colnames(inormal),l] <- inormal[1,]
    AA.var[colnames(inormal),l] <- (inormal[2,])^2
  }

  ## Calculate FDR
  tmp.pval <- NULL
  for(k in 1:nrow(AA.est)){
    nonna.id <- !is.na(AA.est[k,])
    if(sum(nonna.id) > 1){
      m <- metafor::rma(yi = AA.est[k,nonna.id], vi = AA.var[k,nonna.id], method = "EE")
      tmp.pval <- c(tmp.pval, m$QEp)
    }else{
      tmp.pval <- c(tmp.pval, NA)
    }
  }

  ## Combined statistics
  est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
  var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
  meta.coef <- est.statics / var.statics
  meta.var <- 1 / var.statics
  q.coef <- (meta.coef)^2 / meta.var
  pval.sin <- 1 - pchisq(q.coef, df = 1)
  qval.sin <- p.adjust(pval.sin, method = "fdr")

  lmclr.res <- data.frame(features = rownames(AA.est),
                          est = meta.coef,
                          var = meta.var,
                          pval = pval.sin,
                          qval = qval.sin,
                          het.pval = tmp.pval,
                          het.qval = p.adjust(tmp.pval, method = "fdr"))

  lmclr.model <- list(est = AA.est, var = AA.var)

  ## LinDA
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    Linda.model <- MicrobiomeStat::linda(feature.dat = t(data.rel[[l]]$Y),
                                         meta.dat = meta.data %>% dplyr::filter(study == l),
                                         formula = paste0('~labels'),
                                         feature.dat.type = "count",
                                         prev.filter = 0,
                                         adaptive = TRUE,
                                         alpha = 0.05)

    AA.est[rownames(Linda.model$output$labels1),l] <- Linda.model$output$labels1$log2FoldChange
    AA.var[rownames(Linda.model$output$labels1),l] <- (Linda.model$output$labels1$lfcSE)^2
  }

  ## Calculate FDR
  tmp.pval <- NULL
  for(k in 1:nrow(AA.est)){
    nonna.id <- !is.na(AA.est[k,])
    if(sum(nonna.id) > 1){
      m <- metafor::rma(yi = AA.est[k,nonna.id], vi = AA.var[k,nonna.id], method = "EE")
      tmp.pval <- c(tmp.pval, m$QEp)
    }else{
      tmp.pval <- c(tmp.pval, NA)
    }
  }

  ## Combined statistics
  est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
  var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
  meta.coef <- est.statics / var.statics
  meta.var <- 1 / var.statics
  q.coef <- (meta.coef)^2 / meta.var
  pval.sin <- 1 - pchisq(q.coef, df = 1)
  qval.sin <- p.adjust(pval.sin, method = "fdr")

  Linda.res <- data.frame(features = rownames(AA.est),
                          est = meta.coef,
                          var = meta.var,
                          pval = pval.sin,
                          qval = qval.sin,
                          het.pval = tmp.pval,
                          het.qval = p.adjust(tmp.pval, method = "fdr"))

  Linda.model <- list(est = AA.est, var = AA.var)

  ## N form whole data
  load("./CRC/Processed_data/data.org.K849.Rdata")
  N <- list()
  for(l in 1:L){
    N[[paste0("S",l)]] <- rowSums(data.rel[[l]]$Y)
  }

  ## Meta-analysis
  null_obj_RA <- palm.null.model(rel.abd = rel.abd, N = N, prev.filter = 0)

  summary.score.RA <- palm.get.summary(null.obj = null_obj_RA,
                                       covariate.interest = covariate.interest)

  ## original PALM corrected
  PALM.model <- palm.test(summary.stats = summary.score.RA, p.adjust.method = "fdr")

  ## Calculate FDR
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    AA.est[PALM.model$disease$palm_fits[[l]]$feature,l] <- PALM.model$disease$palm_fits[[l]]$coef
    AA.var[PALM.model$disease$palm_fits[[l]]$feature,l] <- (PALM.model$disease$palm_fits[[l]]$stderr)^2
  }

  ## Calculate FDR
  tmp.pval <- NULL
  for(k in 1:nrow(AA.est)){
    nonna.id <- !is.na(AA.est[k,])
    m <- metafor::rma(yi = AA.est[k,nonna.id], vi = AA.var[k,nonna.id], method = "EE")
    tmp.pval <- c(tmp.pval, m$QEp)
  }
  PALM.res <- PALM.model$disease$meta_fits
  PALM.res$het.pval = tmp.pval
  PALM.res$het.qval = p.adjust(tmp.pval, method = "fdr")
  PALM.model <- list(est = AA.est, var = AA.var)

  ## Save output
  save(ANCOMBC2.model, lmclr.model, Maaslin2.model, Linda.model, PALM.model,
       ANCOMBC2.res, lmclr.res, Maaslin2.res, Linda.res, PALM.res,
       file = "./CRC/Output/CRC_output.Rdata")

