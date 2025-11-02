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
  library("LOCOM")
  library("DESeq2")

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
                             stderr = sqrt(meta.var),
                             pval = pval.sin,
                             qval = qval.sin,
                             pval.het = tmp.pval,
                             qval.het = p.adjust(tmp.pval, method = "fdr"))

  ANCOMBC2.model <- list(est = AA.est, var = AA.var)

  ## DESeq2
  ## DESeq2 with sensitivity score filter
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    objs <- DESeqDataSetFromMatrix(countData = feature.table[,meta.data$study == l],
                                   colData = meta.data %>% dplyr::filter(study == l) %>%
                                     dplyr::transmute(labels = labels),
                                   design = ~ labels)

    DESeq2.model <- DESeq(object = objs, sfType = "poscounts")

    DESeq2.test <- results(DESeq2.model)

    AA.est[rownames(DESeq2.test),l] <- DESeq2.test$log2FoldChange
    AA.var[rownames(DESeq2.test),l] <- (DESeq2.test$lfcSE)^2
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

  DESeq2.res <- data.frame(features = rownames(AA.est),
                           est = meta.coef,
                           stderr = sqrt(meta.var),
                           pval = pval.sin,
                           qval = qval.sin,
                           pval.het = tmp.pval,
                           qval.het = p.adjust(tmp.pval, method = "fdr"))

  DESeq2.model <- list(est = AA.est, var = AA.var)

  ## MaAsLin2
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    data.tmp <- feature.table[,meta.data$study == l]
    data.tmp <- data.tmp[rowSums(data.tmp) > 0,]
    taxa.id <- rownames(data.tmp)
    rownames(data.tmp) <- paste0("Tax", 1:length(taxa.id))
    res <- Maaslin2::Maaslin2(input_data = data.tmp,
                              input_metadata = meta.data[meta.data$study == l,],
                              output = "./",
                              min_abundance = 0,
                              min_prevalence = 0,
                              normalization = "TSS",
                              transform = 'AST',
                              analysis_method = "LM",
                              max_significance = 1,
                              random_effects = NULL,
                              fixed_effects = "labels",
                              standardize = FALSE,
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE)

    res.tmp <- res$result[res$results$metadata == "labels",]
    res.tmp <- res.tmp[order(as.numeric(sapply(strsplit(res.tmp$feature, split = "Tax"), "[[", 2, simplify = TRUE))),]
    res.tmp$feature <- taxa.id

    AA.est[res.tmp$feature,l] <- res.tmp$coef
    AA.var[res.tmp$feature,l] <- (res.tmp$stderr)^2
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
                             stderr = sqrt(meta.var),
                             pval = pval.sin,
                             qval = qval.sin,
                             pval.het = tmp.pval,
                             qval.het = p.adjust(tmp.pval, method = "fdr"))

  Maaslin2.model <- list(est = AA.est, var = AA.var)

  ## LM-CLR
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    data.tmp <- feature.table[,meta.data$study == l]
    data.tmp <- data.tmp[rowSums(data.tmp) > 0,]
    Z.pool <- log(t(data.tmp) + 0.5)
    clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

    inormal <- apply(clrx_Z, 2, function(d){
      data.d <- data.frame(labels = meta.data$labels[meta.data$study == l], d = d)
      lm.fit <- lm(d ~ labels, data = data.d)
      sum.fit <- summary(lm.fit)
      if("labels1" %in% rownames(sum.fit$coefficients)){
        return(sum.fit$coefficients["labels1",c("Estimate", "Std. Error")])
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
                          stderr = sqrt(meta.var),
                          pval = pval.sin,
                          qval = qval.sin,
                          pval.het = tmp.pval,
                          qval.het = p.adjust(tmp.pval, method = "fdr"))

  lmclr.model <- list(est = AA.est, var = AA.var)

  ## LinDA
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    data.tmp <- feature.table[,meta.data$study == l]
    data.tmp <- data.tmp[rowSums(data.tmp) > 0,]

    Linda.model <- MicrobiomeStat::linda(feature.dat = data.tmp,
                                         meta.dat = meta.data %>% dplyr::filter(study == l),
                                         formula = '~labels',
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
                          stderr = sqrt(meta.var),
                          pval = pval.sin,
                          qval = qval.sin,
                          pval.het = tmp.pval,
                          qval.het = p.adjust(tmp.pval, method = "fdr"))

  Linda.model <- list(est = AA.est, var = AA.var)

  ## DESeq2
  ## DESeq2 with sensitivity score filter
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  ## DESeq2
  for(l in 1:L){
    objs <- DESeqDataSetFromMatrix(countData = feature.table[,meta.data$study == l],
                                   colData = meta.data %>% dplyr::filter(study == l) %>%
                                     dplyr::transmute(labels = labels),
                                   design = ~ labels)

    DESeq2.model <- DESeq(object = objs, sfType = "poscounts", minReplicatesForReplace = Inf)

    DESeq2.test <- results(DESeq2.model, name = "labels_1_vs_0",
                           independentFiltering = FALSE, cooksCutoff = FALSE)

    AA.est[rownames(DESeq2.test),l] <- DESeq2.test$log2FoldChange
    AA.var[rownames(DESeq2.test),l] <- (DESeq2.test$lfcSE)^2
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

  DESeq2.res <- data.frame(features = rownames(AA.est),
                           est = meta.coef,
                           stderr = sqrt(meta.var),
                           pval = pval.sin,
                           qval = qval.sin,
                           pval.het = tmp.pval,
                           qval.het = p.adjust(tmp.pval, method = "fdr"))

  DESeq2.model <- list(est = AA.est, var = AA.var)

  ## N form whole data
  load("./CRC/Processed_data/data.org.K849.Rdata")
  N <- list()
  for(l in 1:L){
    N[[paste0("S",l)]] <- rowSums(data.rel[[l]]$Y)
  }

  ## Meta-analysis
  null_obj_RA <- PALM::palm.null.model(rel.abd = rel.abd, depth = N, prev.filter = 0)

  summary.score.RA <- PALM::palm.get.summary(null.obj = null_obj_RA,
                                             covariate.interest = covariate.interest)

  ## original PALM corrected
  PALM.model <- palm.meta.summary(summary.stats = summary.score.RA, p.adjust.method = "fdr")

  ## Calculate FDR
  AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))
  AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                   dimnames = list(feature.ID, as.character(1:L)))

  for(l in 1:L){
    AA.est[PALM.model$disease$feature,l] <- PALM.model$disease[[paste0("S", l,"_effect")]]
    AA.var[PALM.model$disease$feature,l] <- (PALM.model$disease[[paste0("S", l,"_stderr")]])^2
  }

  ## Calculate FDR
  tmp.pval <- NULL
  for(k in 1:nrow(AA.est)){
    nonna.id <- !is.na(AA.est[k,])
    m <- metafor::rma(yi = AA.est[k,nonna.id], vi = AA.var[k,nonna.id], method = "EE")
    tmp.pval <- c(tmp.pval, m$QEp)
  }
  PALM.res <- PALM.model$disease
  PALM.model <- list(est = AA.est, var = AA.var)

  ## Save output
  save(ANCOMBC2.model, lmclr.model, Maaslin2.model, Linda.model, PALM.model, DESeq2.model,
       ANCOMBC2.res, lmclr.res, Maaslin2.res, Linda.res, PALM.res, DESeq2.res,
       file = "./CRC/Output/CRC_output.Rdata")

