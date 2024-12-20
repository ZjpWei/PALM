# =============================================== #
# (3) Analyze microbiome-metabolome data (other)  #
# =============================================== #

  # Packages ----
  library("coin")
  library("purrr")
  library('tidyverse')
  library('glmnet')
  library("ANCOMBC")
  library("Maaslin2")

  rm(list = ls())

  ## Load data
  load("./Metabolites/Processed_data/processed_data.Rdata")
  load("./Metabolites/Output/PALM_null.Rdata")
  source("./utility/ancombc.R")

  cluster <- lapply(metadata, function(d){
    d <- d %>% dplyr::transmute(Sample, Subject) %>% data.frame()
    rownames(d) <- NULL
    d <- d %>% tibble::column_to_rownames("Sample")
    return(d)
  })
  Metabolite.ids <- unique(common.pairs$Compound)
  Genera.id <- unique(common.pairs$Taxon)

  covariates.disease <- lapply(metadata[datasets], function(d){
    case.names <- c("Healthy", "Control", "Normal", "H", "0", "nonIBD", "NA")
    return(d %>% tibble::rownames_to_column(var = "IDS") %>%
             dplyr::transmute(Sample, Study.Group =  case_when(Study.Group %in% case.names ~ "0",
                                                               .default = "1"))
    )
  })

  # Choose original compounds
  datasets.cmpd <- c()
  covariates_adjust_lst <- list()
  for(d in datasets){
    if(length(data.for.lm[[d]]$cmpd.name) != length(unique(data.for.lm[[d]]$cmpd.name))){
      otu_data <- matrix(0, nrow = nrow(data.for.lm[[d]]$rel.abd), ncol = length(Genera.id),
                         dimnames = list(rownames(data.for.lm[[d]]$rel.abd), Genera.id))
      datasets.cmpd <- c(datasets.cmpd, d)
      cmpd.name <- names(which(table(data.for.lm[[d]]$cmpd.name) > 1))
      cmpd.dat <- data.for.lm[[d]]$rel.cmpd[,data.for.lm[[d]]$cmpd.name %in% cmpd.name]
      otu_data[,colnames(data.for.lm[[d]]$rel.abd)] <- data.for.lm[[d]]$rel.abd %>% as.matrix()
      sample.kep <- intersect(rownames(otu_data), rownames(cmpd.dat))
    }
  }

  # Process data
  for(d in datasets){
    if(d %in% datasets.cmpd){
      load(paste0("./Metabolites/Processed_data/Metabolite_select/Compound_", d, ".Rdata"))
      remove.cmpd <- setdiff(cmpd_sty_scan$Orig.Compound, cmpd_sty_scan %>%
                               dplyr::filter(selected_num == max(selected_num), .by = Compound) %>%
                               dplyr::filter(row_number() == 1, .by = Compound) %>%
                               pull(Orig.Compound))
      keep.cmpd <- !(colnames(data.for.lm[[d]]$rel.cmpd) %in% remove.cmpd)
      cmpd.dat <- data.for.lm[[d]]$rel.cmpd[,keep.cmpd]
      data.for.lm[[d]]$cmpd.name <- data.for.lm[[d]]$cmpd.name[keep.cmpd]
      colnames(cmpd.dat) <- data.for.lm[[d]]$cmpd.name
      data.for.lm[[d]]$rel.cmpd <- cmpd.dat
    }else{
      colnames(data.for.lm[[d]]$rel.cmpd) <- data.for.lm[[d]]$cmpd.name
    }
  }

  ## function
  meta_analysis <- function(AA.est, AA.var, rm.method = "EE"){
    ## Combined statistics
    est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
    var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
    meta.coef <- est.statics / var.statics
    meta.var <- 1 / var.statics
    q.coef <- (meta.coef)^2 / meta.var

    ## Calculate FDR
    pval.sin <- 1 - pchisq(q.coef, df = 1)
    tmp.pval <- NULL
    for(k in 1:nrow(AA.est)){
      nonna.id <- !is.na(AA.est[k,])
      m <- try(metafor::rma(yi = AA.est[k,nonna.id], vi = AA.var[k,nonna.id],
                            control=list(maxiter=2000), method = rm.method), silent = TRUE)
      if(class(m)[1] != "try-error"){
        tmp.pval <- c(tmp.pval, m$QEp)
      }else{
        tmp.pval <- c(tmp.pval, NA)
      }
    }

    combined.model <- data.frame(feature = Genera.id, pval = pval.sin, pval.het = tmp.pval,
                                 coef = meta.coef, sd = sqrt(meta.var))

    return(combined.model)
  }

  ## MTBL analysis
  for(cmpds.id in Metabolite.ids){
    ## ANCOM-BC2
    ANCOMBC2.meta.pval <- data.frame(feature = Genera.id)
    ANCOMBC2.meta.pval.het <- data.frame(feature = Genera.id)
    ANCOMBC2.meta.coef <- data.frame(feature = Genera.id)
    ANCOMBC2.meta.sd <- data.frame(feature = Genera.id)

    ## MaAsLin2
    Maaslin2.meta.pval <- data.frame(feature = Genera.id)
    Maaslin2.meta.pval.het <- data.frame(feature = Genera.id)
    Maaslin2.meta.coef <- data.frame(feature = Genera.id)
    Maaslin2.meta.sd <- data.frame(feature = Genera.id)

    ## LM-CLR
    LMCLR.meta.pval <- data.frame(feature = Genera.id)
    LMCLR.meta.pval.het <- data.frame(feature = Genera.id)
    LMCLR.meta.coef <- data.frame(feature = Genera.id)
    LMCLR.meta.sd <- data.frame(feature = Genera.id)

    ## LinDA
    Linda.meta.pval <- data.frame(feature = Genera.id)
    Linda.meta.pval.het <- data.frame(feature = Genera.id)
    Linda.meta.coef <- data.frame(feature = Genera.id)
    Linda.meta.sd <- data.frame(feature = Genera.id)

    ## ANCOM-BC2
    ANCOMBC2.est <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                           dimnames = list(Genera.id, datasets))
    ANCOMBC2.var <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                           dimnames = list(Genera.id, datasets))

    ## MaAsLin2
    Maaslin2.est <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                           dimnames = list(Genera.id, datasets))
    Maaslin2.var <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                           dimnames = list(Genera.id, datasets))

    ## LM-CLR
    LMCLR.est <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                        dimnames = list(Genera.id, datasets))
    LMCLR.var <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                        dimnames = list(Genera.id, datasets))

    ## LinDA
    LinDA.est <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                        dimnames = list(Genera.id, datasets))
    LinDA.var <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                        dimnames = list(Genera.id, datasets))

    ## Main analysis
    for(d in datasets){
      if(cmpds.id %in% colnames(data.for.lm[[d]]$rel.cmpd)){
        sample.kep <- intersect(rownames(data.for.lm[[d]]$rel.abd),
                                rownames(data.for.lm[[d]]$rel.cmpd)[!is.na(data.for.lm[[d]]$rel.cmpd[,cmpds.id])])

        tmp.otu <- matrix(0, nrow = length(sample.kep), ncol = length(Genera.id),
                          dimnames = list(sample.kep, Genera.id))
        inter.taxa <- colnames(null.obj[[d]]$Y_I)
        tmp.otu[sample.kep, inter.taxa] <- as.matrix(data.for.lm[[d]]$rel.abd[sample.kep,inter.taxa])
        tmp.cmpd <- matrix(0, nrow = length(sample.kep), ncol = 1, dimnames = list(sample.kep, cmpds.id))
        tmp.cmpd[sample.kep, cmpds.id] <- as.matrix(data.for.lm[[d]]$rel.cmpd[sample.kep,cmpds.id])

        ## remove zero observation taxa
        tmp.otu <- tmp.otu[,colSums(tmp.otu)!=0]

        ## Combined
        meta.data <- data.frame(Subject = (cluster[[d]])[sample.kep,"Subject"],
                                Study.Group = covariates.disease[[d]][match(sample.kep, covariates.disease[[d]]$Sample),] %>%
                                  data.frame(row.names = NULL) %>% tibble::column_to_rownames(var = "Sample"),
                                labels = tmp.cmpd[sample.kep, 1])


        taxa.id <- colnames(tmp.otu)
        tmp.otu.maaslin2 <- tmp.otu
        colnames(tmp.otu.maaslin2) <- paste0("Tax", 1:length(taxa.id))

        if(length(unique(meta.data$Subject)) == length(meta.data$Subject)){
          if(length(unique(meta.data$Study.Group)) == 1){
            ## ANCOM-BC2 model
            ANCOMBC2.model <- ancombc.fun(feature.table = t(round(tmp.otu)),
                                          meta = meta.data,
                                          formula = "labels",
                                          adjust.method = "fdr",
                                          group = NULL,
                                          subject = NULL,
                                          method = "ancombc2")

            ## Maaslin2 model
            Maaslin2.model <- Maaslin2::Maaslin2(input_data = tmp.otu.maaslin2,
                                                 input_metadata = meta.data,
                                                 output = "./",
                                                 min_abundance = 0,
                                                 min_prevalence = 0,
                                                 normalization = "TSS",
                                                 transform = "AST",
                                                 analysis_method = "LM",
                                                 max_significance = 1,
                                                 random_effects = NULL,
                                                 fixed_effects = "labels",
                                                 standardize = FALSE,
                                                 plot_heatmap = FALSE,
                                                 plot_scatter = FALSE)

            ## LM-CLR
            Z.pool <- log(tmp.otu + 0.5)
            clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

            inormal <- apply(clrx_Z, 2, function(d){
              data.d <- data.frame(meta.data, d = d)
              lm.fit <- lm(d ~ labels, data = data.d)
              sum.fit <- summary(lm.fit)
              if("labels" %in% rownames(sum.fit$coefficients)){
                return(sum.fit$coefficients["labels", c("Estimate", "Std. Error")])
              }else{
                return(c(NA, NA))
              }
            })

            ## LinDA
            Linda.model <- MicrobiomeStat::linda(feature.dat = t(tmp.otu),
                                                 meta.dat = meta.data,
                                                 formula = '~labels',
                                                 feature.dat.type = "count",
                                                 prev.filter = 0,
                                                 adaptive = TRUE,
                                                 alpha = 0.05)
          }else{
            ## ANCOM-BC2
            ANCOMBC2.model <- ancombc.fun(feature.table = t(round(tmp.otu)),
                                          meta = meta.data,
                                          formula = "Study.Group + labels",
                                          adjust.method = "fdr",
                                          group = NULL,
                                          subject = NULL,
                                          method = "ancombc2")

            ## Maaslin2
            Maaslin2.model <- Maaslin2::Maaslin2(input_data = tmp.otu.maaslin2,
                                                 input_metadata = meta.data,
                                                 output = "./",
                                                 min_abundance = 0,
                                                 min_prevalence = 0,
                                                 normalization = "TSS",
                                                 transform = "AST",
                                                 analysis_method = "LM",
                                                 max_significance = 1,
                                                 random_effects = NULL,
                                                 fixed_effects = c("Study.Group", "labels"),
                                                 standardize = FALSE,
                                                 plot_heatmap = FALSE,
                                                 plot_scatter = FALSE)

            ## LM-CLR
            Z.pool <- log(tmp.otu + 0.5)
            clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

            inormal <- apply(clrx_Z, 2, function(d){
              data.d <- data.frame(meta.data, d = d)
              lm.fit <- lm(d ~ labels + Study.Group, data = data.d)
              sum.fit <- summary(lm.fit)
              if("labels" %in% rownames(sum.fit$coefficients)){
                return(sum.fit$coefficients["labels",c("Estimate", "Std. Error")])
              }else{
                return(c(NA, NA))
              }
            })

            ## LinDA
            Linda.model <- MicrobiomeStat::linda(feature.dat = t(tmp.otu),
                                                 meta.dat = meta.data,
                                                 formula = '~labels+Study.Group',
                                                 feature.dat.type = "count",
                                                 prev.filter = 0,
                                                 adaptive = TRUE,
                                                 alpha = 0.05)
          }
        }else{
          if(length(unique(meta.data$Study.Group)) == 1){
            ## ANCOM-BC2
            ANCOMBC2.model <- ancombc.fun(feature.table = t(round(tmp.otu)),
                                          meta = meta.data,
                                          formula = "labels",
                                          adjust.method = "fdr",
                                          group = NULL,
                                          subject = "Subject",
                                          method = "ancombc2")

            ## Maaslin2
            Maaslin2.model <- Maaslin2::Maaslin2(input_data = tmp.otu.maaslin2,
                                                 input_metadata = meta.data,
                                                 output = "./",
                                                 min_abundance = 0,
                                                 min_prevalence = 0,
                                                 normalization = "TSS",
                                                 transform = "AST",
                                                 analysis_method = "LM",
                                                 max_significance = 1,
                                                 random_effects = "Subject",
                                                 fixed_effects = "labels",
                                                 standardize = FALSE,
                                                 plot_heatmap = FALSE,
                                                 plot_scatter = FALSE)

            ## LM-CLR
            Z.pool <- log(tmp.otu + 0.5)
            clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

            inormal <- apply(clrx_Z, 2, function(d){
              data.d <- data.frame(meta.data, d = d)
              lm.fit <- lme4::lmer(d ~ labels + (1 | Subject), data = data.d)
              sum.fit <- summary(lm.fit)
              if("labels" %in% rownames(sum.fit$coefficients)){
                return(sum.fit$coefficients["labels",c("Estimate", "Std. Error")])
              }else{
                return(c(NA, NA))
              }
            })

            ## LinDA
            Linda.model <- MicrobiomeStat::linda(feature.dat = t(tmp.otu),
                                                 meta.dat = meta.data,
                                                 formula = '~labels + (1|Subject)',
                                                 feature.dat.type = "count",
                                                 prev.filter = 0,
                                                 adaptive = TRUE,
                                                 alpha = 0.05)
          }else{
            ## ANCOM-BC2
            ANCOMBC2.model <- ancombc.fun(feature.table = t(round(tmp.otu)),
                                          meta = meta.data,
                                          formula = "Study.Group + labels",
                                          adjust.method = "fdr",
                                          group = NULL,
                                          subject = "Subject",
                                          method = "ancombc2")

            ## Maaslin2
            Maaslin2.model <- Maaslin2::Maaslin2(input_data = tmp.otu.maaslin2,
                                                 input_metadata = meta.data,
                                                 output = "./",
                                                 min_abundance = 0,
                                                 min_prevalence = 0,
                                                 normalization = "TSS",
                                                 transform = "AST",
                                                 analysis_method = "LM",
                                                 max_significance = 1,
                                                 random_effects = "Subject",
                                                 fixed_effects = c("Study.Group", "labels"),
                                                 standardize = FALSE,
                                                 plot_heatmap = FALSE,
                                                 plot_scatter = FALSE)

            ## LM-CLR
            Z.pool <- log(tmp.otu + 0.5)
            clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))

            inormal <- apply(clrx_Z, 2, function(d){
              data.d <- data.frame(meta.data, d = d)
              lm.fit <- lme4::lmer(d ~ labels + (1 | Subject) + Study.Group, data = data.d)
              sum.fit <- summary(lm.fit)
              if("labels" %in% rownames(sum.fit$coefficients)){
                return(sum.fit$coefficients["labels",c("Estimate", "Std. Error")])
              }else{
                return(c(NA, NA))
              }
            })

            ## LinDA
            Linda.model <- MicrobiomeStat::linda(feature.dat = t(tmp.otu),
                                                 meta.dat = meta.data,
                                                 formula = '~labels+Study.Group+(1|Subject)',
                                                 feature.dat.type = "count",
                                                 prev.filter = 0,
                                                 adaptive = TRUE,
                                                 alpha = 0.05)
          }
        }
        res.tmp <- Maaslin2.model$result[Maaslin2.model$results$metadata == "labels",]
        res.tmp <- res.tmp[order(as.numeric(sapply(strsplit(res.tmp$feature, split = "Tax"), "[[", 2, simplify = TRUE))),]
        res.tmp$feature <- taxa.id
        Maaslin2.single <- res.tmp

        ## ANCOM-BC2
        ANCOMBC2.est[ANCOMBC2.model$res$taxon[ANCOMBC2.model$res$passed_ss_labels],d] <-
          ANCOMBC2.model$res$lfc_labels[ANCOMBC2.model$res$passed_ss_labels]
        ANCOMBC2.var[ANCOMBC2.model$res$taxon[ANCOMBC2.model$res$passed_ss_labels],d] <-
          (ANCOMBC2.model$res$se_labels[ANCOMBC2.model$res$passed_ss_labels])^2

        ## Maaslin2
        Maaslin2.est[Maaslin2.single$feature,d] <-  Maaslin2.single$coef
        Maaslin2.var[Maaslin2.single$feature,d] <- (Maaslin2.single$stderr)^2

        ## LM-CLR
        LMCLR.est[colnames(inormal),d] <- inormal[1,]
        LMCLR.var[colnames(inormal),d] <- (inormal[2,])^2

        ## LinDA
        LinDA.est[rownames(Linda.model$output$labels),d] <- Linda.model$output$labels$log2FoldChange
        LinDA.var[rownames(Linda.model$output$labels),d] <- (Linda.model$output$labels$lfcSE)^2
      }
    }

    ANCOMBC2.model <- meta_analysis(AA.est = ANCOMBC2.est, AA.var = ANCOMBC2.var, rm.method = "EE")

    Maaslin2.model <- meta_analysis(AA.est = Maaslin2.est, AA.var = Maaslin2.var, rm.method = "EE")

    LMCLR.model <- meta_analysis(AA.est = LMCLR.est, AA.var = LMCLR.var, rm.method = "EE")

    Linda.model <- meta_analysis(AA.est = LinDA.est, AA.var = LinDA.var, rm.method = "EE")

    ## ANCOM-BC2
    ANCOMBC2.meta.pval <- ANCOMBC2.meta.pval %>%
      dplyr::left_join(ANCOMBC2.model %>%
                         dplyr::transmute(feature, pval) %>%
                         dplyr::rename_with(~cmpds.id, pval), by = "feature")

    ANCOMBC2.meta.pval.het <- ANCOMBC2.meta.pval.het %>%
      dplyr::left_join(ANCOMBC2.model %>%
                         dplyr::transmute(feature, pval.het) %>%
                         dplyr::rename_with(~cmpds.id, pval.het), by = "feature")

    ANCOMBC2.meta.coef <- ANCOMBC2.meta.coef %>%
      dplyr::left_join(ANCOMBC2.model %>%
                         dplyr::transmute(feature, coef) %>%
                         dplyr::rename_with(~cmpds.id, coef), by = "feature")

    ANCOMBC2.meta.sd <- ANCOMBC2.meta.sd %>%
      dplyr::left_join(ANCOMBC2.model %>%
                         dplyr::transmute(feature, sd) %>%
                         dplyr::rename_with(~cmpds.id, sd), by = "feature")

    ## Maaslin2
    Maaslin2.meta.pval <- Maaslin2.meta.pval %>%
      dplyr::left_join(Maaslin2.model %>%
                         dplyr::transmute(feature, pval) %>%
                         dplyr::rename_with(~cmpds.id, pval), by = "feature")

    Maaslin2.meta.pval.het <- Maaslin2.meta.pval.het %>%
      dplyr::left_join(Maaslin2.model %>%
                         dplyr::transmute(feature, pval.het) %>%
                         dplyr::rename_with(~cmpds.id, pval.het), by = "feature")

    Maaslin2.meta.coef <- Maaslin2.meta.coef %>%
      dplyr::left_join(Maaslin2.model %>%
                         dplyr::transmute(feature, coef) %>%
                         dplyr::rename_with(~cmpds.id, coef), by = "feature")

    Maaslin2.meta.sd <- Maaslin2.meta.sd %>%
      dplyr::left_join(Maaslin2.model %>%
                         dplyr::transmute(feature, sd) %>%
                         dplyr::rename_with(~cmpds.id, sd), by = "feature")

    ## LM-CLR
    LMCLR.meta.pval <- LMCLR.meta.pval %>%
      dplyr::left_join(LMCLR.model %>%
                         dplyr::transmute(feature, pval) %>%
                         dplyr::rename_with(~cmpds.id, pval), by = "feature")

    LMCLR.meta.pval.het <- LMCLR.meta.pval.het %>%
      dplyr::left_join(LMCLR.model %>%
                         dplyr::transmute(feature, pval.het) %>%
                         dplyr::rename_with(~cmpds.id, pval.het), by = "feature")

    LMCLR.meta.coef <- LMCLR.meta.coef %>%
      dplyr::left_join(LMCLR.model %>%
                         dplyr::transmute(feature, coef) %>%
                         dplyr::rename_with(~cmpds.id, coef), by = "feature")

    LMCLR.meta.sd <- LMCLR.meta.sd %>%
      dplyr::left_join(LMCLR.model %>%
                         dplyr::transmute(feature, sd) %>%
                         dplyr::rename_with(~cmpds.id, sd), by = "feature")

    ## LinDA
    Linda.meta.pval <- Linda.meta.pval %>%
      dplyr::left_join(Linda.model %>%
                         dplyr::transmute(feature, pval) %>%
                         dplyr::rename_with(~cmpds.id, pval), by = "feature")

    Linda.meta.pval.het <- Linda.meta.pval.het %>%
      dplyr::left_join(Linda.model %>%
                         dplyr::transmute(feature, pval.het) %>%
                         dplyr::rename_with(~cmpds.id, pval.het), by = "feature")

    Linda.meta.coef <- Linda.meta.coef %>%
      dplyr::left_join(Linda.model %>%
                         dplyr::transmute(feature, coef) %>%
                         dplyr::rename_with(~cmpds.id, coef), by = "feature")

    Linda.meta.sd <- Linda.meta.sd %>%
      dplyr::left_join(Linda.model %>%
                         dplyr::transmute(feature, sd) %>%
                         dplyr::rename_with(~cmpds.id, sd), by = "feature")

    ## Save data
    save(ANCOMBC2.meta.pval.het, ANCOMBC2.meta.pval, ANCOMBC2.meta.coef, ANCOMBC2.meta.sd,
         Maaslin2.meta.pval.het, Maaslin2.meta.pval, Maaslin2.meta.coef, Maaslin2.meta.sd,
         LMCLR.meta.pval.het, LMCLR.meta.pval, LMCLR.meta.coef, LMCLR.meta.sd,
         Linda.meta.pval.het, Linda.meta.pval, Linda.meta.coef, Linda.meta.sd,

         ANCOMBC2.est, ANCOMBC2.var,
         Maaslin2.est, Maaslin2.var,
         LMCLR.est, LMCLR.var,
         LinDA.est, LinDA.var,
         file = paste0("./Metabolites/Output/Metabolites/Compare_method_", cmpds.id, ".Rdata"))
  }



