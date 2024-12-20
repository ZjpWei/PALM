# =============================================== #
#   (2) Analyze microbiome-metabolome data (PALM) #
# =============================================== #

  # Packages ----
  library('tidyverse')
  library('PALM')

  rm(list = ls())

  ## Load data
  load("./Metabolites/Processed_data/processed_data.Rdata")
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

  # PALM analysis
  otu_data_lst <- list()
  cmpd_data_lst <- list()
  cluster_data_lst <- list()
  covariates_adjust_lst <- list()
  for(d in datasets){
    sample.kep <- intersect(rownames(data.for.lm[[d]]$rel.abd), rownames(data.for.lm[[d]]$rel.cmpd))
    covariates_adjust_lst[[d]] <- covariates.disease[[d]][match(sample.kep, covariates.disease[[d]]$Sample),] %>%
      data.frame(row.names = NULL) %>%
      tibble::column_to_rownames(var = "Sample")
    otu_data_lst[[d]] <- as.matrix(data.for.lm[[d]]$rel.abd %>% dplyr::slice(match(sample.kep, rownames(data.for.lm[[d]]$rel.abd))))
    cmpd_data_lst[[d]] <- as.matrix(data.for.lm[[d]]$rel.cmpd %>% dplyr::slice(match(sample.kep, rownames(data.for.lm[[d]]$rel.cmpd))))
    cluster_data_lst[[d]] <- (cluster[[d]])[sample.kep,]
    names(cluster_data_lst[[d]]) <- sample.kep
  }
  covariates_adjust_lst$POYET_BIO_ML_2019 <- NULL

  ## Save preprocessing
  save.image("./Metabolites/Processed_data/MTBL.RData")

  # PALM meta-analysis
  null.obj <- palm.null.model(rel.abd = otu_data_lst, covariate.adjust = covariates_adjust_lst)

  save(null.obj, file = "./Metabolites/Output/PALM_null.Rdata")

  summary.stat.study.all <- palm.get.summary(null.obj = null.obj, covariate.interest = cmpd_data_lst,
                                             cluster = cluster_data_lst)

  PALM.test <- palm.test(summary.stats = summary.stat.study.all, p.adjust.method = "fdr")

  # PALM model ----
  PALM.meta.pval <- data.frame(feature = Genera.id)
  PALM.meta.pval.het <- data.frame(feature = Genera.id)
  PALM.meta.coef <- data.frame(feature = Genera.id)
  PALM.meta.sd <- data.frame(feature = Genera.id)

  for(cmpds.id in Metabolite.ids){
    PALM.est <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                     dimnames = list(Genera.id, datasets))
    PALM.var <- matrix(NA, nrow = length(Genera.id), ncol = length(datasets),
                     dimnames = list(Genera.id, datasets))

    for(d in datasets){
      if(cmpds.id %in% colnames(data.for.lm[[d]]$rel.cmpd)){
        PALM.est[PALM.test[[cmpds.id]]$palm_fits[[d]]$feature,d] <- PALM.test[[cmpds.id]]$palm_fits[[d]]$coef
        PALM.var[PALM.test[[cmpds.id]]$palm_fits[[d]]$feature,d] <- (PALM.test[[cmpds.id]]$palm_fits[[d]]$stderr)^2
      }
    }

    ## Calculate FDR
    meta.coef <- PALM.test[[cmpds.id]]$meta_fits$coef
    meta.sd <- PALM.test[[cmpds.id]]$meta_fits$stderr
    pval.sin <- PALM.test[[cmpds.id]]$meta_fits$pval
    tmp.pval <- NULL
    for(k in 1:nrow(PALM.est)){
      nonna.id <- !is.na(PALM.est[k,])
      if(sum(nonna.id) > 1){
        m <- metafor::rma(yi = PALM.est[k,nonna.id], vi = PALM.var[k,nonna.id],
                          control=list(maxiter=2000), method = "EE")
        tmp.pval <- c(tmp.pval, m$QEp)
      }else{
        tmp.pval <- c(tmp.pval, NA)
      }
    }

    PALM.model <- data.frame(feature = Genera.id, pval = pval.sin, pval.het = tmp.pval,
                             coef = meta.coef, sd = meta.sd)

    PALM.meta.pval <- PALM.meta.pval %>%
      dplyr::left_join(PALM.model %>%
                         dplyr::transmute(feature, pval) %>%
                         dplyr::rename_with(~cmpds.id, pval), by = "feature")

    PALM.meta.pval.het <- PALM.meta.pval.het %>%
      dplyr::left_join(PALM.model %>%
                         dplyr::transmute(feature, pval.het) %>%
                         dplyr::rename_with(~cmpds.id, pval.het), by = "feature")

    PALM.meta.coef <- PALM.meta.coef %>%
      dplyr::left_join(PALM.model %>%
                         dplyr::transmute(feature, coef) %>%
                         dplyr::rename_with(~cmpds.id, coef), by = "feature")

    PALM.meta.sd <- PALM.meta.sd %>%
      dplyr::left_join(PALM.model %>%
                         dplyr::transmute(feature, sd) %>%
                         dplyr::rename_with(~cmpds.id, sd), by = "feature")

  }

  save(PALM.meta.pval, PALM.meta.pval.het, PALM.meta.coef, PALM.meta.sd,
       PALM.est, PALM.var, PALM.test,
       file = "./Metabolites/Output/MTBL_PALM.Rdata")

