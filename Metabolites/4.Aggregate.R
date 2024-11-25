# =============================================== #
#        (4) Aggregate MTBL analysis              #
# =============================================== #

  # Packages ----
  library('tidyverse')

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

  ## Aggregate results
  for(s in Metabolite.ids){
    load(paste0("./Metabolites/Output/Metabolites/Compare_method_", s, ".Rdata"))
    if(s == Metabolite.ids[1]){
      ## ANCOM-BC2
      ANCOMBC2.meta.coef.all <- ANCOMBC2.meta.coef
      ANCOMBC2.meta.pval.all <- ANCOMBC2.meta.pval
      ANCOMBC2.meta.pval.het.all <- ANCOMBC2.meta.pval.het

      ## LM-INT
      LMCLR.meta.coef.all <- LMINT.meta.coef
      LMCLR.meta.pval.all <- LMINT.meta.pval
      LMCLR.meta.pval.het.all <- LMINT.meta.pval.het

      ## Linda
      Linda.meta.coef.all <- Linda.meta.coef
      Linda.meta.pval.all <- Linda.meta.pval
      Linda.meta.pval.het.all <- Linda.meta.pval.het

      ## Maaslin2
      Maaslin2.meta.coef.all <- Maaslin2.meta.coef
      Maaslin2.meta.pval.all <- Maaslin2.meta.pval
      Maaslin2.meta.pval.het.all <- Maaslin2.meta.pval.het
    }else{
      ## ANCOM-BC2
      ANCOMBC2.meta.coef.all <- ANCOMBC2.meta.coef.all %>%
        dplyr::left_join(ANCOMBC2.meta.coef, by = "feature")
      ANCOMBC2.meta.pval.all <- ANCOMBC2.meta.pval.all %>%
        dplyr::left_join(ANCOMBC2.meta.pval, by = "feature")
      ANCOMBC2.meta.pval.het.all <- ANCOMBC2.meta.pval.het.all %>%
        dplyr::left_join(ANCOMBC2.meta.pval.het, by = "feature")

      ## LM-INT
      LMCLR.meta.coef.all <- LMCLR.meta.coef.all %>%
        dplyr::left_join(LMINT.meta.coef, by = "feature")
      LMCLR.meta.pval.all <- LMCLR.meta.pval.all %>%
        dplyr::left_join(LMINT.meta.pval, by = "feature")
      LMCLR.meta.pval.het.all <- LMCLR.meta.pval.het.all %>%
        dplyr::left_join(LMINT.meta.pval.het, by = "feature")

      ## Linda
      Linda.meta.coef.all <- Linda.meta.coef.all %>%
        dplyr::left_join(Linda.meta.coef, by = "feature")
      Linda.meta.pval.all <- Linda.meta.pval.all %>%
        dplyr::left_join(Linda.meta.pval, by = "feature")
      Linda.meta.pval.het.all <- Linda.meta.pval.het.all %>%
        dplyr::left_join(Linda.meta.pval.het, by = "feature")

      ## Maaslin2
      Maaslin2.meta.coef.all <- Maaslin2.meta.coef.all %>%
        dplyr::left_join(Maaslin2.meta.coef, by = "feature")
      Maaslin2.meta.pval.all <- Maaslin2.meta.pval.all %>%
        dplyr::left_join(Maaslin2.meta.pval, by = "feature")
      Maaslin2.meta.pval.het.all <- Maaslin2.meta.pval.het.all %>%
        dplyr::left_join(Maaslin2.meta.pval.het, by = "feature")
    }
    rm(list = c("ANCOMBC2.meta.coef", "ANCOMBC2.meta.pval", "ANCOMBC2.meta.pval.het",
                "Linda.meta.coef", "Linda.meta.pval", "Linda.meta.pval.het",
                "LMINT.meta.coef", "LMINT.meta.pval", "LMINT.meta.pval.het",
                "Maaslin2.meta.coef", "Maaslin2.meta.pval", "Maaslin2.meta.pval.het"))
  }

  ## ANCOM-BC2
  ANCOMBC2.meta.coef <- ANCOMBC2.meta.coef.all
  ANCOMBC2.meta.pval <- ANCOMBC2.meta.pval.all
  ANCOMBC2.meta.pval.het <- ANCOMBC2.meta.pval.het.all

  ## Linda
  Linda.meta.coef <- Linda.meta.coef.all
  Linda.meta.pval <- Linda.meta.pval.all
  Linda.meta.pval.het <- Linda.meta.pval.het.all

  ## LM-INT
  LMCLR.meta.coef <- LMCLR.meta.coef.all
  LMCLR.meta.pval <- LMCLR.meta.pval.all
  LMCLR.meta.pval.het <- LMCLR.meta.pval.het.all

  ## MaAsLin2
  Maaslin2.meta.coef <- Maaslin2.meta.coef.all
  Maaslin2.meta.pval <- Maaslin2.meta.pval.all
  Maaslin2.meta.pval.het <- Maaslin2.meta.pval.het.all

  save(ANCOMBC2.meta.pval.het, ANCOMBC2.meta.pval, ANCOMBC2.meta.coef,
       Linda.meta.pval.het, Linda.meta.pval, Linda.meta.coef,
       LMCLR.meta.pval.het, LMCLR.meta.pval, LMCLR.meta.coef,
       Maaslin2.meta.pval.het, Maaslin2.meta.pval, Maaslin2.meta.coef,
       file = "./Metabolites/Output/MTBL_others.Rdata")
