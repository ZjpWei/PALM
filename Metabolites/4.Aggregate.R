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
      ANCOMBC2.meta.sd.all <- ANCOMBC2.meta.sd

      ## LM-INT
      LMCLR.meta.coef.all <- LMCLR.meta.coef
      LMCLR.meta.pval.all <- LMCLR.meta.pval
      LMCLR.meta.pval.het.all <- LMCLR.meta.pval.het
      LMCLR.meta.sd.all <- LMCLR.meta.sd

      ## Linda
      Linda.meta.coef.all <- Linda.meta.coef
      Linda.meta.pval.all <- Linda.meta.pval
      Linda.meta.pval.het.all <- Linda.meta.pval.het
      Linda.meta.sd.all <- Linda.meta.sd

      ## Maaslin2
      Maaslin2.meta.coef.all <- Maaslin2.meta.coef
      Maaslin2.meta.pval.all <- Maaslin2.meta.pval
      Maaslin2.meta.pval.het.all <- Maaslin2.meta.pval.het
      Maaslin2.meta.sd.all <- Maaslin2.meta.sd

    }else{
      ## ANCOM-BC2
      ANCOMBC2.meta.coef.all <- ANCOMBC2.meta.coef.all %>%
        dplyr::left_join(ANCOMBC2.meta.coef, by = "feature")
      ANCOMBC2.meta.pval.all <- ANCOMBC2.meta.pval.all %>%
        dplyr::left_join(ANCOMBC2.meta.pval, by = "feature")
      ANCOMBC2.meta.pval.het.all <- ANCOMBC2.meta.pval.het.all %>%
        dplyr::left_join(ANCOMBC2.meta.pval.het, by = "feature")
      ANCOMBC2.meta.sd.all <- ANCOMBC2.meta.sd.all %>%
        dplyr::left_join(ANCOMBC2.meta.sd, by = "feature")
      ## LM-INT
      LMCLR.meta.coef.all <- LMCLR.meta.coef.all %>%
        dplyr::left_join(LMCLR.meta.coef, by = "feature")
      LMCLR.meta.pval.all <- LMCLR.meta.pval.all %>%
        dplyr::left_join(LMCLR.meta.pval, by = "feature")
      LMCLR.meta.pval.het.all <- LMCLR.meta.pval.het.all %>%
        dplyr::left_join(LMCLR.meta.pval.het, by = "feature")
      LMCLR.meta.sd.all <- LMCLR.meta.sd.all %>%
        dplyr::left_join(LMCLR.meta.sd, by = "feature")

      ## Linda
      Linda.meta.coef.all <- Linda.meta.coef.all %>%
        dplyr::left_join(Linda.meta.coef, by = "feature")
      Linda.meta.pval.all <- Linda.meta.pval.all %>%
        dplyr::left_join(Linda.meta.pval, by = "feature")
      Linda.meta.pval.het.all <- Linda.meta.pval.het.all %>%
        dplyr::left_join(Linda.meta.pval.het, by = "feature")
      Linda.meta.sd.all <- Linda.meta.sd.all %>%
        dplyr::left_join(Linda.meta.sd, by = "feature")

      ## Maaslin2
      Maaslin2.meta.coef.all <- Maaslin2.meta.coef.all %>%
        dplyr::left_join(Maaslin2.meta.coef, by = "feature")
      Maaslin2.meta.pval.all <- Maaslin2.meta.pval.all %>%
        dplyr::left_join(Maaslin2.meta.pval, by = "feature")
      Maaslin2.meta.pval.het.all <- Maaslin2.meta.pval.het.all %>%
        dplyr::left_join(Maaslin2.meta.pval.het, by = "feature")
      Maaslin2.meta.sd.all <- Maaslin2.meta.sd.all %>%
        dplyr::left_join(Maaslin2.meta.sd, by = "feature")
    }
    rm(list = c("ANCOMBC2.meta.coef", "ANCOMBC2.meta.pval", "ANCOMBC2.meta.pval.het",
                "Linda.meta.coef", "Linda.meta.pval", "Linda.meta.pval.het",
                "LMCLR.meta.coef", "LMCLR.meta.pval", "LMCLR.meta.pval.het",
                "Maaslin2.meta.coef", "Maaslin2.meta.pval", "Maaslin2.meta.pval.het"))
  }

  ## ANCOM-BC2
  ANCOMBC2.meta.coef <- ANCOMBC2.meta.coef.all
  ANCOMBC2.meta.pval <- ANCOMBC2.meta.pval.all
  ANCOMBC2.meta.pval.het <- ANCOMBC2.meta.pval.het.all
  ANCOMBC2.meta.sd <- ANCOMBC2.meta.sd.all

  ## Linda
  Linda.meta.coef <- Linda.meta.coef.all
  Linda.meta.pval <- Linda.meta.pval.all
  Linda.meta.pval.het <- Linda.meta.pval.het.all
  Linda.meta.sd <- Linda.meta.sd.all

  ## LM-INT
  LMCLR.meta.coef <- LMCLR.meta.coef.all
  LMCLR.meta.pval <- LMCLR.meta.pval.all
  LMCLR.meta.pval.het <- LMCLR.meta.pval.het.all
  LMCLR.meta.sd <- LMCLR.meta.sd.all

  ## MaAsLin2
  Maaslin2.meta.coef <- Maaslin2.meta.coef.all
  Maaslin2.meta.pval <- Maaslin2.meta.pval.all
  Maaslin2.meta.pval.het <- Maaslin2.meta.pval.het.all
  Maaslin2.meta.sd <- Maaslin2.meta.sd.all

  save(ANCOMBC2.meta.pval.het, ANCOMBC2.meta.pval, ANCOMBC2.meta.coef, ANCOMBC2.meta.sd,
       Linda.meta.pval.het, Linda.meta.pval, Linda.meta.coef, Linda.meta.sd,
       LMCLR.meta.pval.het, LMCLR.meta.pval, LMCLR.meta.coef,LMCLR.meta.sd,
       Maaslin2.meta.pval.het, Maaslin2.meta.pval, Maaslin2.meta.coef, Maaslin2.meta.sd,
       file = "./Metabolites/Output/MTBL_others.Rdata")


