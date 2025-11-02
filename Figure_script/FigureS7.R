# =============================================== #
#                   FigureS7                      #
# =============================================== #

  # Packages ----
  library("tidyverse")
  library("ggplot2")

  # Main ----
  rm(list = ls())
  load("./Metabolites/Output/MTBL_PALM.Rdata")
  load("./Metabolites/Output/MTBL_others.Rdata")
  load("./Metabolites/Processed_data/MTBL.RData")

  feature.ids <- PALM.meta.coef$feature
  MTBL.ids <- colnames(PALM.meta.coef)[-1]
  PALM_coef_mat <- PALM.meta.coef %>% tibble::column_to_rownames("feature")
  PALM_pval_mat <- PALM.meta.pval %>% tibble::column_to_rownames(var = "feature")
  PALM_qval_mat <- apply(PALM_pval_mat, 2, p.adjust, method = "fdr")
  PALM_coef_mat[is.na(PALM_coef_mat)] <- 0
  PALM_coef_mat[PALM_qval_mat > 0.05] <- 0

  ANCOMBC2_coef_mat <- ANCOMBC2.meta.coef %>% tibble::column_to_rownames(var = "feature")
  ANCOMBC2_pval_mat <- ANCOMBC2.meta.pval %>% tibble::column_to_rownames(var = "feature")
  ANCOMBC2_qval_mat <- apply(ANCOMBC2_pval_mat, 2, p.adjust, method = "fdr")
  ANCOMBC2_coef_mat[is.na(ANCOMBC2_coef_mat)] <- 0
  ANCOMBC2_coef_mat[ANCOMBC2_qval_mat > 0.05] <- 0

  Linda_coef_mat <- Linda.meta.coef %>% tibble::column_to_rownames(var = "feature")
  Linda_pval_mat <- Linda.meta.pval %>% tibble::column_to_rownames(var = "feature")
  Linda_qval_mat <- apply(Linda_pval_mat, 2, p.adjust, method = "fdr")
  Linda_coef_mat[is.na(Linda_coef_mat)] <- 0
  Linda_coef_mat[Linda_qval_mat > 0.05] <- 0

  LMCLR_coef_mat <- LMCLR.meta.coef %>% tibble::column_to_rownames(var = "feature")
  LMCLR_pval_mat <- LMCLR.meta.pval %>% tibble::column_to_rownames(var = "feature")
  LMCLR_qval_mat <- apply(LMCLR_pval_mat, 2, p.adjust, method = "fdr")
  LMCLR_coef_mat[is.na(LMCLR_coef_mat)] <- 0
  LMCLR_coef_mat[LMCLR_qval_mat > 0.05] <- 0

  Maaslin2_coef_mat <- Maaslin2.meta.coef %>% tibble::column_to_rownames(var = "feature")
  Maaslin2_pval_mat <- Maaslin2.meta.pval %>% tibble::column_to_rownames(var = "feature")
  Maaslin2_qval_mat <- apply(Maaslin2_pval_mat, 2, p.adjust, method = "fdr")
  Maaslin2_coef_mat[is.na(Maaslin2_coef_mat)] <- 0
  Maaslin2_coef_mat[Maaslin2_qval_mat > 0.05] <- 0

  ## PALM v.s. ANCOMBC2
  prop <- NULL
  method <- NULL
  MTBL <- NULL
  for(i in MTBL.ids){
    ## ANCOMBC2
    prop <- c(prop, length(intersect(which(PALM_coef_mat[,i]!=0), which(ANCOMBC2_coef_mat[,i]!=0))) / length(which(PALM_coef_mat[,i]!=0)))
    method <- c(method, "ANCOM-BC2")
    MTBL <- c(MTBL, i)

    ## LinDA
    prop <- c(prop, length(intersect(which(PALM_coef_mat[,i]!=0), which(Linda_coef_mat[,i]!=0))) / length(which(PALM_coef_mat[,i]!=0)))
    method <- c(method, "LinDA")
    MTBL <- c(MTBL, i)

    ## LMCLR
    prop <- c(prop, length(intersect(which(PALM_coef_mat[,i]!=0), which(LMCLR_coef_mat[,i]!=0))) / length(which(PALM_coef_mat[,i]!=0)))
    method <- c(method, "LM-CLR")
    MTBL <- c(MTBL, i)

    ## Maaslin2
    prop <- c(prop, length(intersect(which(PALM_coef_mat[,i]!=0), which(Maaslin2_coef_mat[,i]!=0))) / length(which(PALM_coef_mat[,i]!=0)))
    method <- c(method, "MaAsLin2")
    MTBL <- c(MTBL, i)
  }
  df.boxplot <- data.frame(prop = prop, method = method, MTBL = MTBL)
  df.boxplot <- df.boxplot[!is.na(df.boxplot$prop),]

  ## Generate figure ----
  pdf("./Figure/FigureS7.pdf", width = 6.26, height = 4.76, bg = "white")

  ggplot(df.boxplot, aes(x = method, y = prop)) + geom_boxplot() +
    labs(x = "Method", y = "Proportion of PALM features \n overlapping with other methods") +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          title = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = 0.1, size = 18),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_markdown(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none")

  dev.off()
