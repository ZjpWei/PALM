# =============================================== #
#              MTBL Data: Figure a                #
# =============================================== #

  ## packages ----
  library("tidyverse")

  rm(list = ls())
  load("./Metabolites/Output/MTBL_PALM.Rdata")
  load("./Metabolites/Output/MTBL_others.Rdata")
  target.fdr <- 0.05

  ## PALM
  PALM_coef_mat <- PALM.meta.coef %>% tibble::column_to_rownames(var = "feature")
  PALM_pval_mat <- PALM.meta.pval %>% tibble::column_to_rownames(var = "feature")
  PALM_qval_mat <- apply(PALM_pval_mat, 2, p.adjust, method = "fdr")
  PALM.count <- colSums(PALM_qval_mat <= target.fdr, na.rm = TRUE)

  ## ANCOM-BC2
  ANCOMBC2_coef_mat <- ANCOMBC2.meta.coef %>% tibble::column_to_rownames(var = "feature")
  ANCOMBC2_pval_mat <- ANCOMBC2.meta.pval %>% tibble::column_to_rownames(var = "feature")
  ANCOMBC2_qval_mat <- apply(ANCOMBC2_pval_mat, 2, p.adjust, method = "fdr")
  ANCOMBC2.count <- colSums(ANCOMBC2_qval_mat <= target.fdr, na.rm = TRUE)

  ## LinDA
  Linda_coef_mat <- Linda.meta.coef %>% tibble::column_to_rownames(var = "feature")
  Linda_pval_mat <- Linda.meta.pval %>% tibble::column_to_rownames(var = "feature")
  Linda_qval_mat <- apply(Linda_pval_mat, 2, p.adjust, method = "fdr")
  Linda.count <- colSums(Linda_qval_mat <= target.fdr, na.rm = TRUE)

  ## LM-CLR
  LMCLR_coef_mat <- LMCLR.meta.coef %>% tibble::column_to_rownames(var = "feature")
  LMCLR_pval_mat <- LMCLR.meta.pval %>% tibble::column_to_rownames(var = "feature")
  LMCLR_qval_mat <- apply(LMCLR_pval_mat, 2, p.adjust, method = "fdr")
  LMCLR.count <- colSums(LMCLR_qval_mat <= target.fdr, na.rm = TRUE)

  ## Maaslin2
  Maaslin2_coef_mat <- Maaslin2.meta.coef %>% tibble::column_to_rownames(var = "feature")
  Maaslin2_pval_mat <- Maaslin2.meta.pval %>% tibble::column_to_rownames(var = "feature")
  Maaslin2_qval_mat <- apply(Maaslin2_pval_mat, 2, p.adjust, method = "fdr")
  Maaslin2.count <- colSums(Maaslin2_qval_mat <= target.fdr, na.rm = TRUE)

  ## plot
  pdf("./Figure/Figure_mtbl_a.pdf", width = 12.45, height = 3.13, bg = "white")

  par(mfrow = c(1, 4), mar = c(2, 2.5, 2, 0.2), oma = c(0, 0, 0, 0))

  plot(ANCOMBC2.count, PALM.count, col = "blue", pch = 18, ylab = "", xlab = "",
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0,1,col="brown", lty = 2)

  plot(Linda.count, PALM.count, col = "skyblue", pch = 18, ylab = " ", xlab = "",
       cex.axis = 1.6,  cex.lab = 1.6, yaxt = "n")
  abline(0,1,col="brown", lty = 2)

  plot(LMCLR.count, PALM.count, col = "orange", pch = 18, ylab = " ", xlab = "",
       cex.axis = 1.6,  cex.lab = 1.6, yaxt = "n")
  abline(0,1,col="brown", lty = 2)

  plot(Maaslin2.count, PALM.count,  col = "#4dac26", pch = 18, ylab = " ", xlab = "",
       cex.axis = 1.6,  cex.lab = 1.6, yaxt = "n")
  abline(0,1,col="brown", lty = 2)

  dev.off()
