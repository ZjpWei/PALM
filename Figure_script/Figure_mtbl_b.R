# =============================================== #
#              MTBL Data: Figure a                #
# =============================================== #

  # Package
  library(tidyverse)

  # General ----
  rm(list = ls())
  load("./Metabolites/Output/MTBL_PALM.Rdata")
  load("./Metabolites/Output/MTBL_others.Rdata")
  target.fdr <- 0.05

  ## PALM
  PALM_het_mat <- PALM.meta.pval.het %>% tibble::column_to_rownames(var = "feature")
  PALM.count <- colMeans(apply(PALM_het_mat, 2, p.adjust, method = "fdr") <= 0.1, na.rm = TRUE)

  ## ANCOM-BC2
  ANCOMBC2_het_mat <- ANCOMBC2.meta.pval.het %>% tibble::column_to_rownames(var = "feature")
  ANCOMBC2.count <- colMeans(apply(ANCOMBC2_het_mat, 2, p.adjust, method = "fdr") <= 0.1, na.rm = TRUE)

  ## LinDA
  Linda_het_mat <- Linda.meta.pval.het %>% tibble::column_to_rownames(var = "feature")
  Linda.count <- colMeans(apply(Linda_het_mat, 2, p.adjust, method = "fdr") <= 0.1, na.rm = TRUE)

  ## LM-CLR
  LMCLR_het_mat <- LMCLR.meta.pval.het %>% tibble::column_to_rownames(var = "feature")
  LMCLR.count <- colMeans(apply(LMCLR_het_mat, 2, p.adjust, method = "fdr") <= 0.1, na.rm = TRUE)

  ## Maaslin2
  Maaslin2_het_mat <- Maaslin2.meta.pval.het %>% tibble::column_to_rownames(var = "feature")
  Maaslin2.count <- colMeans(apply(Maaslin2_het_mat, 2, p.adjust, method = "fdr") <= 0.1, na.rm = TRUE)

  ## plot
  pdf("./Figure/Figure_mtbl_b.pdf", width = 12.45, height = 3.13, bg = "white")

  par(mfrow = c(1, 4), mar = c(2, 2.5, 2, 0.2), oma = c(0, 0, 0, 0))

  plot(ANCOMBC2.count, PALM.count, col = "blue", pch = 18, ylab = "", xlab = "",
       xlim = c(0, 1), ylim = c(0, 1), cex.axis = 1.6,  cex.lab = 1.6)

  abline(0,1,col="brown", lty = 2)

  #mtext("Prop. Het. taxa \nin ANCOM-BC2", side = 1, line = 4)

  plot(Linda.count, PALM.count, col = "skyblue", pch = 18, ylab = "", xlab = "",
       xlim = c(0, 1), ylim = c(0, 1), cex.axis = 1.6,  cex.lab = 1.6, yaxt = "n")

  abline(0,1,col="brown", lty = 2)

  #mtext("Prop. Het. taxa \n in LinDA", side = 1, line = 4)

  plot(LMCLR.count, PALM.count, col = "orange", pch = 18, ylab = "", xlab = "",
       xlim = c(0, 1), ylim = c(0, 1), cex.axis = 1.6,  cex.lab = 1.6, yaxt = "n")

  abline(0,1,col="brown", lty = 2)

  #mtext("Prop. Het. taxa \n in LM-CLR", side = 1, line = 4)

  plot(Maaslin2.count, PALM.count,  col = "#4dac26", pch = 18, ylab = "", xlab = "",
       xlim = c(0, 1), ylim = c(0, 1), cex.axis = 1.6,  cex.lab = 1.6, yaxt = "n")

  abline(0,1,col="brown", lty = 2)

  #mtext("Prop. Het. taxa \n MaAsLin2", side = 1, line = 4)

  dev.off()
