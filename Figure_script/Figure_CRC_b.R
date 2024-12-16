# =============================================== #
#       CRC Data: Figure b, scatter plots         #
# =============================================== #

  # Packages
  library("dplyr")

  # General ----
  rm(list = ls())

  # Load models and data
  load("./CRC/Output/CRC_output.Rdata")
  load("./CRC/Processed_data/data.org.K401.Rdata")
  L <- length(data.rel)

  ## Align figures
  pdf("./Figure/Figure_CRC_b.pdf", width = 12.45, height = 3.13, bg = "white")

  par(mfrow = c(1, 4), mar = c(2, 2, 2, 0.2), oma = c(0, 0, 0, 0))

  ## PALM v.s. ANCOM-BC2
  data.plot <- PALM.res %>% dplyr::transmute(feature = feature, pval.y = -log10(het.pval)) %>%
    dplyr::left_join(ANCOMBC2.res %>% dplyr::transmute(feature = features, pval.x = -log10(het.pval)), by = "feature")

  plot(data.plot$pval.x, data.plot$pval.y,
       ylab = "",#expression(atop("Het. test " ~ -log10[10](p-val) ,"in PALM")),
       xlab = "",
       xlim = c(0, 5), ylim = c(0, 5), col = "blue", pch = 18, mgp=c(2,0.5,0),
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0, 1, col = "brown", lty = 2)
  #mtext(expression(atop("Het. test " ~ -log10[10](p-val) ,"in ANCOM-BC2")), side = 1, line = 4)

  ## PALM v.s. LinDA
  data.plot <- PALM.res %>% dplyr::transmute(feature = feature, pval.y = -log10(het.pval)) %>%
    dplyr::left_join(Linda.res %>% dplyr::transmute(feature = features, pval.x = -log10(het.pval)), by = "feature") %>%
    dplyr::filter(!is.na(pval.x))

  plot(data.plot$pval.x, data.plot$pval.y,
       ylab = "",
       xlab = "",
       xlim = c(0, 5), ylim = c(0, 5), col = "skyblue", pch = 18, mgp=c(2,0.5,0), yaxt = "n",
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0, 1, col = "brown", lty = 2)
  # mtext(expression(atop("Het. test " ~ -log10[10](p-val) ,"in LinDA")), side = 1, line = 4)

  ## PALM v.s. LM-CLR
  data.plot <- PALM.res %>% dplyr::transmute(feature = feature, pval.y = -log10(het.pval)) %>%
    dplyr::left_join(lmclr.res %>% dplyr::transmute(feature = features, pval.x = -log10(het.pval)), by = "feature") %>%
    dplyr::filter(!is.na(pval.x))

  plot(data.plot$pval.x, data.plot$pval.y, ylab = " ", xlab = "",
       xlim = c(0, 5), ylim = c(0, 5), col = "orange", pch = 18, mgp=c(2,0.5,0),  yaxt = "n",
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0, 1, col = "brown", lty = 2)
  #mtext(expression(atop("Het. test " ~ -log10[10](p-val) ,"in LM-CLR")), side = 1, line = 4)

  ## PALM v.s. MaAsLin2
  data.plot <- PALM.res %>% dplyr::transmute(feature = feature, pval.y = -log10(het.pval)) %>%
    dplyr::left_join(Maaslin2.res %>% dplyr::transmute(feature = features, pval.x = -log10(het.pval)), by = "feature") %>%
    dplyr::filter(!is.na(pval.x))

  plot(data.plot$pval.x, data.plot$pval.y, ylab = " ", xlab = "",
       xlim = c(0, 5), ylim = c(0, 5), col = "#4dac26", pch = 18, mgp=c(2,0.5,0),  yaxt = "n",
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0, 1, col = "brown", lty = 2)
  #mtext(expression(atop("Het. test " ~ -log10[10](p-val) ,"in MaAsLin2")), side = 1, line = 4)

  dev.off()


