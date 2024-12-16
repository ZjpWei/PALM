# =============================================== #
#       CRC Data: Figure a, Venn diagram          #
# =============================================== #

  # Package
  library("VennDiagram")

  # General ----
  rm(list = ls())

  # Load models and data
  load("./CRC/Processed_data/data.org.K401.Rdata")
  load("./CRC/Output/CRC_output.Rdata")
  L <- length(data.rel)

  ## Match FDR 0.05
  fdr.cut <- 0.05

  # Get selected taxa for each method (other methods match the number selected by Melody)
  index <- list()
  taxa.id <- PALM.res$feature[PALM.res$qval <= fdr.cut]
  len.tax <- length(taxa.id)

  ## Selected taxa
  Venn.set <- list("PALM" = taxa.id,
                   "ANCOM-BC2" = ANCOMBC2.res$features[ANCOMBC2.res$qval <= fdr.cut & !is.na(ANCOMBC2.res$qval)],
                   "MaAsLin2" = Maaslin2.res$features[Maaslin2.res$qval <= fdr.cut & !is.na(Maaslin2.res$qval)],
                   "LM-CLR" = lmclr.res$features[lmclr.res$qval <= fdr.cut & !is.na(lmclr.res$qval)],
                   "LinDA" = Linda.res$features[Linda.res$qval <= fdr.cut & !is.na(Linda.res$qval)])

  venn.diagram(Venn.set, filename = "./Figure/Figure_CRC_a.png",
               category.names = c("" , "" , "", "", ""),
               cex = 1.7,
               width = 3000, height = 3000,
               fill = c("red", "blue", "#4dac26", "orange", "skyblue"),
               cat.dist = c(0.1, 0.06, 0.09, 0.1, 0.1))
