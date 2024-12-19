# =============================================== #
#          Figure 4: CRC data analysis            #
# =============================================== #

  # Package ----
  library("VennDiagram")
  library("plotly")
  library("tidyr")
  library("stringr")
  library("dplyr")
  library("ggplot2")
  library("grid")
  library("ggtext")
  library("cowplot")

  # Figure a ----
  rm(list = ls())

  ## Load models and data
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

  venn.diagram(Venn.set, filename = "./Figure/Figure4_a.png",
               category.names = c("" , "" , "", "", ""),
               cex = 1.7,
               width = 3000, height = 3000,
               fill = c("red", "blue", "#4dac26", "orange", "skyblue"),
               cat.dist = c(0.1, 0.06, 0.09, 0.1, 0.1))

  # Figure b ----
  rm(list = ls())

  # Load models and data
  load("./CRC/Output/CRC_output.Rdata")
  load("./CRC/Processed_data/data.org.K401.Rdata")
  L <- length(data.rel)

  ## Align figures
  pdf("./Figure/Figure4_b.pdf", width = 12.45, height = 3.13, bg = "white")

  par(mfrow = c(1, 4), mar = c(2, 2, 2, 0.2), oma = c(0, 0, 0, 0))

  ## PALM v.s. ANCOM-BC2
  data.plot <- PALM.res %>% dplyr::transmute(feature = feature, pval.y = -log10(het.qval)) %>%
    dplyr::left_join(ANCOMBC2.res %>% dplyr::transmute(feature = features, pval.x = -log10(het.qval)), by = "feature")

  plot(data.plot$pval.x, data.plot$pval.y,
       ylab = "",#expression(atop("Het. test " ~ -log10[10](p-val) ,"in PALM")),
       xlab = "",
       col = "blue", pch = 18, mgp=c(2,0.5,0),
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0, 1, col = "brown", lty = 2)
  #mtext(expression(atop("Het. test " ~ -log10[10](p-val) ,"in ANCOM-BC2")), side = 1, line = 4)

  ## PALM v.s. LinDA
  data.plot <- PALM.res %>% dplyr::transmute(feature = feature, pval.y = -log10(het.qval)) %>%
    dplyr::left_join(Linda.res %>% dplyr::transmute(feature = features, pval.x = -log10(het.qval)), by = "feature") %>%
    dplyr::filter(!is.na(pval.x))

  plot(data.plot$pval.x, data.plot$pval.y,
       ylab = "",
       xlab = "",
       col = "skyblue", pch = 18, mgp=c(2,0.5,0), yaxt = "n",
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0, 1, col = "brown", lty = 2)
  # mtext(expression(atop("Het. test " ~ -log10[10](p-val) ,"in LinDA")), side = 1, line = 4)

  ## PALM v.s. LM-CLR
  data.plot <- PALM.res %>% dplyr::transmute(feature = feature, pval.y = -log10(het.qval)) %>%
    dplyr::left_join(lmclr.res %>% dplyr::transmute(feature = features, pval.x = -log10(het.qval)), by = "feature") %>%
    dplyr::filter(!is.na(pval.x))

  plot(data.plot$pval.x, data.plot$pval.y, ylab = " ", xlab = "",
       col = "orange", pch = 18, mgp=c(2,0.5,0),  yaxt = "n",
       xlim =c(0, 4) ,
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0, 1, col = "brown", lty = 2)
  #mtext(expression(atop("Het. test " ~ -log10[10](p-val) ,"in LM-CLR")), side = 1, line = 4)

  ## PALM v.s. MaAsLin2
  data.plot <- PALM.res %>% dplyr::transmute(feature = feature, pval.y = -log10(het.qval)) %>%
    dplyr::left_join(Maaslin2.res %>% dplyr::transmute(feature = features, pval.x = -log10(het.qval)), by = "feature") %>%
    dplyr::filter(!is.na(pval.x))

  plot(data.plot$pval.x, data.plot$pval.y, ylab = " ", xlab = "",
       col = "#4dac26", pch = 18, mgp=c(2,0.5,0),  yaxt = "n",
       cex.axis = 1.6,  cex.lab = 1.6)
  abline(0, 1, col = "brown", lty = 2)
  #mtext(expression(atop("Het. test " ~ -log10[10](p-val) ,"in MaAsLin2")), side = 1, line = 4)

  dev.off()

  # Figure c ----
  rm(list = ls())

  load("./CRC/Processed_data/data.org.K401.Rdata")
  load("./CRC/Output/CRC_output.Rdata")
  load("./Data/CRC_data/data/tax.Rdata")
  target.fdr <- 0.05
  target.pval.fdr <- 0.1

  taxa.id.sig <- intersect(
    intersect(
      intersect(
        intersect(PALM.res$feature[PALM.res$qval <= target.fdr & !is.na(PALM.res$qval)],
                  ANCOMBC2.res$features[ANCOMBC2.res$qval <= target.fdr & !is.na(ANCOMBC2.res$qval)]),
        Maaslin2.res$features[Maaslin2.res$qval <= target.fdr & !is.na(Maaslin2.res$qval)]),
      lmclr.res$features[lmclr.res$qval <= target.fdr & !is.na(lmclr.res$qval)]),
    Linda.res$features[Linda.res$qval <= target.fdr & !is.na(Linda.res$qval)])

  taxa.id.het <-  unique(c(ANCOMBC2.res$features[ANCOMBC2.res$het.qval <= target.pval.fdr & !is.na(ANCOMBC2.res$het.qval)],
                           PALM.res$feature[PALM.res$het.qval <= target.pval.fdr & !is.na(PALM.res$het.qval)],
                           Maaslin2.res$features[Maaslin2.res$het.qval <= target.pval.fdr & !is.na(Maaslin2.res$het.qval)],
                           lmclr.res$features[lmclr.res$het.qval <= target.pval.fdr & !is.na(lmclr.res$het.qval)],
                           Linda.res$features[Linda.res$het.qval <= target.pval.fdr & !is.na(Linda.res$het.qval)]))


  taxa.id.all <- intersect(taxa.id.sig, taxa.id.het)
  taxa.id.all <- taxa.id.all[order(PALM.res$coef[PALM.res$feature %in% taxa.id.all])]

  taxa.id.palm <- setdiff(PALM.res$feature[PALM.res$qval <= target.fdr & !is.na(PALM.res$qval)],
                          unique(c(ANCOMBC2.res$features[ANCOMBC2.res$qval <= target.fdr & !is.na(ANCOMBC2.res$qval)],
                                   Maaslin2.res$features[Maaslin2.res$qval <= target.fdr & !is.na(Maaslin2.res$qval)],
                                   lmclr.res$features[lmclr.res$qval <= target.fdr & !is.na(lmclr.res$qval)],
                                   Linda.res$features[Linda.res$qval <= target.fdr & !is.na(Linda.res$qval)])))

  taxa.id.palm <- taxa.id.palm[order(PALM.res$coef[PALM.res$feature %in% taxa.id.palm])]

  # Rename the taxa name
  rownames(tax) <- tax$OTUID
  renames <- sub("^\\S+\\s+",
                 '', tax[gsub("[][]", "",unlist(regmatches(taxa.id.all,
                                                           gregexpr("\\[.*?\\]",taxa.id.all)))),"mOTU"])
  # 1. "sp.." to "species" and remove following part
  renames <- gsub(" sp..*"," species",renames)
  renames <- gsub(" subsp..*", " subspecies", renames)
  renames <- gsub(" 3_1_57FAA_CT1", "", renames)
  # 2. remove [C] and remove []
  renames <- gsub("\\[C\\]","",renames)
  renames <- gsub("\\[|\\]","",renames)
  # 3. "gen." to "genus"
  renames <- gsub("gen..*","genus",renames)
  # 4. shorten Fusobacterium
  renames <- gsub("Fusobacterium","F.",renames)
  # Upper unknown
  renames <- str_trim(sub("unknown", "Unknown", renames))
  sp_strs <- str_split(renames, " ")
  for(l in 1:length(sp_strs)){
    strs <- sp_strs[[l]]
    for(ll in 1:length(strs)){
      if(strs[[ll]] %in% c("Unknown", "species", "genus", "subspecies")){
        strs[[ll]] <- paste0("</i>",strs[[ll]]," <i>")
      }
    }
    strs[[1]] <- paste0(" <i>",strs[[1]])
    strs[[length(strs)]] <- paste0(strs[[length(strs)]],"</i>")
    sp_strs[[l]] <- strs
  }
  renames <- sub("<i></i>", "", unlist(lapply(sp_strs, paste0, collapse=" ")))
  renames <- paste0(renames, " (",sub("]",")", sub(".*v2_", "", taxa.id.all),")"))
  renames[renames == " <i>Clostridium</i> (0860)"] <- " <i>Clostridium</i> species (0860)"
  names(renames) <- taxa.id.all

  ## ANCOM-BC2
  colnames(ANCOMBC2.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(ANCOMBC2.res[,"het.pval"], method = "fdr")
  names(qvals) <- ANCOMBC2.res$feature
  qvals[is.na(qvals)] <- 1
  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.all){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(ANCOMBC2.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = ANCOMBC2.model$est[l,],
                                     AA.lower = ANCOMBC2.model$est[l,] - sqrt(ANCOMBC2.model$var[l,]),
                                     AA.upper = ANCOMBC2.model$est[l,] + sqrt(ANCOMBC2.model$var[l,])))

    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       AA = ANCOMBC2.res[l,"est"]))
    }
  }

  g.ANCOMBC2.c <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("ANCOM-BC2\nCoefficients") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() + theme_minimal() + ylim(-1.6, 3.0) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_markdown(size = 15),
          axis.text.x =  element_blank(),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.ANCOMBC2.c <- g.ANCOMBC2.c + geom_rect(data = df.area,
                                             aes(xmin = as.numeric(species) - 0.5,
                                                 xmax = as.numeric(species) + 0.5,
                                                 ymin = -Inf, ymax = Inf),
                                             fill = "yellow", alpha = 0.2)
  }

  ## LinDA
  colnames(Linda.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(Linda.res[,"het.pval"], method = "fdr")
  names(qvals) <- Linda.res$feature
  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.all){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(Linda.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = Linda.model$est[l,],
                                     AA.lower = Linda.model$est[l,] - sqrt(Linda.model$var[l,]),
                                     AA.upper = Linda.model$est[l,] + sqrt(Linda.model$var[l,])))
    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       AA = Linda.res[l,"est"]))
    }
  }

  g.LinDA.c <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("LinDA\nCoefficients") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() + geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    ylim(-2.6, 2.4) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x =  element_blank(),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.LinDA.c <- g.LinDA.c + geom_rect(data = df.area,
                                       aes(xmin = as.numeric(species) - 0.5,
                                           xmax = as.numeric(species) + 0.5,
                                           ymin = -Inf, ymax = Inf),
                                       fill = "yellow", alpha = 0.2)
  }

  ## LM-CLR
  colnames(lmclr.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(lmclr.res[,"het.pval"], method = "fdr")
  names(qvals) <- lmclr.res$feature

  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.all){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(lmclr.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = lmclr.model$est[l,],
                                     AA.lower = lmclr.model$est[l,] - sqrt(lmclr.model$var[l,]),
                                     AA.upper = lmclr.model$est[l,] + sqrt(lmclr.model$var[l,])))
    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       AA = lmclr.res[l,"est"]))
    }
  }

  g.lmclr.c <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("LM-CLR\nCoefficients") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() +  geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    ylim(-1.8, 1.9) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.lmclr.c <- g.lmclr.c + geom_rect(data = df.area,
                                       aes(xmin = as.numeric(species) - 0.5,
                                           xmax = as.numeric(species) + 0.5,
                                           ymin = -Inf, ymax = Inf),
                                       fill = "yellow", alpha = 0.2)
  }

  ## Maaslin2
  colnames(Maaslin2.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(Maaslin2.res[,"het.pval"], method = "fdr")
  names(qvals) <- Maaslin2.res$feature

  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.all){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(Maaslin2.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = Maaslin2.model$est[l,],
                                     AA.lower = Maaslin2.model$est[l,] - sqrt(Maaslin2.model$var[l,]),
                                     AA.upper = Maaslin2.model$est[l,] + sqrt(Maaslin2.model$var[l,])))
    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       qval = qvals[l],
                                       AA = Maaslin2.res[l,"est"]))
    }
  }

  g.Maaslin2.c <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("MaAsLin2\nCoefficients") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() + geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    ylim(-0.05, 0.05) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.Maaslin2.c <- g.Maaslin2.c + geom_rect(data = df.area,
                                             aes(xmin = as.numeric(species) - 0.5,
                                                 xmax = as.numeric(species) + 0.5,
                                                 ymin = -Inf, ymax = Inf),
                                             fill = "yellow", alpha = 0.2)
  }

  ## PALM
  colnames(PALM.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(PALM.res[,"het.pval"], method = "fdr")
  names(qvals) <- PALM.res$feature

  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.all){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(PALM.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = PALM.model$est[l,],
                                     AA.lower = PALM.model$est[l,] - sqrt(PALM.model$var[l,]),
                                     AA.upper = PALM.model$est[l,] + sqrt(PALM.model$var[l,])))
    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       AA = PALM.res[l,"coef"]))
    }
  }

  g.PALM.c <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("PALM\nCoefficients") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    ylim(-3, 3.8) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x =  element_blank(),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = c(0.2, 0.85),
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.PALM.c <- g.PALM.c + geom_rect(data = df.area,
                                     aes(xmin = as.numeric(species) - 0.5,
                                         xmax = as.numeric(species) + 0.5,
                                         ymin = -Inf, ymax = Inf),
                                     fill = "yellow", alpha = 0.2)
  }

  ## Unique pattern for PALM
  # Rename the taxa name
  rownames(tax) <- tax$OTUID
  renames <- sub("^\\S+\\s+",
                 '', tax[gsub("[][]", "",unlist(regmatches(taxa.id.palm,
                                                           gregexpr("\\[.*?\\]",taxa.id.palm)))),"mOTU"])
  # 1. "sp.." to "species" and remove following part
  renames <- gsub(" sp..*"," species",renames)
  renames <- gsub(" subsp..*", " subspecies", renames)
  renames <- gsub(" 3_1_57FAA_CT1", "", renames)
  # 2. remove [C] and remove []
  renames <- gsub("\\[C\\]","",renames)
  renames <- gsub("\\[|\\]","",renames)
  # 3. "gen." to "genus"
  renames <- gsub("gen..*","genus",renames)
  # 4. shorten Fusobacterium
  renames <- gsub("Fusobacterium","F.",renames)
  # Upper unknown
  renames <- str_trim(sub("unknown", "Unknown", renames))
  sp_strs <- str_split(renames, " ")
  for(l in 1:length(sp_strs)){
    strs <- sp_strs[[l]]
    for(ll in 1:length(strs)){
      if(strs[[ll]] %in% c("Unknown", "species", "genus", "subspecies")){
        strs[[ll]] <- paste0("</i>",strs[[ll]]," <i>")
      }
    }
    strs[[1]] <- paste0(" <i>",strs[[1]])
    strs[[length(strs)]] <- paste0(strs[[length(strs)]],"</i>")
    sp_strs[[l]] <- strs
  }
  renames <- sub("<i></i>", "", unlist(lapply(sp_strs, paste0, collapse=" ")))
  renames <- paste0(renames, " (",sub("]",")", sub(".*v2_", "", taxa.id.palm),")"))
  renames[renames == " <i>Clostridium</i> (0860)"] <- " <i>Clostridium</i> species (0860)"
  names(renames) <- taxa.id.palm

  ## ANCOM-BC2
  colnames(ANCOMBC2.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(ANCOMBC2.res[,"het.pval"], method = "fdr")
  names(qvals) <- ANCOMBC2.res$feature
  qvals[is.na(qvals)] <- 1
  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.palm){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(ANCOMBC2.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = ANCOMBC2.model$est[l,],
                                     AA.lower = ANCOMBC2.model$est[l,] - sqrt(ANCOMBC2.model$var[l,]),
                                     AA.upper = ANCOMBC2.model$est[l,] + sqrt(ANCOMBC2.model$var[l,])))

    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       AA = ANCOMBC2.res[l,"est"]))
    }
  }

  g.ANCOMBC2.d <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("ANCOM-BC2\nCoefficients") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() + theme_minimal() + ylim(-1.6, 3.0) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_markdown(size = 15),
          axis.text.x =  element_text(size = 15),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.ANCOMBC2.d <- g.ANCOMBC2.d + geom_rect(data = df.area,
                                             aes(xmin = as.numeric(species) - 0.5,
                                                 xmax = as.numeric(species) + 0.5,
                                                 ymin = -Inf, ymax = Inf),
                                             fill = "yellow", alpha = 0.2)
  }

  ## LinDA
  colnames(Linda.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(Linda.res[,"het.pval"], method = "fdr")
  names(qvals) <- Linda.res$feature
  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.palm){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(Linda.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = Linda.model$est[l,],
                                     AA.lower = Linda.model$est[l,] - sqrt(Linda.model$var[l,]),
                                     AA.upper = Linda.model$est[l,] + sqrt(Linda.model$var[l,])))
    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       AA = Linda.res[l,"est"]))
    }
  }

  g.LinDA.d <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("LinDA\nCoefficients") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() + geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    ylim(-2.6, 2.4) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x =  element_text(size = 15),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.LinDA.d <- g.LinDA.d + geom_rect(data = df.area,
                                       aes(xmin = as.numeric(species) - 0.5,
                                           xmax = as.numeric(species) + 0.5,
                                           ymin = -Inf, ymax = Inf),
                                       fill = "yellow", alpha = 0.2)
  }

  ## LM-CLR
  colnames(lmclr.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(lmclr.res[,"het.pval"], method = "fdr")
  names(qvals) <- lmclr.res$feature

  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.palm){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(lmclr.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = lmclr.model$est[l,],
                                     AA.lower = lmclr.model$est[l,] - sqrt(lmclr.model$var[l,]),
                                     AA.upper = lmclr.model$est[l,] + sqrt(lmclr.model$var[l,])))
    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       AA = lmclr.res[l,"est"]))
    }
  }

  g.lmclr.d <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("LM-CLR\nCoefficients") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() +  geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    ylim(-1.8, 1.9) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x =  element_text(size = 15),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.lmclr.d <- g.lmclr.d + geom_rect(data = df.area,
                                       aes(xmin = as.numeric(species) - 0.5,
                                           xmax = as.numeric(species) + 0.5,
                                           ymin = -Inf, ymax = Inf),
                                       fill = "yellow", alpha = 0.2)
  }

  ## Maaslin2
  colnames(Maaslin2.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(Maaslin2.res[,"het.pval"], method = "fdr")
  names(qvals) <- Maaslin2.res$feature

  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.palm){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(Maaslin2.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = Maaslin2.model$est[l,],
                                     AA.lower = Maaslin2.model$est[l,] - sqrt(Maaslin2.model$var[l,]),
                                     AA.upper = Maaslin2.model$est[l,] + sqrt(Maaslin2.model$var[l,])))
    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       qval = qvals[l],
                                       AA = Maaslin2.res[l,"est"]))
    }
  }

  g.Maaslin2.d <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("MaAsLin2\nCoefficients") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() + geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    ylim(-0.05, 0.05) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x =  element_text(size = 15),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.Maaslin2.d <- g.Maaslin2.d + geom_rect(data = df.area,
                                             aes(xmin = as.numeric(species) - 0.5,
                                                 xmax = as.numeric(species) + 0.5,
                                                 ymin = -Inf, ymax = Inf),
                                             fill = "yellow", alpha = 0.2)
  }

  ## PALM
  colnames(PALM.model$est) <- paste0("CRC", 1:5)
  qvals <- p.adjust(PALM.res[,"het.pval"], method = "fdr")
  names(qvals) <- PALM.res$feature

  df.plot <- NULL
  df.area <- NULL
  for(l in taxa.id.palm){
    df.plot <- rbind(df.plot, tibble(species = factor(renames[l], levels= renames),
                                     Study = factor(colnames(PALM.model$est),
                                                    levels = paste0("CRC", 5:1)),
                                     AA = PALM.model$est[l,],
                                     AA.lower = PALM.model$est[l,] - sqrt(PALM.model$var[l,]),
                                     AA.upper = PALM.model$est[l,] + sqrt(PALM.model$var[l,])))
    if(qvals[l] <= target.pval.fdr){
      df.area <- rbind(df.area, tibble(species = factor(renames[l], levels= renames),
                                       AA = PALM.res[l,"coef"]))
    }
  }

  g.PALM.d <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("PALM\nCoefficients") +
    scale_color_manual(values = c("#F8766D","#A3A500","#00BF7D","#00B0F6","#E76BF3"),
                       breaks = c("CRC5", "CRC4", "CRC3", "CRC2", "CRC1")) +
    coord_flip() +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    ylim(-3, 3.8) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x =  element_text(size = 15),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.PALM.d <- g.PALM.d + geom_rect(data = df.area,
                                     aes(xmin = as.numeric(species) - 0.5,
                                         xmax = as.numeric(species) + 0.5,
                                         ymin = -Inf, ymax = Inf),
                                     fill = "yellow", alpha = 0.2)
  }

  ## plots
  pdf("./Figure/Figure4_c.pdf", width = 20.75, height = 10, bg = "white")

  plot_grid(plot_grid(g.ANCOMBC2.c, g.ANCOMBC2.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
            plot_grid(g.LinDA.c, g.LinDA.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
            plot_grid(g.lmclr.c, g.lmclr.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
            plot_grid(g.Maaslin2.c, g.Maaslin2.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
            plot_grid(g.PALM.c, g.PALM.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
            ncol = 5, align = 'h', rel_widths = c(0.2, 0.1, 0.1, 0.1, 0.1))


  dev.off()
