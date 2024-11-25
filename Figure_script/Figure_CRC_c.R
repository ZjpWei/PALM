# =============================================== #
#       CRC Data: Figure b, scatter plots         #
# =============================================== #

  # Packages ----
  library("plotly")
  library("tidyr")
  library("stringr")
  library("dplyr")
  library("ggplot2")
  library("cluster")
  library("ggdendro")
  library("grid")
  library("randomcoloR")
  library("ggtext")
  library("RColorBrewer")
  library("cowplot")

  # General ----
  rm(list = ls())
  load("./CRC/Processed_data/data.org.K401.Rdata")
  load("./CRC/Output/CRC_output.Rdata")

  target.fdr <- 0.05
  target.pval.fdr <- 0.1

  taxa.id.all <- intersect(
                 intersect(
                 intersect(
                 intersect(PALM.res$feature[PALM.res$qval <= target.fdr],
                           ANCOMBC2.res$features[ANCOMBC2.res$het.qval <= target.fdr & !is.na(ANCOMBC2.res$het.qval)]),
                           Maaslin2.res$features[Maaslin2.res$qval <= target.fdr]),
                           lmclr.res$features[lmclr.res$qval <= target.fdr]),
                           Linda.res$features[Linda.res$qval <= target.fdr])

  taxa.id.all <- taxa.id.all[order(PALM.res$coef[PALM.res$feature %in% taxa.id.all])]

  ## taxa rename
  renames <- gsub(" sp..*"," species", taxa.id.all)
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
  sp_strs <- lapply(sp_strs, function(d){d[1:(length(d)-1)]})
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

  g.ANCOMBC2 <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("ANCOM-BC2\nCoefficients") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    coord_flip() + theme_minimal() +
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
    g.ANCOMBC2 <- g.ANCOMBC2 + geom_rect(data = df.area,
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

  g.LinDA <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("LinDA\nCoefficients") +
    coord_flip() +     geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
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
    g.LinDA <- g.LinDA + geom_rect(data = df.area,
                                         aes(xmin = as.numeric(species) - 0.5,
                                             xmax = as.numeric(species) + 0.5,
                                             ymin = -Inf, ymax = Inf),
                                         fill = "yellow", alpha = 0.2)
  }

  ## LM-INT
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

  g.lmclr <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("LM-CLR\nCoefficients") +
    coord_flip() +  geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
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
    g.lmclr <- g.lmclr + geom_rect(data = df.area,
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
                                       AA = Maaslin2.res[l,"est"]))
    }
  }

  g.Maaslin2 <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("MaAsLin2\nCoefficients") +
    coord_flip() + geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
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
    g.Maaslin2 <- g.Maaslin2 + geom_rect(data = df.area,
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

  g.PALM <- df.plot %>%
    ggplot(aes(x=species, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    theme_minimal() + ylab("PALM\nCoefficients") +
    coord_flip() +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x =  element_text(size = 15),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = c(0.2, 0.85),
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(!is.null(df.area)){
    g.PALM <- g.PALM + geom_rect(data = df.area,
                                   aes(xmin = as.numeric(species) - 0.5,
                                       xmax = as.numeric(species) + 0.5,
                                       ymin = -Inf, ymax = Inf),
                                   fill = "yellow", alpha = 0.2)
  }

  ## plots
  pdf("./Figure/Figure_rel_c.pdf", width = 19.75, height = 7, bg = "white")

  plot_grid(g.ANCOMBC2, g.LinDA, g.lmclr, g.Maaslin2, g.PALM,
            nrow = 1, align = 'h', rel_widths = c(0.2, 0.10, 0.10, 0.10, 0.10))

  dev.off()

