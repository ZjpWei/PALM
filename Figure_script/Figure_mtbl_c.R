# =============================================== #
#       Metabolites meta-analysis: Figure c       #
#                  HMDB0000039                    #
# =============================================== #

  # Package ----
  library(plotly)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(cluster)
  library(ggdendro)
  library(grid)
  library(randomcoloR)
  library(ggtext)
  library(RColorBrewer)
  library(cowplot)

  rm(list = ls())
  # Selection ----
  selected.lst <- "HMDB0000039"
  load("./Metabolites/Output/MTBL_PALM.Rdata")
  load(paste0("./Metabolites/Output/Metabolites/Compare_method_", selected.lst, ".Rdata"))
  study.id <- which(apply(ANCOMBC2.est, 2, function(d){return(!all(is.na(d)))}))
  study.names <- names(study.id)
  target.fdr <- 0.05
  target.get.fdr <- 0.1

  ## PALM
  PALM.df <- PALM.meta.coef %>% dplyr::select(feature, coef = selected.lst) %>%
    dplyr::left_join(
      PALM.meta.pval %>% dplyr::select(feature, pval = selected.lst) %>%
        dplyr::mutate(qval = p.adjust(pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      PALM.meta.pval.het %>% dplyr::select(feature, het.pval = selected.lst) %>%
        dplyr::mutate(het.qval = p.adjust(het.pval, method = "fdr")), by = "feature")
  for(st.id in names(PALM.test[[selected.lst]]$palm_fits)){
    PALM.df <- PALM.df %>% dplyr::left_join(
      PALM.test[[selected.lst]]$palm_fits[[st.id]] %>%
        dplyr::select(feature, !!paste0(st.id,":est") := coef, !!paste0(st.id,":std") := stderr), by = "feature")
  }
  PALM.df <- PALM.df %>% dplyr::mutate(feature = gsub(".*;g__","",feature)) %>%
    tibble::column_to_rownames("feature")

  ## ANCOM-BC2
  ANCOMBC2.df <- ANCOMBC2.meta.coef %>% dplyr::select(feature, coef = selected.lst) %>%
    dplyr::left_join(
      ANCOMBC2.meta.pval %>% dplyr::select(feature, pval = selected.lst) %>%
        dplyr::mutate(qval = p.adjust(pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      ANCOMBC2.meta.pval.het %>% dplyr::select(feature, het.pval = selected.lst) %>%
        dplyr::mutate(het.qval = p.adjust(het.pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      data.frame(ANCOMBC2.est[,study.id]) %>%
        dplyr::rename_with(~ paste0(., ":est")) %>% tibble::rownames_to_column("feature")
    ) %>% dplyr::left_join(
      data.frame(sqrt(ANCOMBC2.var)[,apply(ANCOMBC2.var, 2, function(d){return(!all(is.na(d)))})]) %>%
        dplyr::rename_with(~ paste0(., ":std")) %>% tibble::rownames_to_column("feature")
    ) %>% dplyr::mutate(feature = gsub(".*;g__","",feature)) %>% tibble::column_to_rownames("feature")

  ## MaAsLin2
  Maaslin2.df <- Maaslin2.meta.coef %>% dplyr::select(feature, coef = selected.lst) %>%
    dplyr::left_join(
      Maaslin2.meta.pval %>% dplyr::select(feature, pval = selected.lst) %>%
        dplyr::mutate(qval = p.adjust(pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      Maaslin2.meta.pval.het %>% dplyr::select(feature, het.pval = selected.lst) %>%
        dplyr::mutate(het.qval = p.adjust(het.pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      data.frame(Maaslin2.est[,study.id]) %>%
        dplyr::rename_with(~ paste0(., ":est")) %>% tibble::rownames_to_column("feature")
    ) %>% dplyr::left_join(
      data.frame(sqrt(Maaslin2.var)[,apply(Maaslin2.var, 2, function(d){return(!all(is.na(d)))})]) %>%
        dplyr::rename_with(~ paste0(., ":std")) %>% tibble::rownames_to_column("feature")
    ) %>% dplyr::mutate(feature = gsub(".*;g__","",feature)) %>% tibble::column_to_rownames("feature")

  ## LM-CLR
  LMCLR.df <- LMCLR.meta.coef %>% dplyr::select(feature, coef = selected.lst) %>%
    dplyr::left_join(
      LMCLR.meta.pval %>% dplyr::select(feature, pval = selected.lst) %>%
        dplyr::mutate(qval = p.adjust(pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      LMCLR.meta.pval.het %>% dplyr::select(feature, het.pval = selected.lst) %>%
        dplyr::mutate(het.qval = p.adjust(het.pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      data.frame(LMCLR.est[,study.id]) %>%
        dplyr::rename_with(~ paste0(., ":est")) %>% tibble::rownames_to_column("feature")
    ) %>% dplyr::left_join(
      data.frame(sqrt(LMCLR.var)[,apply(LMCLR.var, 2, function(d){return(!all(is.na(d)))})]) %>%
        dplyr::rename_with(~ paste0(., ":std")) %>% tibble::rownames_to_column("feature")
    ) %>% dplyr::mutate(feature = gsub(".*;g__","",feature)) %>% tibble::column_to_rownames("feature")

  ## LinDA
  Linda.df <- Linda.meta.coef %>% dplyr::select(feature, coef = selected.lst) %>%
    dplyr::left_join(
      Linda.meta.pval %>% dplyr::select(feature, pval = selected.lst) %>%
        dplyr::mutate(qval = p.adjust(pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      Linda.meta.pval.het %>% dplyr::select(feature, het.pval = selected.lst) %>%
        dplyr::mutate(het.qval = p.adjust(het.pval, method = "fdr")), by = "feature") %>%
    dplyr::left_join(
      data.frame(LinDA.est[,study.id]) %>%
        dplyr::rename_with(~ paste0(., ":est")) %>% tibble::rownames_to_column("feature")
    ) %>% dplyr::left_join(
      data.frame(sqrt(LinDA.var)[,apply(LinDA.var, 2, function(d){return(!all(is.na(d)))})]) %>%
        dplyr::rename_with(~ paste0(., ":std")) %>% tibble::rownames_to_column("feature")
    ) %>% dplyr::mutate(feature = gsub(".*;g__","",feature)) %>% tibble::column_to_rownames("feature")

  ##@ Metabolite information
  source("./utility/hmdb_utils.R")
  hmdb.ids <- selected.lst
  hmdb.data <- get.hmdb.data.by.ids(hmdb.ids)
  hmdb.data <- hmdb.data %>% dplyr::rename(Compound = HMDB)

  genus.lst <- intersect(rownames(PALM.df)[PALM.df$qval <= target.fdr & ANCOMBC2.df$qval <= target.fdr &
                                 Maaslin2.df$qval <= target.fdr & Linda.df$qval <= target.fdr & LMCLR.df$qval <= target.fdr],

                         rownames(PALM.df)[PALM.df$het.qval <= target.get.fdr | ANCOMBC2.df$het.qval <= target.get.fdr |
                                           Maaslin2.df$het.qval <= target.get.fdr | Linda.df$het.qval <= target.get.fdr |
                                           LMCLR.df$het.qval <= target.get.fdr])

  genus.lst <- genus.lst[order(PALM.df[genus.lst, "coef"])]

  # PALM plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     Study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = PALM.df[genus.lst, paste0(l,":est")],
                                     AA.lower = PALM.df[genus.lst, paste0(l,":est")] - PALM.df[genus.lst, paste0(l,":std")],
                                     AA.upper = PALM.df[genus.lst, paste0(l,":est")] + PALM.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[PALM.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)

  g.PALM.c <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = Study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = Study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("PALM\nCoefficients") + xlab("Genus") + ylim(-2.2, 1) +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_blank(),
          axis.text.x =  element_blank(),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = c(0.25, 0.8),
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(nrow(df.area) != 0){
    g.PALM.c <- g.PALM.c + geom_rect(data = df.area,
                                 aes(xmin = as.numeric(genus) - 0.5,
                                     xmax = as.numeric(genus) + 0.5,
                                     ymin = -Inf, ymax = Inf),
                                 fill = "yellow", alpha = 0.2)
  }

  # MaAsLin2 plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = Maaslin2.df[genus.lst, paste0(l,":est")],
                                     AA.lower = Maaslin2.df[genus.lst, paste0(l,":est")] - Maaslin2.df[genus.lst, paste0(l,":std")],
                                     AA.upper = Maaslin2.df[genus.lst, paste0(l,":est")] + Maaslin2.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[Maaslin2.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)

  g.MaAsLin2.c <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("MaAsLin2\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-0.02, 0.04) +
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

  if(nrow(df.area) != 0){
    g.MaAsLin2.c <- g.MaAsLin2.c + geom_rect(data = df.area,
                                 aes(xmin = as.numeric(genus) - 0.5,
                                     xmax = as.numeric(genus) + 0.5,
                                     ymin = -Inf, ymax = Inf),
                                 fill = "yellow", alpha = 0.2)
  }

  # ANCOM-BC2 ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = ANCOMBC2.df[genus.lst, paste0(l,":est")],
                                     AA.lower = ANCOMBC2.df[genus.lst, paste0(l,":est")] - ANCOMBC2.df[genus.lst, paste0(l,":std")],
                                     AA.upper = ANCOMBC2.df[genus.lst, paste0(l,":est")] + ANCOMBC2.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[ANCOMBC2.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)
  g.ANCOMBC2.c <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("ANCOM-BC2\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-0.8, 0.7) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black', linewidth = 1),
          axis.text.y = element_markdown(size = 15),
          axis.text.x = element_blank(),

          legend.title = element_text(hjust = 0.5, size = 16),
          legend.text = element_text(size = 13),
          legend.position = "none",
          legend.direction = "vertical",
          legend.box = "vertical",
          strip.text = element_blank()) +
    guides(color = guide_legend(reverse=T))

  if(nrow(df.area) != 0){
    g.ANCOMBC2.c <- g.ANCOMBC2.c + geom_rect(data = df.area,
                                         aes(xmin = as.numeric(genus) - 0.5,
                                             xmax = as.numeric(genus) + 0.5,
                                             ymin = -Inf, ymax = Inf),
                                         fill = "yellow", alpha = 0.2)
  }

  # LinDA plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = Linda.df[genus.lst, paste0(l,":est")],
                                     AA.lower = Linda.df[genus.lst, paste0(l,":est")] - Linda.df[genus.lst, paste0(l,":std")],
                                     AA.upper = Linda.df[genus.lst, paste0(l,":est")] + Linda.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[Linda.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)
  g.Linda.c <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("LinDA\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-1.5, 1) +
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

  if(nrow(df.area) != 0){
    g.Linda.c <- g.Linda.c + geom_rect(data = df.area,
                                         aes(xmin = as.numeric(genus) - 0.5,
                                             xmax = as.numeric(genus) + 0.5,
                                             ymin = -Inf, ymax = Inf),
                                         fill = "yellow", alpha = 0.2)
  }

  # LM-CLR plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = LMCLR.df[genus.lst, paste0(l,":est")],
                                     AA.lower = LMCLR.df[genus.lst, paste0(l,":est")] - LMCLR.df[genus.lst, paste0(l,":std")],
                                     AA.upper = LMCLR.df[genus.lst, paste0(l,":est")] + LMCLR.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[LMCLR.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)
  g.LMCLR.c <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("LM-CLR\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-1, 0.8) +
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

  if(nrow(df.area) != 0){
    g.LMCLR.c <- g.LMCLR.c + geom_rect(data = df.area,
                                   aes(xmin = as.numeric(genus) - 0.5,
                                       xmax = as.numeric(genus) + 0.5,
                                       ymin = -Inf, ymax = Inf),
                                   fill = "yellow", alpha = 0.2)
  }


  ## Figure d
  genus.lst <- c("Lachnospira", "Faecalibacterium", "Eubacterium_I")
  genus.lst <- genus.lst[order(PALM.df[genus.lst, "coef"])]

  # PALM plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = PALM.df[genus.lst, paste0(l,":est")],
                                     AA.lower = PALM.df[genus.lst, paste0(l,":est")] - PALM.df[genus.lst, paste0(l,":std")],
                                     AA.upper = PALM.df[genus.lst, paste0(l,":est")] + PALM.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[PALM.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)

  g.PALM.d <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("PALM\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-2.2, 1) +
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

  if(nrow(df.area) != 0){
    g.PALM.d <- g.PALM.d + geom_rect(data = df.area,
                                 aes(xmin = as.numeric(genus) - 0.5,
                                     xmax = as.numeric(genus) + 0.5,
                                     ymin = -Inf, ymax = Inf),
                                 fill = "yellow", alpha = 0.2)
  }


  # MaAsLin2 plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = Maaslin2.df[genus.lst, paste0(l,":est")],
                                     AA.lower = Maaslin2.df[genus.lst, paste0(l,":est")] - Maaslin2.df[genus.lst, paste0(l,":std")],
                                     AA.upper = Maaslin2.df[genus.lst, paste0(l,":est")] + Maaslin2.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[Maaslin2.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)

  g.MaAsLin2.d <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("MaAsLin2\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-0.02, 0.04) +
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

  if(nrow(df.area) != 0){
    g.MaAsLin2.d <- g.MaAsLin2.d + geom_rect(data = df.area,
                                         aes(xmin = as.numeric(genus) - 0.5,
                                             xmax = as.numeric(genus) + 0.5,
                                             ymin = -Inf, ymax = Inf),
                                         fill = "yellow", alpha = 0.2)
  }

  # ANCOM-BC2 plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = ANCOMBC2.df[genus.lst, paste0(l,":est")],
                                     AA.lower = ANCOMBC2.df[genus.lst, paste0(l,":est")] - ANCOMBC2.df[genus.lst, paste0(l,":std")],
                                     AA.upper = ANCOMBC2.df[genus.lst, paste0(l,":est")] + ANCOMBC2.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[ANCOMBC2.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)
  g.ANCOMBC2.d <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("ANCOM-BC2\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-0.8, 0.7) +
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_blank(),
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

  if(nrow(df.area) != 0){
    g.ANCOMBC2.d <- g.ANCOMBC2.d + geom_rect(data = df.area,
                                         aes(xmin = as.numeric(genus) - 0.5,
                                             xmax = as.numeric(genus) + 0.5,
                                             ymin = -Inf, ymax = Inf),
                                         fill = "yellow", alpha = 0.2)
  }

  # LinDA plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = Linda.df[genus.lst, paste0(l,":est")],
                                     AA.lower = Linda.df[genus.lst, paste0(l,":est")] - Linda.df[genus.lst, paste0(l,":std")],
                                     AA.upper = Linda.df[genus.lst, paste0(l,":est")] + Linda.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[Linda.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)
  g.Linda.d <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("LinDA\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-1.5, 1) +
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

  if(nrow(df.area) != 0){
    g.Linda.d <- g.Linda.d + geom_rect(data = df.area,
                                   aes(xmin = as.numeric(genus) - 0.5,
                                       xmax = as.numeric(genus) + 0.5,
                                       ymin = -Inf, ymax = Inf),
                                   fill = "yellow", alpha = 0.2)
  }

  # LM-CLR plots ----
  df.plot <- NULL
  for(l in study.names){
    df.plot <- rbind(df.plot, tibble(genus = factor(genus.lst, levels= genus.lst),
                                     study = factor(rep(l, length(genus.lst)), levels = study.names,
                                                    labels = paste0("MTBL", rev(study.id))),
                                     AA = LMCLR.df[genus.lst, paste0(l,":est")],
                                     AA.lower = LMCLR.df[genus.lst, paste0(l,":est")] - LMCLR.df[genus.lst, paste0(l,":std")],
                                     AA.upper = LMCLR.df[genus.lst, paste0(l,":est")] + LMCLR.df[genus.lst, paste0(l,":std")]))
  }
  df.area <- tibble(genus = factor(genus.lst[LMCLR.df[genus.lst, "het.qval"] <= target.get.fdr],
                                   levels= genus.lst), AA = Inf)
  g.LMCLR.d <- df.plot %>%
    ggplot(aes(x=genus, y=AA)) +
    geom_errorbar(aes(ymin = AA.lower, ymax = AA.upper, group = study),
                  position = position_dodge(width = 0.6),
                  width = 0, linewidth = 0.5) +
    geom_point(pch = 18, aes(color = study),
               size = 2.5, position = position_dodge(width = 0.6)) +
    ylab("LM-CLR\nCoefficients") + xlab("Genus") +
    geom_hline(aes(yintercept = 0),colour="#990000", linetype="dashed") +
    scale_color_manual(values = c("#F8766D","#CD9600","#7CAE00","#00BE67",
                                  "#00BFC4","#00A9FF","#C77CFF","#FF61CC"),
                       breaks = c("MTBL8", "MTBL7", "MTBL6", "MTBL5",
                                  "MTBL4", "MTBL3", "MTBL2", "MTBL1")) +
    coord_flip() + theme_minimal() + ylim(-1, 0.8) +
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

  if(nrow(df.area) != 0){
    g.LMCLR.d <- g.LMCLR.d + geom_rect(data = df.area,
                                   aes(xmin = as.numeric(genus) - 0.5,
                                       xmax = as.numeric(genus) + 0.5,
                                       ymin = -Inf, ymax = Inf),
                                   fill = "yellow", alpha = 0.2)
  }

  pp_d <- plot_grid(plot_grid(g.ANCOMBC2.c, g.ANCOMBC2.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
                    plot_grid(g.Linda.c, g.Linda.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
                    plot_grid(g.LMCLR.c, g.LMCLR.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
                    plot_grid(g.MaAsLin2.c, g.MaAsLin2.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
                    plot_grid(g.PALM.c, g.PALM.d, nrow = 2, align = 'v', rel_heights = c(0.7, 0.3)),
                    ncol = 5, align = 'h', rel_widths = c(0.2, 0.12, 0.12, 0.12, 0.12))

  pdf("./Figure/Figure_mtbl_c.pdf", width = 17.60, height = 7.7, bg = "white")

  pp_d

  dev.off()

