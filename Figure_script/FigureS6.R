# =============================================== #
#                    FigureS6                     #
# =============================================== #

  # Packages ----
  library("plotly")
  library("tidyr")
  library("ggplot2")
  library("grid")
  library("randomcoloR")
  library("ggtext")

  # Figure a ----
  rm(list = ls())

  load("./Metabolites/Output/MTBL_PALM.Rdata")
  load("./Metabolites/Output/MTBL_others.Rdata")
  load("./Metabolites/Processed_data/MTBL.RData")

  feature.ids <- PALM.meta.coef$feature
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

  sample.id <- as.vector(unlist(lapply(data.for.lm, function(d){rownames(d$rel.abd)})))
  taxa <- NULL
  for(d in names(data.for.lm)){
    taxa <- c(taxa, colnames(data.for.lm[[d]]$rel.abd))
  }
  taxa <- unique(taxa)
  Otu_big_mat <- matrix(0, nrow = length(sample.id), ncol = length(taxa),
                        dimnames = list(sample.id, taxa))
  for(d in names(data.for.lm)){
    Otu_big_mat[rownames(data.for.lm[[d]]$rel.abd),colnames(data.for.lm[[d]]$rel.abd)] <- as.matrix(data.for.lm[[d]]$rel.abd)
  }
  tax.prevalence <- colMeans(Otu_big_mat / rowSums(Otu_big_mat))
  names(tax.prevalence) <- gsub(".*;g__","g__", names(tax.prevalence))

  ## ANCOMBC2
  df_ANCOMBC2 <- data.frame(PALM = rowSums(PALM_coef_mat!=0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(ANCOMBC2 = rowSums(ANCOMBC2_coef_mat!=0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Phocaeicola"] <- -1.5
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Alistipes"] <- 0
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Agathobacter"] <- 0.2
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Agathobacter"] <- -0.1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Bacteroides"] <- 0.6
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Bacteroides"] <- -0.2
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Faecalibacterium"] <- 0.3
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Faecalibacterium"] <- -0.7
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Bifidobacterium"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Bifidobacterium"] <- -0.2
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Prevotella"] <- -2
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Prevotella"] <- -0.5
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Lachnospira"] <- -3
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Lachnospira"] <- 1
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Collinsella"] <- -3
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Collinsella"] <- 2
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Mediterraneibacter"] <- -2
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Mediterraneibacter"] <- 4.5
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Parabacteroides"] <- 0
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Fusicatenibacter"] <- 0
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Parabacteroides"] <- 8
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Fusicatenibacter"] <- 8
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Gemmiger"] <- -1
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Ruminococcus_E"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Gemmiger"] <- 6
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Ruminococcus_E"] <- 6

  p_ANCOMBC2 <- df_ANCOMBC2 %>% ggplot(aes(x = ANCOMBC2, y = PALM, label = label)) +
      geom_jitter(alpha = 0.7, aes(size = prevalence, color = group)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in ANCOMBC2 results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(0, 255) + xlim(0, 255) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_ANCOMBC2$hjust, vjust=df_ANCOMBC2$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = c(0.85,0.2)) + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_a_ANCOMBC2.pdf", width = 8.15, height = 6.65, bg = "white")

  p_ANCOMBC2

  dev.off()

  ## LinDA ----
  df_Linda <- data.frame(PALM = rowSums(PALM_coef_mat!=0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(Linda = rowSums(Linda_coef_mat!=0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_Linda$vjust[df_Linda$label == "Phocaeicola"] <- -1.5
  df_Linda$vjust[df_Linda$label == "Bacteroides"] <- -1.5
  df_Linda$vjust[df_Linda$label == "Agathobacter"] <- -2
  df_Linda$vjust[df_Linda$label == "Faecalibacterium"] <- -2.2
  df_Linda$vjust[df_Linda$label == "Prevotella"] <- -4
  df_Linda$vjust[df_Linda$label == "Mediterraneibacter"] <- -1.5
  df_Linda$hjust[df_Linda$label == "Mediterraneibacter"] <- 1.3
  df_Linda$vjust[df_Linda$label == "Lachnospira"] <- 0
  df_Linda$hjust[df_Linda$label == "Lachnospira"] <- -0.2
  df_Linda$vjust[df_Linda$label == "Ruminococcus_E"] <- 3
  df_Linda$vjust[df_Linda$label == "Parabacteroides"] <- 9
  df_Linda$hjust[df_Linda$label == "Bifidobacterium"] <- 0.9

  p_Linda <- df_Linda %>% ggplot(aes(x = Linda, y = PALM, color = group, label = label)) +
      geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in Linda results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(0, 255) + xlim(0, 255) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_Linda$hjust, vjust=df_Linda$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_a_Linda.pdf", width = 8.15, height = 6.65, bg = "white")

  p_Linda

  dev.off()

  ## LM-CLR ----
  df_LMCLR <- data.frame(PALM = rowSums(PALM_coef_mat!=0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(LMCLR = rowSums(LMCLR_coef_mat!=0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_LMCLR$vjust[df_LMCLR$label == "Phocaeicola"] <- -1.5
  df_LMCLR$vjust[df_LMCLR$label == "Agathobacter"] <- -2
  df_LMCLR$vjust[df_LMCLR$label == "Bacteroides"] <- -1.5
  df_LMCLR$hjust[df_LMCLR$label == "Bacteroides"] <- 1
  df_LMCLR$hjust[df_LMCLR$label == "Bifidobacterium"] <- 1
  df_LMCLR$vjust[df_LMCLR$label == "Prevotella"] <- -2
  df_LMCLR$vjust[df_LMCLR$label == "Parabacteroides"] <- 9
  df_LMCLR$vjust[df_LMCLR$label == "Fusicatenibacter"] <- 9
  df_LMCLR$vjust[df_LMCLR$label == "Mediterraneibacter"] <- -1
  df_LMCLR$vjust[df_LMCLR$label == "Lachnospira"] <- 6
  df_LMCLR$vjust[df_LMCLR$label == "Ruminococcus_E"] <- 6

  p_LMCLR <- df_LMCLR %>% ggplot(aes(x = LMCLR, y = PALM, color = group, label = label)) +
      geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in LMCLR results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(0, 255) + xlim(0, 255) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_LMCLR$hjust, vjust=df_LMCLR$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_a_LMCLR.pdf", width = 8.15, height = 6.65, bg = "white")

  p_LMCLR

  dev.off()

  ## MaAsLin2 ----
  df_Maaslin2 <- data.frame(PALM = rowSums(PALM_coef_mat!=0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(Maaslin2 = rowSums(Maaslin2_coef_mat!=0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_Maaslin2$vjust[df_Maaslin2$label == "Phocaeicola"] <- -1.5
  df_Maaslin2$hjust[df_Maaslin2$label == "Agathobacter"] <- 1.2
  df_Maaslin2$vjust[df_Maaslin2$label == "Agathobacter"] <- 0
  df_Maaslin2$hjust[df_Maaslin2$label == "Bifidobacterium"] <- 0.2
  df_Maaslin2$vjust[df_Maaslin2$label == "Bifidobacterium"] <- -0.5
  df_Maaslin2$hjust[df_Maaslin2$label == "Lachnospira"] <- -1
  df_Maaslin2$vjust[df_Maaslin2$label == "Bacteroides"] <- -1.5
  df_Maaslin2$hjust[df_Maaslin2$label == "Bacteroides"] <- 0.2
  df_Maaslin2$vjust[df_Maaslin2$label == "Prevotella"] <- -2
  df_Maaslin2$vjust[df_Maaslin2$label == "Mediterraneibacter"] <- 6
  df_Maaslin2$hjust[df_Maaslin2$label == "Collinsella"] <- 2
  df_Maaslin2$hjust[df_Maaslin2$label == "Parabacteroides"] <- 0
  df_Maaslin2$vjust[df_Maaslin2$label == "Gemmiger"] <- 2
  df_Maaslin2$hjust[df_Maaslin2$label == "Faecalibacterium"] <- 0

  p_Maaslin2 <- df_Maaslin2 %>% ggplot(aes(x = Maaslin2, y = PALM, color = group, label = label)) +
      geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in Maaslin2 results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(0, 255) + xlim(0, 255) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_Maaslin2$hjust, vjust=df_Maaslin2$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_a_Maaslin2.pdf", width = 8.15, height = 6.65, bg = "white")

  p_Maaslin2

  dev.off()

  # Figure b ----
  rm(list = ls())

  load("./Metabolites/Output/MTBL_PALM.Rdata")
  load("./Metabolites/Output/MTBL_others.Rdata")
  load("./Metabolites/Processed_data/MTBL.RData")

  feature.ids <- PALM.meta.coef$feature
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

  sample.id <- as.vector(unlist(lapply(data.for.lm, function(d){rownames(d$rel.abd)})))
  taxa <- NULL
  for(d in names(data.for.lm)){
    taxa <- c(taxa, colnames(data.for.lm[[d]]$rel.abd))
  }
  taxa <- unique(taxa)
  Otu_big_mat <- matrix(0, nrow = length(sample.id), ncol = length(taxa),
                        dimnames = list(sample.id, taxa))
  for(d in names(data.for.lm)){
    Otu_big_mat[rownames(data.for.lm[[d]]$rel.abd),colnames(data.for.lm[[d]]$rel.abd)] <- as.matrix(data.for.lm[[d]]$rel.abd)
  }
  tax.prevalence <- colMeans(Otu_big_mat / rowSums(Otu_big_mat))
  names(tax.prevalence) <- gsub(".*;g__","g__", names(tax.prevalence))

  ## ANCOMBC2
  df_ANCOMBC2 <- data.frame(PALM = rowSums(PALM_coef_mat>0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(ANCOMBC2 = rowSums(ANCOMBC2_coef_mat>0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Phocaeicola"] <- -4
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Blautia_A"] <- -4
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Blautia_A"] <- 0
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Agathobacter"] <- 1
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Roseburia"] <- 0
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Bacteroides"] <- -3.5
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Prevotella"] <- 0
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Prevotella"] <- -4
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Bifidobacterium"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Bifidobacterium"] <- 1.5
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Collinsella"] <- -2
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Collinsella"] <- 3
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Lachnospira"] <- -2
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Lachnospira"] <- 3
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Parabacteroides"] <- -1.5
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Parabacteroides"] <- 1
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Mediterraneibacter"] <- 0.1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Mediterraneibacter"] <- 11
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Fusicatenibacter"] <- 0.3
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Fusicatenibacter"] <- 12
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Alistipes"] <- 6
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Alistipes"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Gemmiger"] <- 6
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Gemmiger"] <- -2
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Ruminococcus_E"] <- 6
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Ruminococcus_E"] <- -1

  p_ANCOMBC2 <- df_ANCOMBC2 %>% ggplot(aes(x = ANCOMBC2, y = PALM, label = label)) +
    geom_jitter(alpha = 0.7, aes(size = prevalence, color = group)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in ANCOMBC2 results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(-5, 160) + xlim(-5, 160) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_ANCOMBC2$hjust, vjust=df_ANCOMBC2$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = c(0.85,0.2)) + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_b_ANCOMBC2.pdf", width = 8.15, height = 6.65, bg = "white")

  p_ANCOMBC2

  dev.off()

  ## LinDA ----
  df_Linda <- data.frame(PALM = rowSums(PALM_coef_mat>0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(Linda = rowSums(Linda_coef_mat>0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_Linda$vjust[df_Linda$label == "Phocaeicola"] <- -1.5
  df_Linda$vjust[df_Linda$label == "Bacteroides"] <- -1.5
  df_Linda$vjust[df_Linda$label == "Bifidobacterium"] <- -5
  df_Linda$vjust[df_Linda$label == "Collinsella"] <- -4
  df_Linda$vjust[df_Linda$label == "Prevotella"] <- -5
  df_Linda$hjust[df_Linda$label == "Prevotella"] <- 1
  df_Linda$hjust[df_Linda$label == "Faecalibacterium"] <- 0
  df_Linda$vjust[df_Linda$label == "Faecalibacterium"] <- -1
  df_Linda$vjust[df_Linda$label == "Mediterraneibacter"] <- -1
  df_Linda$hjust[df_Linda$label == "Mediterraneibacter"] <- 1.3
  df_Linda$vjust[df_Linda$label == "Lachnospira"] <- -4
  df_Linda$hjust[df_Linda$label == "Lachnospira"] <- 2
  df_Linda$hjust[df_Linda$label == "Fusicatenibacter"] <- 1.1
  df_Linda$vjust[df_Linda$label == "Fusicatenibacter"] <- 1
  df_Linda$vjust[df_Linda$label == "Gemmiger"] <- 5
  df_Linda$vjust[df_Linda$label == "Ruminococcus_E"] <- 5
  df_Linda$hjust[df_Linda$label == "Ruminococcus_E"] <- 0
  df_Linda$vjust[df_Linda$label == "Parabacteroides"] <- 6
  df_Linda$hjust[df_Linda$label == "Parabacteroides"] <- 1

  p_Linda <- df_Linda %>% ggplot(aes(x = Linda, y = PALM, color = group, label = label)) +
    geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in Linda results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(-5, 160) + xlim(-5, 160) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_Linda$hjust, vjust=df_Linda$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_b_Linda.pdf", width = 8.15, height = 6.65, bg = "white")

  p_Linda

  dev.off()

  ## LM-CLR ----
  df_LMCLR <- data.frame(PALM = rowSums(PALM_coef_mat>0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(LMCLR = rowSums(LMCLR_coef_mat>0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_LMCLR$vjust[df_LMCLR$label == "Phocaeicola"] <- -1.5
  df_LMCLR$vjust[df_LMCLR$label == "Bacteroides"] <- -1.5
  df_LMCLR$vjust[df_LMCLR$label == "Bifidobacterium"] <- 10
  df_LMCLR$vjust[df_LMCLR$label == "Collinsella"] <- 6
  df_LMCLR$hjust[df_LMCLR$label == "Collinsella"] <- -1
  df_LMCLR$vjust[df_LMCLR$label == "Ruminococcus_E"] <- 6
  df_LMCLR$hjust[df_LMCLR$label == "Parabacteroides"] <- 1
  df_LMCLR$hjust[df_LMCLR$label == "Alistipes"] <- -1
  df_LMCLR$hjust[df_LMCLR$label == "Faecalibacterium"] <- 0
  df_LMCLR$vjust[df_LMCLR$label == "Mediterraneibacter"] <- -1
  df_LMCLR$hjust[df_LMCLR$label == "Mediterraneibacter"] <- 1.5
  df_LMCLR$vjust[df_LMCLR$label == "Lachnospira"] <- -2
  df_LMCLR$hjust[df_LMCLR$label == "Lachnospira"] <- 2
  df_LMCLR$hjust[df_LMCLR$label == "Fusicatenibacter"] <- 1.1
  df_LMCLR$vjust[df_LMCLR$label == "Fusicatenibacter"] <- 1



  p_LMCLR <- df_LMCLR %>% ggplot(aes(x = LMCLR, y = PALM, color = group, label = label)) +
    geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in LMCLR results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(-5, 160) + xlim(-5, 160) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_LMCLR$hjust, vjust=df_LMCLR$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_b_LMCLR.pdf", width = 8.15, height = 6.65, bg = "white")

  p_LMCLR

  dev.off()

  ## MaAsLin2 ----
  df_Maaslin2 <- data.frame(PALM = rowSums(PALM_coef_mat>0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(Maaslin2 = rowSums(Maaslin2_coef_mat>0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_Maaslin2$vjust[df_Maaslin2$label == "Phocaeicola"] <- -1.5
  df_Maaslin2$vjust[df_Maaslin2$label == "Bacteroides"] <- -1.5
  df_Maaslin2$hjust[df_Maaslin2$label == "Bacteroides"] <- 0.2
  df_Maaslin2$hjust[df_Maaslin2$label == "Agathobacter"] <- 1.2
  df_Maaslin2$vjust[df_Maaslin2$label == "Agathobacter"] <- 0
  df_Maaslin2$vjust[df_Maaslin2$label == "Roseburia"] <- -2
  df_Maaslin2$hjust[df_Maaslin2$label == "Bifidobacterium"] <- -0.2
  df_Maaslin2$vjust[df_Maaslin2$label == "Bifidobacterium"] <- 1
  df_Maaslin2$hjust[df_Maaslin2$label == "Faecalibacterium"] <- -0.2
  df_Maaslin2$vjust[df_Maaslin2$label == "Faecalibacterium"] <- 0
  df_Maaslin2$vjust[df_Maaslin2$label == "Prevotella"] <- -5
  df_Maaslin2$vjust[df_Maaslin2$label == "Escherichia"] <- -5
  df_Maaslin2$hjust[df_Maaslin2$label == "Escherichia"] <- 1.5
  df_Maaslin2$vjust[df_Maaslin2$label == "Fusicatenibacter"] <- -5
  df_Maaslin2$hjust[df_Maaslin2$label == "Fusicatenibacter"] <- 0.8
  df_Maaslin2$vjust[df_Maaslin2$label == "Mediterraneibacter"] <- 2
  df_Maaslin2$vjust[df_Maaslin2$label == "Ruminococcus_E"] <- 4
  df_Maaslin2$hjust[df_Maaslin2$label == "Ruminococcus_E"] <- -1
  df_Maaslin2$vjust[df_Maaslin2$label == "Parabacteroides"] <- 8
  df_Maaslin2$vjust[df_Maaslin2$label == "Gemmiger"] <- 4

  p_Maaslin2 <- df_Maaslin2 %>% ggplot(aes(x = Maaslin2, y = PALM, color = group, label = label)) +
    geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in Maaslin2 results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(-5, 160) + xlim(0, 160) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_Maaslin2$hjust, vjust=df_Maaslin2$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_b_Maaslin2.pdf", width = 8.15, height = 6.65, bg = "white")

  p_Maaslin2

  dev.off()

  # Figure c ----
  rm(list = ls())

  load("./Metabolites/Output/MTBL_PALM.Rdata")
  load("./Metabolites/Output/MTBL_others.Rdata")
  load("./Metabolites/Processed_data/MTBL.RData")

  feature.ids <- PALM.meta.coef$feature
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

  sample.id <- as.vector(unlist(lapply(data.for.lm, function(d){rownames(d$rel.abd)})))
  taxa <- NULL
  for(d in names(data.for.lm)){
    taxa <- c(taxa, colnames(data.for.lm[[d]]$rel.abd))
  }
  taxa <- unique(taxa)
  Otu_big_mat <- matrix(0, nrow = length(sample.id), ncol = length(taxa),
                        dimnames = list(sample.id, taxa))
  for(d in names(data.for.lm)){
    Otu_big_mat[rownames(data.for.lm[[d]]$rel.abd),colnames(data.for.lm[[d]]$rel.abd)] <- as.matrix(data.for.lm[[d]]$rel.abd)
  }
  tax.prevalence <- colMeans(Otu_big_mat / rowSums(Otu_big_mat))
  names(tax.prevalence) <- gsub(".*;g__","g__", names(tax.prevalence))

  ## ANCOMBC2
  df_ANCOMBC2 <- data.frame(PALM = rowSums(PALM_coef_mat<0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(ANCOMBC2 = rowSums(ANCOMBC2_coef_mat<0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Phocaeicola"] <- -8
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Ruminococcus_E"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Gemmiger"] <- 8
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Gemmiger"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Escherichia"] <- 5
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Escherichia"] <- -1.3
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Fusicatenibacter"] <- 0.3
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Fusicatenibacter"] <- 4
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Prevotella"] <- 0
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Prevotella"] <- 3
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Agathobacter"] <- 0.6
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Collinsella"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Collinsella"] <- 6
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Bacteroides"] <- 0.7
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Bacteroides"] <- -10
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Enterocloster"] <- -9
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Anaerostipes"] <- -8
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Anaerostipes"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Mediterraneibacter"] <- -9
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Mediterraneibacter"] <- 0.5
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Blautia_A"] <- -2.5
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Blautia_A"] <- 0
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Roseburia"] <- 0
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Bifidobacterium"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Bifidobacterium"] <- 1.5
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Lachnospira"] <- -2
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Lachnospira"] <- 4
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Parabacteroides"] <- -1.5
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Parabacteroides"] <- 0.5
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Faecalibacterium"] <- -4
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Faecalibacterium"] <- -1
  df_ANCOMBC2$vjust[df_ANCOMBC2$label == "Roseburia"] <- -4
  df_ANCOMBC2$hjust[df_ANCOMBC2$label == "Roseburia"] <- 0.8

  p_ANCOMBC2 <- df_ANCOMBC2 %>% ggplot(aes(x = ANCOMBC2, y = PALM, label = label)) +
    geom_jitter(alpha = 0.7, aes(size = prevalence, color = group)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in ANCOMBC2 results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(-5, 200) + xlim(-5, 200) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_ANCOMBC2$hjust, vjust=df_ANCOMBC2$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = c(0.85,0.2)) + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_c_ANCOMBC2.pdf", width = 8.15, height = 6.65, bg = "white")

  p_ANCOMBC2

  dev.off()

  ## LinDA ----
  df_Linda <- data.frame(PALM = rowSums(PALM_coef_mat<0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(Linda = rowSums(Linda_coef_mat<0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_Linda$vjust[df_Linda$label == "Phocaeicola"] <- -8
  df_Linda$hjust[df_Linda$label == "Phocaeicola"] <- 1.5
  df_Linda$vjust[df_Linda$label == "Bacteroides"] <- -8
  df_Linda$vjust[df_Linda$label == "Anaerostipes"] <- -9
  df_Linda$vjust[df_Linda$label == "Enterocloster"] <- -8
  df_Linda$hjust[df_Linda$label == "Enterocloster"] <- -1
  df_Linda$hjust[df_Linda$label == "Ruminococcus_E"] <- -1
  df_Linda$vjust[df_Linda$label == "Collinsella"] <- 4
  df_Linda$vjust[df_Linda$label == "Fusicatenibacter"] <- 3
  df_Linda$hjust[df_Linda$label == "Fusicatenibacter"] <- 0.6
  df_Linda$vjust[df_Linda$label == "Prevotella"] <- 3
  df_Linda$hjust[df_Linda$label == "Prevotella"] <- 0.2
  df_Linda$hjust[df_Linda$label == "Bifidobacterium"] <- 1.2
  df_Linda$vjust[df_Linda$label == "Bifidobacterium"] <- -5
  df_Linda$hjust[df_Linda$label == "Roseburia"] <- 1.2
  df_Linda$hjust[df_Linda$label == "Faecalibacterium"] <- 0
  df_Linda$vjust[df_Linda$label == "Faecalibacterium"] <- 4
  df_Linda$vjust[df_Linda$label == "Gemmiger"] <- -4
  df_Linda$hjust[df_Linda$label == "Gemmiger"] <- -1
  df_Linda$vjust[df_Linda$label == "Escherichia"] <- -6
  df_Linda$hjust[df_Linda$label == "Escherichia"] <- -1
  df_Linda$vjust[df_Linda$label == "Blautia_A"] <- -4
  df_Linda$hjust[df_Linda$label == "Blautia_A"] <- 0.2
  df_Linda$vjust[df_Linda$label == "Mediterraneibacter"] <- -2
  df_Linda$hjust[df_Linda$label == "Mediterraneibacter"] <- 1.3
  df_Linda$vjust[df_Linda$label == "Lachnospira"] <- -4
  df_Linda$hjust[df_Linda$label == "Lachnospira"] <- 2
  df_Linda$vjust[df_Linda$label == "Parabacteroides"] <- -7
  df_Linda$hjust[df_Linda$label == "Parabacteroides"] <- 1

  p_Linda <- df_Linda %>% ggplot(aes(x = Linda, y = PALM, color = group, label = label)) +
    geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in Linda results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(-5, 200) + xlim(-5, 200) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_Linda$hjust, vjust=df_Linda$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_c_Linda.pdf", width = 8.15, height = 6.65, bg = "white")

  p_Linda

  dev.off()

  ## LM-CLR ----
  df_LMCLR <- data.frame(PALM = rowSums(PALM_coef_mat<0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(LMCLR = rowSums(LMCLR_coef_mat<0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_LMCLR$vjust[df_LMCLR$label == "Phocaeicola"] <- -8
  df_LMCLR$hjust[df_LMCLR$label == "Phocaeicola"] <- 1.5
  df_LMCLR$vjust[df_LMCLR$label == "Bacteroides"] <- -8
  df_LMCLR$vjust[df_LMCLR$label == "Anaerostipes"] <- -9
  df_LMCLR$vjust[df_LMCLR$label == "Enterocloster"] <- -8
  df_LMCLR$hjust[df_LMCLR$label == "Enterocloster"] <- -1
  df_LMCLR$hjust[df_LMCLR$label == "Ruminococcus_E"] <- -1
  df_LMCLR$vjust[df_LMCLR$label == "Collinsella"] <- 4
  df_LMCLR$hjust[df_LMCLR$label == "Collinsella"] <- 0
  df_LMCLR$vjust[df_LMCLR$label == "Fusicatenibacter"] <- 3
  df_LMCLR$hjust[df_LMCLR$label == "Fusicatenibacter"] <- 0.6
  df_LMCLR$vjust[df_LMCLR$label == "Prevotella"] <- 3
  df_LMCLR$hjust[df_LMCLR$label == "Prevotella"] <- 0.2
  df_LMCLR$hjust[df_LMCLR$label == "Bifidobacterium"] <- 1.2
  df_LMCLR$vjust[df_LMCLR$label == "Bifidobacterium"] <- -5
  df_LMCLR$hjust[df_LMCLR$label == "Roseburia"] <- 1.3
  df_LMCLR$vjust[df_LMCLR$label == "Roseburia"] <- -0.5
  df_LMCLR$hjust[df_LMCLR$label == "Faecalibacterium"] <- 0
  df_LMCLR$vjust[df_LMCLR$label == "Faecalibacterium"] <- 4
  df_LMCLR$vjust[df_LMCLR$label == "Gemmiger"] <- -4
  df_LMCLR$hjust[df_LMCLR$label == "Gemmiger"] <- -1
  df_LMCLR$vjust[df_LMCLR$label == "Escherichia"] <- -6
  df_LMCLR$hjust[df_LMCLR$label == "Escherichia"] <- -1
  df_LMCLR$vjust[df_LMCLR$label == "Blautia_A"] <- -4
  df_LMCLR$hjust[df_LMCLR$label == "Blautia_A"] <- 0.2
  df_LMCLR$vjust[df_LMCLR$label == "Mediterraneibacter"] <- -2
  df_LMCLR$hjust[df_LMCLR$label == "Mediterraneibacter"] <- 1.3
  df_LMCLR$vjust[df_LMCLR$label == "Lachnospira"] <- -4
  df_LMCLR$hjust[df_LMCLR$label == "Lachnospira"] <- 2
  df_LMCLR$vjust[df_LMCLR$label == "Parabacteroides"] <- -7
  df_LMCLR$hjust[df_LMCLR$label == "Parabacteroides"] <- 1

  p_LMCLR <- df_LMCLR %>% ggplot(aes(x = LMCLR, y = PALM, color = group, label = label)) +
    geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in LMCLR results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(-5, 200) + xlim(-5, 200) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_LMCLR$hjust, vjust=df_LMCLR$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_c_LMCLR.pdf", width = 8.15, height = 6.65, bg = "white")

  p_LMCLR

  dev.off()

  ## MaAsLin2 ----
  df_Maaslin2 <- data.frame(PALM = rowSums(PALM_coef_mat<0)) %>%
    tibble::rownames_to_column(var = "Genera") %>% dplyr::left_join(
      data.frame(Maaslin2 = rowSums(Maaslin2_coef_mat<0)) %>%
        tibble::rownames_to_column(var = "Genera"), by = "Genera") %>%
    dplyr::mutate(Genera = gsub(".*;g__","g__",Genera)) %>%
    dplyr::arrange(desc(PALM)) %>%
    dplyr::mutate(vjust = rep(-1, 101)) %>% dplyr::mutate(hjust = rep(0.5, 101)) %>%
    dplyr::left_join(data.frame(tax.prevalence = tax.prevalence) %>% tibble::rownames_to_column("Genera"), by = "Genera") %>%
    dplyr::mutate(prevalence = dplyr::case_when(tax.prevalence >=  0.1 ~ "1 \u00D7 10 <sup>-1</sup>",
                                                tax.prevalence < 0.1 & tax.prevalence >= 0.01 ~ "1 \u00D7 10 <sup>-2</sup>",
                                                tax.prevalence < 0.01 & tax.prevalence >= 0.001 ~ "1 \u00D7 10 <sup>-3</sup>",
                                                tax.prevalence < 0.001 & tax.prevalence >= 0.0001 ~ "1 \u00D7 10 <sup>-4</sup>",
                                                TRUE ~ "")) %>%
    dplyr::mutate(group = dplyr::case_when(tax.prevalence >=  0.01 ~ "abundant genera",
                                           TRUE ~ "rare genera")) %>%
    dplyr::mutate(label = dplyr::case_when(group == "abundant genera" ~ gsub("g__","", Genera),
                                           TRUE ~ ""))

  df_Maaslin2$vjust[df_Maaslin2$label == "Phocaeicola"] <- -8
  df_Maaslin2$hjust[df_Maaslin2$label == "Phocaeicola"] <- 1.5
  df_Maaslin2$vjust[df_Maaslin2$label == "Bacteroides"] <- -8
  df_Maaslin2$hjust[df_Maaslin2$label == "Bacteroides"] <- 0
  df_Maaslin2$vjust[df_Maaslin2$label == "Anaerostipes"] <- -9
  df_Maaslin2$vjust[df_Maaslin2$label == "Enterocloster"] <- -8
  df_Maaslin2$hjust[df_Maaslin2$label == "Enterocloster"] <- -1
  df_Maaslin2$hjust[df_Maaslin2$label == "Ruminococcus_E"] <- -1
  df_Maaslin2$vjust[df_Maaslin2$label == "Collinsella"] <- 4
  df_Maaslin2$hjust[df_Maaslin2$label == "Collinsella"] <- -0.1
  df_Maaslin2$vjust[df_Maaslin2$label == "Fusicatenibacter"] <- 3
  df_Maaslin2$hjust[df_Maaslin2$label == "Fusicatenibacter"] <- 0.6
  df_Maaslin2$vjust[df_Maaslin2$label == "Prevotella"] <- 3
  df_Maaslin2$hjust[df_Maaslin2$label == "Prevotella"] <- 0.2
  df_Maaslin2$hjust[df_Maaslin2$label == "Bifidobacterium"] <- 1.2
  df_Maaslin2$vjust[df_Maaslin2$label == "Bifidobacterium"] <- -5
  df_Maaslin2$hjust[df_Maaslin2$label == "Roseburia"] <- 1
  df_Maaslin2$hjust[df_Maaslin2$label == "Faecalibacterium"] <- 0
  df_Maaslin2$vjust[df_Maaslin2$label == "Faecalibacterium"] <- 4
  df_Maaslin2$vjust[df_Maaslin2$label == "Gemmiger"] <- -4
  df_Maaslin2$hjust[df_Maaslin2$label == "Gemmiger"] <- -1
  df_Maaslin2$vjust[df_Maaslin2$label == "Escherichia"] <- -6
  df_Maaslin2$hjust[df_Maaslin2$label == "Escherichia"] <- -1.5
  df_Maaslin2$vjust[df_Maaslin2$label == "Blautia_A"] <- -4
  df_Maaslin2$hjust[df_Maaslin2$label == "Blautia_A"] <- -1
  df_Maaslin2$vjust[df_Maaslin2$label == "Mediterraneibacter"] <- -2.5
  df_Maaslin2$hjust[df_Maaslin2$label == "Mediterraneibacter"] <- 1.4
  df_Maaslin2$vjust[df_Maaslin2$label == "Lachnospira"] <- -4
  df_Maaslin2$hjust[df_Maaslin2$label == "Lachnospira"] <- 2
  df_Maaslin2$vjust[df_Maaslin2$label == "Parabacteroides"] <- -7
  df_Maaslin2$hjust[df_Maaslin2$label == "Parabacteroides"] <- 1

  p_Maaslin2 <- df_Maaslin2 %>% ggplot(aes(x = Maaslin2, y = PALM, color = group, label = label)) +
    geom_jitter(alpha = 0.7, aes(size = prevalence)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    xlab("Count in Maaslin2 results")+ ylab("Count in PALM results") +
    theme_minimal() + ylim(-5, 200) + xlim(0, 200) + ggtitle("") +
    scale_color_manual(breaks = c("rare genera", "abundant genera"),
                       values=c("gray", "red")) +
    scale_size_manual(values = c(8,4,2,1),
                      breaks = c("1 \u00D7 10 <sup>-1</sup>",
                                 "1 \u00D7 10 <sup>-2</sup>",
                                 "1 \u00D7 10 <sup>-3</sup>",
                                 "1 \u00D7 10 <sup>-4</sup>")) +
    geom_text(hjust=df_Maaslin2$hjust, vjust=df_Maaslin2$vjust, size = 4, color = "black") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill=NULL, colour='black'),
          plot.margin = unit(c(0.2, 0.3, 0.2, 0.2), 'cm'),
          axis.text.y = element_text(size = 15),
          axis.text.x = element_markdown(size = 15),
          legend.title = element_markdown(hjust = 0.5, size = 15),
          legend.text = element_markdown(size = 15),
          legend.position = "none") + labs(size="Trimmed average <br> proportion") +
    guides(size = guide_legend(override.aes = list(color = c("red", "red", "gray","gray"))),
           color = "none")

  pdf("./Figure/FigureS6_c_Maaslin2.pdf", width = 8.15, height = 6.65, bg = "white")

  p_Maaslin2

  dev.off()
