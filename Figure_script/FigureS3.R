# =============================================== #
#                    Figure S3                    #
# =============================================== #

  # Packages ----
  library("phyloseq")
  library("ggplot2")
  library("MIDASim")
  library("vegan")

  # General ----
  rm(list = ls())

  ## Original ----
  set.seed(2022)
  load("./Data/CRC_data/data/meta.Rdata")
  load("./CRC/Processed_data/data.org.K401.Rdata")
  studymeta <- unique(meta$Study)

  ## Bray dostance
  g.plots <- list()
  for(l in 1:length(data.rel)){
    non.id <- which(colSums(data.rel[[l]]$Y) !=0)

    ## Real data
    feature.table <- data.rel[[l]]$Y[,non.id]/rowSums(data.rel[[l]]$Y[,non.id])
    outcome <- data.rel[[l]]$X
    Group <- rep("Original data", length(data.rel[[l]]$X))

    ## Sim data
    count.ibd.setup <- MIDASim.setup(data.rel[[l]]$Y[,non.id], mode = 'nonparametric', n.break.ties = 100)
    count.ibd.modified <- MIDASim.modify(count.ibd.setup)
    simulated.data <- MIDASim(count.ibd.modified)

    feature.table <- rbind(feature.table,  simulated.data$sim_count/rowSums( simulated.data$sim_count))
    outcome <- c(outcome, data.rel[[l]]$X)
    Group <- c(Group, rep("Simulated data", length(data.rel[[l]]$X)))

    ## add labels
    rownames(feature.table) <- paste0("sa", as.character(1:length(outcome)))
    names(outcome) <- rownames(feature.table)
    meta.data = data.frame(labels = outcome)
    rownames(feature.table) <- rownames(meta.data)

    # Perform PERMANOVA (optional)
    # distance_matrix <- vegdist(feature.table, method = "bray")
    #
    # adonis_result <- adonis2(formula = distance_matrix ~ Group,
    #                          data = data.frame(Group = Group),
    #                          permutations = 999)

    ## PCoA plot
    feature.table = t(feature.table)
    OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
    META = sample_data(meta.data)
    PHYSEQ = phyloseq(OTU, META)

    # Calculate ordination
    iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")

    PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"],
                       PCo2 = iMDS$vectors[,"Axis.2"], Data = Group,
                       outcome = as.character(outcome))

    g.plots[[l]] <-  ggplot(PCoA, aes(x=PCo1, y=PCo2, color=Data)) + geom_point() +
      ggtitle(paste0("CRC ",l)) +
      theme(text = element_text(size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_rect(fill = 'white'),
            panel.border = element_rect(colour = "black", fill=NA),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.key = element_rect(fill = "white"),
            legend.title = element_blank(),
            legend.text = element_text(size = 13),
            legend.position = "right")
  }

  ## Generate figures ----
  pdf("./Figure/FigureS3_a.pdf", width = 17.25, height = 4, bg = "white")

  ggpubr::ggarrange(g.plots[[1]], g.plots[[2]], g.plots[[3]],
                    g.plots[[4]], g.plots[[5]], nrow = 1, ncol = 5,
                    common.legend = TRUE, legend = "bottom")

  dev.off()


  ## Jaccard distance
  set.seed(2022)
  g.plots <- list()
  for(l in 1:length(data.rel)){
    non.id <- which(colSums(data.rel[[l]]$Y) !=0)

    ## Real data
    feature.table <- data.rel[[l]]$Y[,non.id]/rowSums(data.rel[[l]]$Y[,non.id])
    outcome <- data.rel[[l]]$X
    Group <- rep("Original data", length(data.rel[[l]]$X))

    ## Sim data
    count.ibd.setup <- MIDASim.setup(data.rel[[l]]$Y[,non.id], mode = 'nonparametric', n.break.ties = 100)
    count.ibd.modified <- MIDASim.modify(count.ibd.setup)
    simulated.data <- MIDASim(count.ibd.modified)

    feature.table <- rbind(feature.table,  simulated.data$sim_count/rowSums( simulated.data$sim_count))
    outcome <- c(outcome, data.rel[[l]]$X)
    Group <- c(Group, rep("Simulated data", length(data.rel[[l]]$X)))

    ## add labels
    rownames(feature.table) <- paste0("sa", as.character(1:length(outcome)))
    names(outcome) <- rownames(feature.table)
    meta.data = data.frame(labels = outcome)
    rownames(feature.table) <- rownames(meta.data)

    # Perform PERMANOVA (optional)
    distance_matrix <- vegdist(feature.table, method = "jaccard")

    adonis_result <- adonis2(formula = distance_matrix ~ Group,
                             data = data.frame(Group = Group),
                             permutations = 999)

    print(adonis_result)

    ## PCoA plot
    feature.table = t(feature.table)
    OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
    META = sample_data(meta.data)
    PHYSEQ = phyloseq(OTU, META)

    # Calculate ordination
    iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "jaccard")

    PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"],
                       PCo2 = iMDS$vectors[,"Axis.2"], Data = Group,
                       outcome = as.character(outcome))

    g.plots[[l]] <-  ggplot(PCoA, aes(x=PCo1, y=PCo2, color=Data)) + geom_point() +
      ggtitle(paste0("CRC ",l)) +
      theme(text = element_text(size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank(),
            panel.background = element_rect(fill = 'white'),
            panel.border = element_rect(colour = "black", fill=NA),
            plot.title = element_text(hjust = 0.5, size = 15),
            legend.key = element_rect(fill = "white"),
            legend.title = element_blank(),
            legend.text = element_text(size = 13),
            legend.position = "right")
  }

  ## Generate figures ----
  pdf("./Figure/FigureS3_b.pdf", width = 17.25, height = 4, bg = "white")

  ggpubr::ggarrange(g.plots[[1]], g.plots[[2]], g.plots[[3]],
                    g.plots[[4]], g.plots[[5]], nrow = 1, ncol = 5,
                    common.legend = TRUE, legend = "bottom")

  dev.off()


  ## Violin plot 1
  set.seed(2022)
  g.plots <- list()
  for(l in 1:length(data.rel)){
    non.id <- which(colSums(data.rel[[l]]$Y) !=0)

    ## Real data
    feature.table <- data.rel[[l]]$Y[,non.id]
    seq_real <- rowSums(feature.table)

    ## Sim data
    count.ibd.setup <- MIDASim.setup(data.rel[[l]]$Y[,non.id], mode = 'nonparametric', n.break.ties = 100)
    count.ibd.modified <- MIDASim.modify(count.ibd.setup)
    simulated.data <- MIDASim(count.ibd.modified)
    seq_sim <- rowSums(simulated.data$sim_count)

    # Combine into a tidy data frame
    df.compare <- data.frame(
      Depth = c(seq_real, seq_sim),
      Source = rep(c("Original data", "Simulated data"), each = length(seq_real))
    )

    g.plots[[l]] <- ggplot(df.compare, aes(x = Source, y = Depth, fill = Source)) +
      geom_violin(trim = FALSE, alpha = 0.7) +
      geom_boxplot(width = 0.1, color = "black", alpha = 0.6, outlier.shape = NA) +  # optional overlay
      theme_bw() +
      theme(
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
      ) +
      scale_y_continuous(
        labels = scales::scientific,      # ← 显示为 1e5 格式
        expand = expansion(mult = c(0, 0.05))
      ) +
      labs(
        x = "",
        y = "",
        title = paste0("CRC ", l)
      )
  }

  ## Generate figures ----
  pdf("./Figure/FigureS3_c.pdf", width = 20, height = 4.5, bg = "white")

  ggpubr::ggarrange(g.plots[[1]], g.plots[[2]], g.plots[[3]],
                    g.plots[[4]], g.plots[[5]], nrow = 1, ncol = 5,
                    common.legend = TRUE, legend = "bottom")

  dev.off()


  ## zero-proportion violin plot
  set.seed(2022)
  g.plots <- list()
  for(l in 1:length(data.rel)){
    non.id <- which(colSums(data.rel[[l]]$Y) !=0)

    ## Real data
    feature.table <- data.rel[[l]]$Y[,non.id]
    seq_real <- colMeans(feature.table!=0)

    ## Sim data
    count.ibd.setup <- MIDASim.setup(data.rel[[l]]$Y[,non.id], mode = 'nonparametric', n.break.ties = 100)
    count.ibd.modified <- MIDASim.modify(count.ibd.setup)
    simulated.data <- MIDASim(count.ibd.modified)
    seq_sim <- colMeans(simulated.data$sim_count!=0)

    # Combine into a tidy data frame
    df.compare <- data.frame(
      Depth = c(seq_real, seq_sim),
      Source = rep(c("Original data", "Simulated data"), each = length(seq_real))
    )

    g.plots[[l]] <- ggplot(df.compare, aes(x = Source, y = Depth, fill = Source)) +
      geom_violin(trim = TRUE, alpha = 0.7) +
      geom_boxplot(width = 0.1, color = "black", alpha = 0.6, outlier.shape = NA) +  # optional overlay
      coord_cartesian(ylim = c(0, 1)) +          # clamp the visible range
      theme_bw() +
      theme(
        axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
      ) +
      labs(
        x = "",
        y = "",
        title = paste0("CRC ", l)
      )
  }

  ## Generate figures ----
  pdf("./Figure/FigureS3_d.pdf", width = 20, height = 4.5, bg = "white")

  ggpubr::ggarrange(g.plots[[1]], g.plots[[2]], g.plots[[3]],
                    g.plots[[4]], g.plots[[5]], nrow = 1, ncol = 5,
                    common.legend = TRUE, legend = "bottom")

  dev.off()
