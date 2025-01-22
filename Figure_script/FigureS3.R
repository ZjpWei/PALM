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
  pdf("./Figure/FigureS3_seed2022.pdf", width = 10.35, height = 6.95, bg = "white")

  ggpubr::ggarrange(g.plots[[1]], g.plots[[2]], g.plots[[3]],
                    g.plots[[4]], g.plots[[5]], nrow = 2, ncol = 3,
                    common.legend = TRUE, legend = "bottom")

  dev.off()

