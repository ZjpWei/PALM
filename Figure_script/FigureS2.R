# =============================================== #
#                   Figure S2                     #
# =============================================== #

  # Packages ----
  library("phyloseq")
  library("ggplot2")
  library("vegan")

  # Figure a ----
  rm(list = ls())

  load("./Data/CRC_data/data/meta.Rdata")
  load("./CRC/Processed_data/data.org.K849.Rdata")
  studymeta <- unique(meta$Study)
  feature.table <- NULL
  outcome <- NULL
  Study <- NULL
  for(l in 1:length(data.rel)){
    feature.table <- rbind(feature.table, data.rel[[l]]$Y/rowSums(data.rel[[l]]$Y) )
    outcome <- c(outcome, data.rel[[l]]$X)
    Study <- c(Study, rep(paste0("CRC",as.character(l)), length(data.rel[[l]]$X)))
  }
  rownames(feature.table) <- paste0("sa", as.character(1:length(outcome)))
  names(outcome) <- rownames(feature.table)
  meta.data = data.frame(labels = outcome)
  rownames(feature.table) <- rownames(meta.data)

  # Perform PERMANOVA (optional)
  # distance_matrix <- vegdist(feature.table, method = "bray")
  #
  # adonis_result <- adonis2(formula = distance_matrix ~ Study,
  #                          data = data.frame(Group = Study),
  #                          permutations = 100000)

  ## PCoA plot
  feature.table = t(feature.table)
  OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
  META = sample_data(meta.data)
  PHYSEQ = phyloseq(OTU, META)

  # Calculate ordination
  iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")

  PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"], PCo2 = iMDS$vectors[,"Axis.2"], Study = Study,
                     outcome = as.character(outcome))

  ### Generate figures ----
  pdf("./Figure/FigureS2_a.pdf", width = 6.13, height = 4.39, bg = "white")

  ggplot(PCoA, aes(x=PCo1, y=PCo2, color=Study)) + geom_point() +
    scale_color_manual(values = c("#E76BF3","#00B0F6","#00BF7D","#A3A500","#F8766D"),
                       breaks = c("CRC1", "CRC2", "CRC3", "CRC4", "CRC5")) +
    theme(text = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.key = element_rect(fill = "white"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.position = "right")

  dev.off()

  # Figure b ----
  rm(list = ls())

  load("./MTBL.RData")
  feature.table <- NULL
  Study <- NULL
  outcome <- NULL
  for(k in 1:length(datasets)){
    d <- datasets[k]
    dims <- dim(data.for.lm[[d]]$rel.abd)
    feature.table.single <- matrix(0, nrow = dims[1], ncol = length(Genera.id),
                                   dimnames = list(rownames(data.for.lm[[d]]$rel.abd), Genera.id))

    feature.table.single[,colnames(data.for.lm[[d]]$rel.abd)] <- as.matrix(data.for.lm[[d]]$rel.abd/rowSums(data.for.lm[[d]]$rel.abd))
    ## Real data
    feature.table <- rbind(feature.table, feature.table.single)
    Study <- c(Study, rep(paste0("MTBL", k), nrow(data.for.lm[[d]]$rel.abd)))
  }
  meta.data = data.frame(Study = Study)
  rownames(meta.data) <- rownames(feature.table)

  # Perform PERMANOVA (optional)
  # distance_matrix <- vegdist(feature.table, method = "bray")
  #
  # adonis_result <- adonis2(formula = distance_matrix ~ Study,
  #                          data = data.frame(Study = Study),
  #                          permutations = 100000)

  ## PCoA plot
  feature.table = t(feature.table)
  OTU = otu_table(as.matrix(feature.table), taxa_are_rows = TRUE)
  META = sample_data(meta.data)
  PHYSEQ = phyloseq(OTU, META)

  # Calculate ordination
  iMDS <- phyloseq::ordinate(PHYSEQ, "PCoA", distance = "bray")

  PCoA <- data.frame(PCo1 = iMDS$vectors[,"Axis.1"], PCo2 = iMDS$vectors[,"Axis.2"], Study = Study)

  ### Generate Figure
  pdf("./Figure/FigureS2_b.pdf", width = 6.13, height = 4.39, bg = "white")

  ggplot(PCoA, aes(x=PCo1, y=PCo2, color=Study)) + geom_point() +
    scale_color_manual(values = c("#FF61CC","#C77CFF","#00A9FF","#00BFC4",
                                  "#00BE67","#7CAE00","#CD9600","#F8766D"),
                       breaks = c("MTBL1", "MTBL2", "MTBL3", "MTBL4",
                                  "MTBL5", "MTBL6", "MTBL7", "MTBL8")) +
    theme(text = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill = 'white'),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.key = element_rect(fill = "white"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.position = "right")

  dev.off()
