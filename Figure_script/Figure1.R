# =============================================== #
#                   Figure 1                      #
# =============================================== #

  # Packages ----
  library('ggplot2')
  library("tidyverse")
  library("latex2exp")
  library("ggtext")

  # Main ----
  rm(list = ls())

  # Setup path.
  data.loc <- "./Simulation/Independent/"

  PRC_all <- NULL
  for(pos.lst in c(0.5, 1)){
    for(u.lst in c(0, 1)){
      for(tag in c(0.05, 0.1, 0.15, 0.2)){
        PRC_tmp <- NULL
        for(s in 1:100){
          if(file.exists(paste0(data.loc, "Sim_Ka", tag, "_Pos", pos.lst, "_mu", u.lst, "_", as.character(s), ".Rdata"))){
            load(paste0(data.loc, "Sim_Ka", tag, "_Pos", pos.lst, "_mu", u.lst, "_", as.character(s), ".Rdata"))
            PRC$ep.fdr[is.na(PRC$ep.fdr)] <- 0
            PRC_tmp <- rbind(PRC_tmp, PRC)
          }
        }
        tmp.prc <- data.frame(PRC_tmp) %>% group_by(method, Settings, tax.type, Study) %>%
          summarize(FDR = mean(ep.fdr), sd.FDR = sd(ep.fdr),
                    Power = mean(ep.power), sd.Power = sd(ep.power),
                    het = mean(ep.het.fdr), sd.het = sd(ep.het.fdr))
        tmp.prc$x.label <- tag
        if(pos.lst == 0.5){
          tmp.prc$pos <- paste0("Balanced +/-")
        }else if(pos.lst == 1){
          tmp.prc$pos <- paste0("Dominant +")
        }
        if(u.lst == 0){
          tmp.prc$u <- paste0("Even depth")
        }else if(u.lst == 1){
          tmp.prc$u <- paste0("Uneven depth")
        }
        PRC_all <- rbind(PRC_all, tmp.prc)
      }
    }
  }
  PRC_all$x.label <- factor(PRC_all$x.label, levels = unique(PRC_all$x.label), ordered = TRUE)
  PRC_all$Method <- factor(PRC_all$method, levels = c("ANCOM-BC2", "LinDA", "LM-CLR", "MaAsLin2", "PALM"), ordered = TRUE)

  ## Generate figures large/species
  p.large.fdr.species <- PRC_all %>% dplyr::filter(Settings == "large", tax.type == "species", Study == "meta") %>%
    ggplot(aes(x=x.label, y=FDR, group = Method)) +
    geom_hline(aes(yintercept = 0.05),colour="#990000", linetype="dashed") +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(FDR - sd.FDR, 0), ymax = pmin(FDR + sd.FDR, 1), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylim(0, 1) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  p.large.pw.species <- PRC_all %>% dplyr::filter(Settings == "large", tax.type == "species", Study == "meta") %>%
    ggplot(aes(x=x.label, y=Power, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(Power - sd.Power, 0.5), ymax = pmin(Power + sd.Power, 1), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylim(0.5, 1) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  p.large.het.species <- PRC_all %>% dplyr::filter(Settings == "large", tax.type == "species", Study == "meta") %>%
    ggplot(aes(x=x.label, y=het, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(het - sd.het, 0), ymax = pmin(het + sd.het, 0.4), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylab("Proportion of het. features") + ylim(0, 0.4) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  ## Generate figures small/species
  p.small.fdr.species <- PRC_all %>% dplyr::filter(Settings == "small", tax.type == "species", Study == "meta") %>%
    ggplot(aes(x=x.label, y=FDR, group = Method)) +
    geom_hline(aes(yintercept = 0.05),colour="#990000", linetype="dashed") +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(FDR - sd.FDR, 0), ymax = pmin(FDR + sd.FDR, 1), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylim(0, 1) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  p.small.pw.species <- PRC_all %>% dplyr::filter(Settings == "small", tax.type == "species", Study == "meta") %>%
    ggplot(aes(x=x.label, y=Power, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(Power - sd.Power, 0.5), ymax = pmin(Power + sd.Power, 1), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylim(0.5, 1) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  p.small.het.species <- PRC_all %>% dplyr::filter(Settings == "small", tax.type == "species", Study == "meta") %>%
    ggplot(aes(x=x.label, y=het, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(het - sd.het, 0), ymax = pmin(het + sd.het, 0.4), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylab("Proportion of het. features") + ylim(0, 0.4) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  ## Generate figures large/genus
  p.large.fdr.genus <- PRC_all %>% dplyr::filter(Settings == "large", tax.type == "genus", Study == "meta") %>%
    ggplot(aes(x=x.label, y=FDR, group = Method)) +
    geom_hline(aes(yintercept = 0.05),colour="#990000", linetype="dashed") +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(FDR - sd.FDR, 0), ymax = pmin(FDR + sd.FDR, 1), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylim(0, 1) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  p.large.pw.genus <- PRC_all %>% dplyr::filter(Settings == "large", tax.type == "genus", Study == "meta") %>%
    ggplot(aes(x=x.label, y=Power, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(Power - sd.Power, 0.5), ymax = pmin(Power + sd.Power, 1), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylim(0.5, 1) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  p.large.het.genus <- PRC_all %>% dplyr::filter(Settings == "large", tax.type == "genus", Study == "meta") %>%
    ggplot(aes(x=x.label, y=het, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(het - sd.het, 0), ymax = pmin(het + sd.het, 0.4), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylab("Proportion of het. features") + ylim(0, 0.4) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  ## Generate figures small/genus
  p.small.fdr.genus <- PRC_all %>% dplyr::filter(Settings == "small", tax.type == "genus", Study == "meta") %>%
    ggplot(aes(x=x.label, y=FDR, group = Method)) +
    geom_hline(aes(yintercept = 0.05),colour="#990000", linetype="dashed") +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(FDR - sd.FDR, 0), ymax = pmin(FDR + sd.FDR, 1), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylim(0, 1) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  p.small.pw.genus <- PRC_all %>% dplyr::filter(Settings == "small", tax.type == "genus", Study == "meta") %>%
    ggplot(aes(x=x.label, y=Power, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(Power - sd.Power, 0.5), ymax = pmin(Power + sd.Power, 1), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) + ylim(0.5, 1) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="none",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  p.small.het.genus <- PRC_all %>% dplyr::filter(Settings == "small", tax.type == "genus", Study == "meta") %>%
    ggplot(aes(x=x.label, y=het, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method),
              position = position_dodge(width = 0.3)) +
    geom_point(size = 2,aes(color = Method), pch = 18,
               position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = pmax(het - sd.het, 0), ymax = pmin(het + sd.het, 0.4), color = Method), width = 0.3,
                  position = position_dodge(width = 0.3)) +  ylab("Proportion of het. features") + ylim(0, 0.4) +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position ="bottom",
          legend.box="vertical",
          axis.text = element_text(size = 18),
                    legend.title = element_blank(),
          legend.text = element_text(size = 18),
          strip.text = element_markdown(size = 18),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") + facet_grid(u ~ pos) +
    guides(color = guide_legend(order = 1, nrow = 1)) +
    labs(linetype="Data")

  # Generate figures ----
  pp1_1 <- ggpubr::annotate_figure(
    ggpubr::ggarrange(p.large.fdr.species, p.large.fdr.genus, p.small.fdr.species, p.small.fdr.genus,
                      nrow = 1, ncol = 4,  common.legend = TRUE, legend = "none", label.y = "AUPRC"),
    bottom = ggpubr::text_grob("Proportion of differential AA features", size = 20, face = "bold"),
    left = ggpubr::text_grob("FDR", size = 20, rot = 90, hjust = 0.2, face = "bold"),
  )

  pp1_2 <- ggpubr::annotate_figure(
    ggpubr::ggarrange(p.large.pw.species, p.large.pw.genus, p.small.pw.species, p.small.pw.genus,
                      nrow = 1, ncol = 4,  common.legend = TRUE, legend = "none", label.y = "AUPRC"),
    bottom = ggpubr::text_grob("Proportion of differential AA features", size = 20, face = "bold"),
    left = ggpubr::text_grob("Power", size = 20, rot = 90, hjust = 0.2, face = "bold"),
  )

  pp1_3 <- ggpubr::annotate_figure(
    ggpubr::ggarrange(p.large.het.species, p.large.het.genus, p.small.het.species, p.small.het.genus,
                      nrow = 1, ncol = 4,  common.legend = TRUE, legend = "none", label.y = "AUPRC") ,
    bottom = ggpubr::text_grob("Proportion of differential AA features", size = 20, face = "bold"),
    left = ggpubr::text_grob("Proportion of het. features", size = 20, rot = 90, hjust = 0.4, face = "bold"),
  )

  pdf("./Figure/Figure1_FDR.pdf", width = 21, height = 5.6, bg = "white")

  pp1_1

  dev.off()

  pdf("./Figure/Figure1_Power.pdf", width = 21, height = 5.6, bg = "white")

  pp1_2

  dev.off()

  pdf("./Figure/Figure1_het.pdf", width = 21, height = 5.6, bg = "white")

  pp1_3

  dev.off()

