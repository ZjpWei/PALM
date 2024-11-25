# =============================================== #
#             Method analysis time                #
# =============================================== #

  ## K = 401, large
  run_time_ll <- data.frame(time = c(c(1, 10, 100, 1000) * 19639.2,
                                     c(1, 10, 100, 1000) * 480.73,
                                     c(1, 10, 100, 1000) * 60.641,
                                     c(1, 10, 100, 1000) * 49.843,
                                     c(1, 10, 100, 1000) * 75.180) / 60 / 60,
  num_cov = rep(c(2, 3, 4, 5), 5),
  Method = c(rep("ANCOM-BC2", 4),
             rep("MaAsLin2", 4),
             rep("LM-CLR", 4),
             rep("LinDA", 4),
             rep("PALM", 4)),
  dim = c(rep("species K = 401/large", 20)))

  ## K = 401, small
  run_time_ls <- data.frame(time = c(c(1, 10, 100, 1000) * 20012.5,
                                     c(1, 10, 100, 1000) * 441.26,
                                     c(1, 10, 100, 1000) * 55.453,
                                     c(1, 10, 100, 1000) * 44.959,
                                     c(1, 10, 100, 1000) * 11.521) / 60 / 60,

  num_cov = rep(c(2, 3, 4, 5), 5),
  Method = c(rep("ANCOM-BC2", 4),
             rep("MaAsLin2", 4),
             rep("LM-CLR", 4),
             rep("LinDA", 4),
             rep("PALM", 4)),
  dim = c(rep("species K = 401/small", 20)))

  ## K = 92, large
  run_time_sl <- data.frame(time = c(c(1, 10, 100, 1000) * 4464.1,
                                  c(1, 10, 100, 1000) * 126.51,
                                  c(1, 10, 100, 1000) * 14.111,
                                  c(1, 10, 100, 1000) * 15.208,
                                  c(1, 10, 100, 1000) * 17.073) / 60 / 60,
  num_cov = rep(c(2, 3, 4, 5), 5),
  Method = c(rep("ANCOM-BC2", 4),
             rep("MaAsLin2", 4),
             rep("LM-CLR", 4),
             rep("LinDA", 4),
             rep("PALM", 4)),
  dim = c(rep("species K = 92/large", 20)))

  ## K = 92, small
  run_time_ss <- data.frame(time = c(c(1, 10, 100, 1000) * 4708.8,
                                  c(1, 10, 100, 1000) * 122.18,
                                  c(1, 10, 100, 1000) * 13.270,
                                  c(1, 10, 100, 1000) * 13.908,
                                  c(1, 10, 100, 1000) * 2.562) / 60 / 60,
  num_cov = rep(c(2, 3, 4, 5), 5),
  Method = c(rep("ANCOM-BC2", 4),
             rep("MaAsLin2", 4),
             rep("LM-CLR", 4),
             rep("LinDA", 4),
             rep("PALM", 4)),
  dim = c(rep("genus K = 92/small", 20)))


  library("ggplot2")
  library("scales")

  p.ll <- run_time_ll %>%
    ggplot(aes(x=num_cov, y=time, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method)) +
    geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_x_continuous(breaks = 2:5,
                       labels = label_math(expr = 10^.x))  +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    xlab("Number of covariates") +
    ylab("Time (hours)") +  coord_cartesian(ylim = c(0, 50), clip = "off") +  # Prevents truncation
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(1, 10, 10, 10),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position = "none",
          legend.box="vertical",
          axis.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          strip.text = element_text(size = 16),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") +
    guides(colour = guide_legend(ncol = 1))

  p.ls <- run_time_ls %>%
    ggplot(aes(x=num_cov, y=time, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method)) +
    geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_x_continuous(breaks = 2:5,
                       labels = label_math(expr = 10^.x))  +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    xlab("Number of covariates") +
    ylab("Time (hours)") +  coord_cartesian(ylim = c(0, 50), clip = "off") +  # Prevents truncation
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(1, 10, 10, 10),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position = "none",
          legend.box="vertical",
          axis.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          strip.text = element_text(size = 16),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") +
    guides(colour = guide_legend(ncol = 1))


  p.sl <- run_time_sl %>%
    ggplot(aes(x=num_cov, y=time, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method)) +
    geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_x_continuous(breaks = 2:5,
                       labels = label_math(expr = 10^.x))  +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    xlab("Number of covariates") +
    ylab("Time (hours)") +  coord_cartesian(ylim = c(0, 50), clip = "off") +  # Prevents truncation
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(1, 10, 10, 10),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position = "none",
          legend.box="vertical",
          axis.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          strip.text = element_text(size = 16),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") +
    guides(colour = guide_legend(ncol = 1))


  p.ss <- run_time_ss %>%
    ggplot(aes(x=num_cov, y=time, group = Method)) +
    geom_line(linewidth = 0.7, aes(color=Method)) +
    geom_point(size = 2,aes(color = Method), pch = 18) +
    scale_x_continuous(breaks = 2:5,
                       labels = label_math(expr = 10^.x))  +
    scale_color_manual(
      breaks = c("ANCOM-BC2", "LinDA", "LM-CLR",  "MaAsLin2", "PALM"),
      values = c("blue", "skyblue", "orange", "#4dac26", "red")) +
    xlab("Number of covariates") +
    ylab("Time (hours)") +  coord_cartesian(ylim = c(0, 50), clip = "off") +  # Prevents truncation
    theme_bw() +
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(1, 10, 10, 10),
          plot.title = element_blank(),
          axis.title = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = 'white'),
          legend.position = "none",
          legend.box="vertical",
          axis.text = element_text(size = 15),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          strip.text = element_text(size = 16),
          legend.key.width = unit(1.5,"cm")) +
    labs(linetype="Data") +
    guides(colour = guide_legend(ncol = 1))

  pp <- ggpubr::annotate_figure(
    ggpubr::ggarrange(p.ll, p.sl, p.ls, p.ss,
                      nrow = 1, ncol = 4,  common.legend = TRUE, legend = "none", label.y = "AUPRC"),
    bottom = ggpubr::text_grob("Number of covariates", size = 15, face = "bold"),
    left = ggpubr::text_grob("Time (hours)", size = 15, rot = 90, hjust = 0.2, face = "bold"),
  )
