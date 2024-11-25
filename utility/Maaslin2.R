Maaslin2_meta <- function(feature.table,
                          meta.data,
                          fixed_effects,
                          random_effects,
                          study = study,
                          rma.method = "EE"){

  study.id <- unique(meta.data[[study]])

  Maaslin2.lst <- list()
  for(l in study.id){
    data.tmp <- feature.table[,meta.data[[study]] == l]
    data.tmp <- data.tmp[rowSums(data.tmp) > 0,]
    taxa.id <- rownames(data.tmp)
    rownames(data.tmp) <- paste0("Tax", 1:length(taxa.id))
    res <- Maaslin2::Maaslin2(input_data = data.tmp,
                              input_metadata = meta.data[meta.data[[study]] == l,],
                              output = "./",
                              min_abundance = 0,
                              min_prevalence = 0,
                              normalization = "TSS",
                              transform = "AST",
                              analysis_method = "LM",
                              max_significance = 1,
                              random_effects = random_effects,
                              fixed_effects = fixed_effects,
                              standardize = FALSE,
                              plot_heatmap = FALSE,
                              plot_scatter = FALSE)

    res.tmp <- res$result[res$results$metadata == "labels",]
    res.tmp <- res.tmp[order(as.numeric(sapply(strsplit(res.tmp$feature, split = "Tax"), "[[", 2, simplify = TRUE))),]
    res.tmp$feature <- taxa.id
    Maaslin2.lst[[l]] <- res.tmp
  }

  meta.fit <- NULL
  meta.pval <- NULL
  het.pval <- NULL
  Q <- NULL
  if(length(study.id) > 0){
    for(l in rownames(feature.table)){
      ests <- NULL
      vars <- NULL
      for(s in study.id){
        if(l %in% Maaslin2.lst[[s]]$feature){
          ests <- c(ests, Maaslin2.lst[[s]]$coef[l == Maaslin2.lst[[s]]$feature])
          vars <- c(vars, (Maaslin2.lst[[s]]$stderr[l == Maaslin2.lst[[s]]$feature])^2)
        }
      }
      m <- metafor::rma(yi = ests, vi = vars, method = rma.method)
      meta.fit <- c(meta.fit, m$beta)
      meta.pval <- c(meta.pval, m$pval)
      het.pval <- c(het.pval, m$QEp)
      Q <- c(Q, m$QE)
    }
    meta.lst <- data.frame(feature = rownames(feature.table),
                           coef = meta.fit,
                           pval = meta.pval,
                           het.pval = het.pval,
                           Q = Q)
  }

  result.all <- list(Maaslin2.lst = Maaslin2.lst, meta.lst = meta.lst)
}




