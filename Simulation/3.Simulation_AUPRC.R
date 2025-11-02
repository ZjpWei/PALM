#================ PALM simulation based on MIDASim (AUPRC) ================#

  # Packages ----
  library("coin")
  library("purrr")
  library('tidyverse')
  library('glmnet')
  library('MicrobiomeStat')
  library("Maaslin2")
  library("MIDASim")
  library("PALM")

  rm(list = ls())
  # Simulation settings ----
  ## Replicate number: integer range from 1 to 100
  ## Parameters

  AUPRC = S = method = batch_type = data_type = Study <- NULL
  for(Setting in c("large", "small")){
    for(tax.type in c("species", "genus")){

      if(Setting == "small"){
        ## Sample size for small setting
        n.sample <- c(20, 40, 60, 80, 100)
        ## Signature effect size
        effect.sz <- 10
      }else if(Setting == "large"){
        ## Sample size for large setting
        n.sample <- c(100, 120, 140, 160, 180)
        ## Signature effect size
        effect.sz <- 5
      }
      L <- length(n.sample)

      ## Replicate number: from 1 to 100
      s <- 1
      ## Signature effect percentage {5%, 10%, 15%, 20%}
      Ka.per <- 0.05
      ## Signature effect direction {0.5, 1}
      pos.pt <- 0.5
      ## Case/control sequence depth unevenness {0, 1}
      mu <- 0
      ## Signature sparsity is fixed to 0.5
      abd.pt <- 0.5

      # Simulate data ----
      ## random seed
      seed <- 2023
      set.seed(s + seed)
      data.loc <- paste0("./Simulation/AUPRC/Sim_Ka", Ka.per, "_Pos", pos.pt, "_mu", mu, "_", s,".Rdata")

      ## Data processing, remove samples with sequence depth < 2000, apply 0.2 prevalence filter.
      load("./Data/CRC_data/data/count.Rdata")
      load("./Data/CRC_data/data/meta.Rdata")
      load("./Data/CRC_data/data/tax.Rdata")

      ## Check data
      if(!all(gsub("[][]", "",unlist(regmatches(colnames(count), gregexpr("\\[.*?\\]",colnames(count))))) == tax$OTUID)){
        stop("species cannot align well between tax and count.\n")
      }
      meta <- as.data.frame(meta)
      study <- c("AT-CRC", "CN-CRC", "DE-CRC", "FR-CRC", "US-CRC")
      rownames(meta) <- meta$Sample_ID
      meta$Group <- factor(meta$Group, level = c("CTR", "CRC"))
      meta$Study <- factor(meta$Study, levels = study)
      sample.id.kp <- names(which(rowSums(count) >= 2000))
      meta <- meta[sample.id.kp,]
      count <- count[sample.id.kp,]
      pre.filter <- colMeans(count != 0) >= 0.2
      count <- count[,pre.filter]
      tax <- tax[pre.filter,]

      ## genus or species data
      if(tax.type == "genus"){
        unique.genus <- unique(tax$genus)
        count.genus <- NULL
        for(gn in unique.genus){
          count.genus <- cbind(count.genus, rowSums(count[,gn == tax$genus,drop=FALSE]))
        }
        colnames(count.genus) <- unique.genus
        ave.prop <- colMeans(count.genus / rowSums(count.genus))
      }else{
        ave.prop <- colMeans(count / rowSums(count))
      }

      ## Detect most abundant taxa and less abundant taxa
      signal.most.abund <- names(which(ave.prop >= 1e-3))
      signal.less.abund <- names(which(ave.prop  < 1e-3))

      ## Format data
      Y.study <- list()
      X.study <- list()
      for(l in 1:L){
        id.set <- meta$Study == study[l]
        condition <- meta$Group[id.set]
        if(tax.type == "genus"){
          Y.study[[l]] <- count.genus[id.set,]
        }else{
          Y.study[[l]] <- count[id.set,]
        }
        X.study[[l]] <- c(condition == 'CRC') + 0
      }

      ## Signal add
      G <- rep(0, ncol(Y.study[[1]]))
      sign.num <- round(Ka.per * length(G))
      signal.m.abd <- sample(signal.most.abund)[1:round(sign.num/2)]
      signal.l.abd <- sample(signal.less.abund)[1:(sign.num - round(sign.num/2))]
      names(G) <- names(ave.prop)
      G[c(signal.m.abd, signal.l.abd)] <- rbernoulli(n = round(Ka.per * length(G)), p = pos.pt) * 2 - 1

      beta <- runif(length(G), min = 1, max = effect.sz)
      names(beta) <- names(ave.prop)
      beta.pos <- beta
      beta.neg <- beta
      beta.pos[G != 1] <- 1
      beta.neg[G != -1] <- 1

      ## sequence depth group
      seq.grp <- rbernoulli(n = 1, p = 0.5)

      ## Original data
      data.rel <- list()
      for(l in 1:L){
        n <- n.sample[l]
        non.filter <- colSums(Y.study[[l]]) > 0
        Y <- (Y.study[[l]])[,non.filter]
        count.ibd.setup = MIDASim.setup(Y, mode = 'nonparametric', n.break.ties = 100)

        ## control data
        sample.id <- sample(1:nrow(Y), size = n/2, replace = TRUE)
        if(!seq.grp){
          lib.size <- rowSums(Y)[sample.id] * (mu + 1)
        }else{
          lib.size <- rowSums(Y)[sample.id]
        }
        mean.avg.prop <- colMeans(Y/rowSums(Y))

        # ## add signal  version 2 (positive on control)
        mean.avg.prop <- mean.avg.prop * beta.neg[non.filter]
        mean.avg.prop <- mean.avg.prop / sum(mean.avg.prop)

        ## formula from GitHub
        obs.sample.1.ct <- rowSums(Y > 0)
        xvar <- log10(rowSums(Y))
        scamfit.non0 = scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
        sample.1.ct = 10^(predict(scamfit.non0, newdata = data.frame(xvar = log10(lib.size) )) )
        n.taxa = ncol(Y)
        input.sample.prop = sample.1.ct/n.taxa
        input.taxa.prop <- count.ibd.setup$taxa.1.prop * (sum(input.sample.prop) * ncol(Y) / (n/2)) / sum(count.ibd.setup$taxa.1.prop)

        ## Simulate data
        count.ibd.modified <- MIDASim.modify(count.ibd.setup,
                                             lib.size = lib.size,
                                             mean.rel.abund = mean.avg.prop,
                                             sample.1.prop = input.sample.prop,
                                             taxa.1.prop = input.taxa.prop)

        simulated.data <- MIDASim(count.ibd.modified)
        control.data <- simulated.data$sim_count

        ## case data
        sample.id <- sample(1:nrow(Y), size = n/2, replace = TRUE)
        if(seq.grp){
          lib.size <- rowSums(Y)[sample.id] * (mu + 1)
        }else{
          lib.size <- rowSums(Y)[sample.id]
        }
        mean.avg.prop <- colMeans(Y/rowSums(Y))

        ## add signal  version 1 (only on case)
        # mean.avg.prop <- mean.avg.prop * (beta ^ G)[non.filter]
        # mean.avg.prop <- mean.avg.prop / sum(mean.avg.prop)

        # ## add signal  version 2 (positive on case)
        mean.avg.prop <- mean.avg.prop * beta.pos[non.filter]
        mean.avg.prop <- mean.avg.prop / sum(mean.avg.prop)

        ## formula from GitHub
        obs.sample.1.ct <- rowSums(Y > 0)
        xvar <- log10(rowSums(Y))
        scamfit.non0 = scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
        sample.1.ct = 10^(predict(scamfit.non0, newdata = data.frame(xvar = log10(lib.size) )) )
        n.taxa = ncol(Y)
        input.sample.prop = sample.1.ct/n.taxa
        input.taxa.prop <- count.ibd.setup$taxa.1.prop * (sum(input.sample.prop) * ncol(Y) / (n/2)) / sum(count.ibd.setup$taxa.1.prop)

        ## Simulate data
        count.ibd.modified <- MIDASim.modify(count.ibd.setup,
                                             lib.size = lib.size,
                                             mean.rel.abund = mean.avg.prop,
                                             sample.1.prop = input.sample.prop,
                                             taxa.1.prop = input.taxa.prop)

        simulated.data <- MIDASim(count.ibd.modified)
        case.data <- simulated.data$sim_count

        ## summarize data
        tmp.data <- matrix(0, nrow = n, ncol = ncol(Y.study[[l]]),
                           dimnames = list(paste0("C",l,"S",1:n), colnames(Y.study[[l]])))
        tmp.data[,non.filter] <- rbind(control.data, case.data)

        case.count <- tmp.data
        X.lst<- rep(c(0,1), each = n/2)

        data.rel[[l]] <- list(Y = case.count, X = X.lst)
      }

      ## Summarize data
      signal.names <- c(signal.l.abd, signal.m.abd)

      #======================================= main analysis ===========================================#
      source("./utility/roc.R")

      ## Align data
      rel.abd <- list()
      covariate.interest <- list()
      study <- NULL
      Y.pool <- NULL
      X.pool <- NULL
      for(l in 1:L){
        study <- c(study, rep(l, length(data.rel[[l]]$X)) )
        Y.pool <- rbind(Y.pool, data.rel[[l]]$Y)
        X.pool <- c(X.pool, data.rel[[l]]$X)
        rel.abd[[paste0("S",l)]] <- data.rel[[l]]$Y
        covariate.interest[[paste0("S",l)]] <- data.frame(disease = data.rel[[l]]$X)
      }

      ## Combine data
      outcome <- X.pool
      feature.table <- Y.pool
      names(outcome) <- rownames(feature.table)
      feature.table = data.frame(t(feature.table))
      meta.data = data.frame(labels = factor(outcome), study = factor(study))
      colnames(feature.table) <- rownames(meta.data)
      feature.ID <- rownames(feature.table)

      ## function use for ANCOM-BC2
      source("./utility/ancombc.R")

      ## ANCOM-BC2 with sensitivity score filter
      AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))
      AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))

      for(l in 1:L){
        tmp.data <- feature.table[,meta.data$study == l]
        ANCOMBC2.model <- ancombc.fun(feature.table = tmp.data[rowSums(tmp.data)!=0,],
                                      meta = meta.data %>% dplyr::filter(study == l),
                                      formula = "labels",
                                      adjust.method = "fdr",
                                      group = NULL,
                                      subject = NULL,
                                      method = "ancombc2")

        AA.est[ANCOMBC2.model$res$taxon[ANCOMBC2.model$res$passed_ss_labels1],l] <-
          ANCOMBC2.model$res$lfc_labels1[ANCOMBC2.model$res$passed_ss_labels1]
        AA.var[ANCOMBC2.model$res$taxon[ANCOMBC2.model$res$passed_ss_labels1],l] <-
          (ANCOMBC2.model$res$se_labels1[ANCOMBC2.model$res$passed_ss_labels1])^2
      }

      ## Combined statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var
      pval.sin <- 1 - pchisq(q.coef, df = 1)
      pval.sin[is.na(pval.sin)] <- 1

      ## ---- Calculate AUPRC
      true.active <- rep(0, length(feature.ID))
      names(true.active) <- feature.ID
      true.active[signal.names] <- 1
      rank.fdr <- rank(pval.sin[feature.ID])
      Paths <- list()
      for(ll in 1:length(rank.fdr)){
        Paths[[ll]] <- rank.fdr <= ll
      }
      hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = TRUE)

      ## Align identity
      AUPRC <- c(AUPRC, hg.pr$AUPRC)
      S <- c(S, s)
      batch_type <- c(batch_type, Setting)
      data_type <- c(data_type, tax.type)
      method <- c(method, "ANCOM-BC2")
      Study <- c(Study, "meta")

      ## MaAsLin2
      AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))
      AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))

      for(l in 1:L){
        data.tmp <- feature.table[,meta.data$study == l]
        data.tmp <- data.tmp[rowSums(data.tmp) > 0,]
        taxa.id <- rownames(data.tmp)
        rownames(data.tmp) <- paste0("Tax", 1:length(taxa.id))
        res <- Maaslin2::Maaslin2(input_data = data.tmp,
                                  input_metadata = meta.data[meta.data$study == l,],
                                  output = "./",
                                  min_abundance = 0,
                                  min_prevalence = 0,
                                  normalization = "TSS",
                                  transform = "AST",
                                  analysis_method = "LM",
                                  max_significance = 1,
                                  random_effects = NULL,
                                  fixed_effects = "labels",
                                  standardize = FALSE,
                                  plot_heatmap = FALSE,
                                  plot_scatter = FALSE)

        res.tmp <- res$result[res$results$metadata == "labels",]
        res.tmp <- res.tmp[order(as.numeric(sapply(strsplit(res.tmp$feature, split = "Tax"), "[[", 2, simplify = TRUE))),]
        res.tmp$feature <- taxa.id

        AA.est[res.tmp$feature,l] <- res.tmp$coef
        AA.var[res.tmp$feature,l] <- (res.tmp$stderr)^2
      }

      ## Combined statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var
      pval.sin <- 1 - pchisq(q.coef, df = 1)
      pval.sin[is.na(pval.sin)] <- 1

      ## ---- Calculate AUPRC
      true.active <- rep(0, length(feature.ID))
      names(true.active) <- feature.ID
      true.active[signal.names] <- 1
      rank.fdr <- rank(pval.sin[feature.ID])
      Paths <- list()
      for(ll in 1:length(rank.fdr)){
        Paths[[ll]] <- rank.fdr <= ll
      }
      hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = TRUE)

      ## Align identity
      AUPRC <- c(AUPRC, hg.pr$AUPRC)
      S <- c(S, s)
      batch_type <- c(batch_type, Setting)
      data_type <- c(data_type, tax.type)
      method <- c(method, "MaAsLin2")
      Study <- c(Study, "meta")

      ## LM-CLR
      ## `Taxon ~ (Intercept) + disease`
      AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))
      AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))

      for(l in 1:L){
        data.tmp <- feature.table[,meta.data$study == l]
        data.tmp <- data.tmp[rowSums(data.tmp) > 0,]
        Z.pool <- log(t(data.tmp) + 0.5)
        clrx_Z <- apply(Z.pool, 2, function(x) x - rowMeans(Z.pool))
        inormal <- apply(clrx_Z, 2, function(d){
          data.d <- data.frame(labels = meta.data$labels[meta.data$study == l], d = d)
          lm.fit <- lm(d ~ labels, data = data.d)
          sum.fit <- summary(lm.fit)
          if("labels1" %in% rownames(sum.fit$coefficients)){
            return(sum.fit$coefficients["labels1", c("Estimate", "Std. Error")])
          }else{
            return(c(NA, NA))
          }
        })
        AA.est[colnames(inormal),l] <- inormal[1,]
        AA.var[colnames(inormal),l] <- (inormal[2,])^2
      }

      ## Combined statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var
      pval.sin <- 1 - pchisq(q.coef, df = 1)
      pval.sin[is.na(pval.sin)] <- 1

      ## ---- Calculate AUPRC
      true.active <- rep(0, length(feature.ID))
      names(true.active) <- feature.ID
      true.active[signal.names] <- 1
      rank.fdr <- rank(pval.sin[feature.ID])
      Paths <- list()
      for(ll in 1:length(rank.fdr)){
        Paths[[ll]] <- rank.fdr <= ll
      }
      hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = TRUE)

      ## Align identity
      AUPRC <- c(AUPRC, hg.pr$AUPRC)
      S <- c(S, s)
      batch_type <- c(batch_type, Setting)
      data_type <- c(data_type, tax.type)
      method <- c(method, "LM-CLR")
      Study <- c(Study, "meta")

      ## LinDA
      AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))
      AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))

      for(l in 1:L){
        data.tmp <- feature.table[,meta.data$study == l]
        data.tmp <- data.tmp[rowSums(data.tmp) > 0,]

        Linda.model <- MicrobiomeStat::linda(feature.dat = data.tmp,
                                             meta.dat = meta.data %>% dplyr::filter(study == l),
                                             formula = '~ labels',
                                             feature.dat.type = "count",
                                             prev.filter = 0,
                                             adaptive = TRUE,
                                             alpha = 0.05)

        AA.est[rownames(Linda.model$output$labels1),l] <- Linda.model$output$labels1$log2FoldChange
        AA.var[rownames(Linda.model$output$labels1),l] <- (Linda.model$output$labels1$lfcSE)^2
      }

      ## Combined statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var
      pval.sin <- 1 - pchisq(q.coef, df = 1)
      pval.sin[is.na(pval.sin)] <- 1

      ## ---- Calculate AUPRC
      true.active <- rep(0, length(feature.ID))
      names(true.active) <- feature.ID
      true.active[signal.names] <- 1
      rank.fdr <- rank(pval.sin[feature.ID])
      Paths <- list()
      for(ll in 1:length(rank.fdr)){
        Paths[[ll]] <- rank.fdr <= ll
      }
      hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = TRUE)

      ## Align identity
      AUPRC <- c(AUPRC, hg.pr$AUPRC)
      S <- c(S, s)
      batch_type <- c(batch_type, Setting)
      data_type <- c(data_type, tax.type)
      method <- c(method, "LinDA")
      Study <- c(Study, "meta")

      ## DESeq2
      AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))
      AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))

      ## DESeq2
      for(l in 1:L){
        objs <- DESeqDataSetFromMatrix(countData = feature.table[,meta.data$study == l],
                                       colData = meta.data %>% dplyr::filter(study == l) %>%
                                         dplyr::transmute(labels = labels),
                                       design = ~ labels)

        DESeq2.model <- DESeq(object = objs, sfType = "poscounts", minReplicatesForReplace = Inf)

        DESeq2.test <- results(DESeq2.model, name = "labels_1_vs_0",
                               independentFiltering = FALSE, cooksCutoff = FALSE)

        AA.est[rownames(DESeq2.test),l] <- DESeq2.test$log2FoldChange
        AA.var[rownames(DESeq2.test),l] <- (DESeq2.test$lfcSE)^2
      }

      ## Combined statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var
      pval.sin <- 1 - pchisq(q.coef, df = 1)
      pval.sin[is.na(pval.sin)] <- 1

      ## ---- Calculate AUPRC
      true.active <- rep(0, length(feature.ID))
      names(true.active) <- feature.ID
      true.active[signal.names] <- 1
      rank.fdr <- rank(pval.sin[feature.ID])
      Paths <- list()
      for(ll in 1:length(rank.fdr)){
        Paths[[ll]] <- rank.fdr <= ll
      }
      hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = TRUE)

      ## Align identity
      AUPRC <- c(AUPRC, hg.pr$AUPRC)
      S <- c(S, s)
      batch_type <- c(batch_type, Setting)
      data_type <- c(data_type, tax.type)
      method <- c(method, "DESeq2")
      Study <- c(Study, "meta")

      ## PALM
      null_obj_RA <- PALM::palm.null.model(rel.abd = rel.abd, prev.filter = 0)

      summary.score.RA <- PALM::palm.get.summary(null.obj = null_obj_RA,
                                                 covariate.interest = covariate.interest)

      ## original PALM corrected
      PALM.model <- palm.meta.summary(summary.stats = summary.score.RA, p.adjust.method = "fdr")

      ## Calculate FDR
      AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))
      AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))

      for(l in 1:L){
        AA.est[PALM.model$disease$feature,l] <- PALM.model$disease[[paste0("S", l,"_effect")]]
        AA.var[PALM.model$disease$feature,l] <- (PALM.model$disease[[paste0("S", l,"_stderr")]])^2
      }

      ## Combined statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var
      pval.sin <- 1 - pchisq(q.coef, df = 1)
      pval.sin[is.na(pval.sin)] <- 1

      ## ---- Calculate AUPRC
      true.active <- rep(0, length(feature.ID))
      names(true.active) <- feature.ID
      true.active[signal.names] <- 1
      rank.fdr <- rank(pval.sin[feature.ID])
      Paths <- list()
      for(ll in 1:length(rank.fdr)){
        Paths[[ll]] <- rank.fdr <= ll
      }
      hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = TRUE)

      ## Align identity
      AUPRC <- c(AUPRC, hg.pr$AUPRC)
      S <- c(S, s)
      batch_type <- c(batch_type, Setting)
      data_type <- c(data_type, tax.type)
      method <- c(method, "PALM")
      Study <- c(Study, "meta")

      ## PALM tuned
      null_obj_RA <- PALM::palm.null.model(rel.abd = rel.abd, prev.filter = 0)

      summary.score.RA <- PALM::palm.get.summary(null.obj = null_obj_RA,
                                                 covariate.interest = covariate.interest,
                                                 correct = FALSE)

      ## PALM tuned
      source("./utility/palm.hom.unadjust.R")
      PALM_tune <- palm.hom.unadjust(summary.stats = summary.score.RA,
                                     output.best.one = TRUE,
                                     verbose = TRUE)

      names(PALM_tune$disease$delta) <- names(summary.score.RA)

      for(d in names(summary.score.RA)) {
        min.delta <- PALM_tune$disease$delta[d]
        non.na <- !is.na(summary.score.RA[[d]]$est[,1])
        summary.score.RA[[d]]$est[non.na, 1] <- summary.score.RA[[d]]$est[non.na,1] + min.delta
        adjust.part <- sum(non.na)/(2 * sum(1/sqrt(2 * pi)/summary.score.RA[[d]]$stderr[non.na, 1]))^2
        summary.score.RA[[d]]$stderr[non.na, 1] <- sqrt((summary.score.RA[[d]]$stderr[non.na, 1])^2 + adjust.part)
      }

      ## original PALM corrected
      PALM.model <- palm.meta.summary(summary.stats = summary.score.RA, p.adjust.method = "fdr")


      ## Calculate FDR
      AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))
      AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))

      for(l in 1:L){
        AA.est[PALM.model$disease$feature,l] <- PALM.model$disease[[paste0("S", l,"_effect")]]
        AA.var[PALM.model$disease$feature,l] <- (PALM.model$disease[[paste0("S", l,"_stderr")]])^2
      }

      ## Combined statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var
      pval.sin <- 1 - pchisq(q.coef, df = 1)
      pval.sin[is.na(pval.sin)] <- 1

      ## ---- Calculate AUPRC
      true.active <- rep(0, length(feature.ID))
      names(true.active) <- feature.ID
      true.active[signal.names] <- 1
      rank.fdr <- rank(pval.sin[feature.ID])
      Paths <- list()
      for(ll in 1:length(rank.fdr)){
        Paths[[ll]] <- rank.fdr <= ll
      }
      hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = TRUE)

      ## Align identity
      AUPRC <- c(AUPRC, hg.pr$AUPRC)
      S <- c(S, s)
      batch_type <- c(batch_type, Setting)
      data_type <- c(data_type, tax.type)
      method <- c(method, "PALM.tuned")
      Study <- c(Study, "meta")

      ## PALM Wald
      source("./utility/PALM_wald.R")
      summary.score.RA <- melody.get.summary.wald.poi(rel.abd = rel.abd,
                                                      covariate.interest = covariate.interest,
                                                      prev.filter = 0)

      ## original PALM corrected
      ## Calculate FDR
      AA.est <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))
      AA.var <- matrix(NA, nrow = length(feature.ID), ncol = L,
                       dimnames = list(feature.ID, as.character(1:L)))

      for(l in 1:L){
        AA.est[rownames(summary.score.RA[[l]]$est),l] <- summary.score.RA[[l]]$est
        AA.var[rownames(summary.score.RA[[l]]$est),l] <- (summary.score.RA[[l]]$stderr)^2
      }

      ## Combined statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var
      pval.sin <- 1 - pchisq(q.coef, df = 1)
      pval.sin[is.na(pval.sin)] <- 1

      ## ---- Calculate AUPRC
      true.active <- rep(0, length(feature.ID))
      names(true.active) <- feature.ID
      true.active[signal.names] <- 1
      rank.fdr <- rank(pval.sin[feature.ID])
      Paths <- list()
      for(ll in 1:length(rank.fdr)){
        Paths[[ll]] <- rank.fdr <= ll
      }
      hg.pr <- huge.pr(Paths, as.vector(true.active), verbose=FALSE, plot = TRUE)

      ## Align identity
      AUPRC <- c(AUPRC, hg.pr$AUPRC)
      S <- c(S, s)
      batch_type <- c(batch_type, Setting)
      data_type <- c(data_type, tax.type)
      method <- c(method, "PALM.wald")
      Study <- c(Study, "meta")

    }
  }

  PRC <- data.frame(AUPRC = AUPRC, S = S, method = method, Study, Settings = batch_type, tax.type = data_type)

  save(PRC, file = data.loc)
