melody.get.summary.wald.poi <- function(rel.abd,
                                        covariate.interest,
                                        covariate.adjust = NULL,
                                        cluster = NULL,
                                        depth.filter = 0,
                                        prev.filter = 0.1,
                                        depth = NULL,
                                        correct = TRUE,
                                        verbose = FALSE) {

  #=== Check input data ===#
  study.ID <- names(rel.abd)
  if(is.null(study.ID)){
    stop("Please check the study name in rel.data.")
  }else{
    if(length(study.ID) != length(unique(study.ID))){
      stop("Multiple studies has the same names in rel.data, please check the input data.")
    }
  }

  ### match rel.data and covariate.adjust
  for(d in study.ID){
    if(!all(!is.na(rel.abd[[d]]))){
      stop("Detect NA in relative abundant counts.\n")
    }
    if(min(rel.abd[[d]]) < 0){
      stop("Detect negative value in relative abundant counts.\n")
    }
    if(!is.null(covariate.adjust[[d]])){
      if(!is.data.frame(covariate.adjust[[d]])){
        stop("covariate.adjust is not a list of data frames.\n")
      }
      if(nrow(covariate.adjust[[d]]) != nrow(rel.abd[[d]])){
        stop("The sample size doesn't match between rel.abd and covariate.adjust, please check the input data.")
      }else{
        rownames(covariate.adjust[[d]]) <- rownames(rel.abd[[d]])
      }
      if(!all(!is.na(covariate.adjust[[d]]))){
        stop("Detect NA in covariate.adjust.\n")
      }
    }
  }

  #=== Filter the samples using depth.filter ===#
  rm.sample.idx <- list()
  for(d in study.ID){
    depth.kp <- which(rowSums(rel.abd[[d]]) > depth.filter)
    rm.sample.idx[[d]] <- which(rowSums(rel.abd[[d]]) <= depth.filter)
    rel.abd[[d]] <- rel.abd[[d]][depth.kp,]
    covariate.interest[[d]] <- data.frame(covariate.interest[[d]]) %>% dplyr::slice(depth.kp)
    if(!is.null(covariate.adjust[[d]])){
      covariate.adjust[[d]] <- data.frame(covariate.adjust[[d]]) %>% dplyr::slice(depth.kp)
    }
  }

  #=== Match samples in relative abundant counts and sample data ===#
  dat <- list()
  for(d in study.ID){
    Y.pool <- rel.abd[[d]]
    if(length(Y.pool) == 0){
      warning(paste0("Less than 20 samples in study ", d, ", the summary statistics is not stable.",
                     " Remove study ", d, "\n"))
    }else if(nrow(Y.pool) < 20){
      warning(paste0("Less than 20 samples in study ", d, ", the summary statistics is not stable.",
                     " Remove study ", d, "\n"))
    }else{
      X.pool <- matrix((covariate.interest[[d]])[[1]], nrow = nrow(Y.pool))
      if(is.null(covariate.adjust)){
        dat[[d]] <- list(Y = Y.pool, X = X.pool)
      }else{
        cov.adjust <- covariate.adjust[[d]]
        for(cov_name in colnames(cov.adjust)){
          if(all(!is.na(cov.adjust[[cov_name]]))){
            if(is.factor(cov.adjust[[cov_name]])){
              class( cov.adjust[cov_name])
              dummys <- as.data.frame(model.matrix(formula(paste("~", cov_name)), data = cov.adjust[cov_name]))
              X.pool <- as.matrix(cbind(dummys[,-1], X.pool))
            }else if(is.character(cov.adjust[[cov_name]])){
              dummys <- as.data.frame(model.matrix(formula(paste("~", cov_name)), data = cov.adjust[cov_name]))
              X.pool <- as.matrix(cbind(dummys[,-1], X.pool))
            }else{
              X.pool <- cbind(cov.adjust[[cov_name]], X.pool)
            }
          }else{
            warning(paste0("NA presents in `", cov_name, "`, please check sample.data."))
          }
        }
        dat[[d]] <- list(Y = Y.pool, X = X.pool)
      }
    }
  }

  #=== check maximum C.V. 1.5 output warning if too large ===#
  feature.ID <- NULL
  for(d in study.ID){
    feature.ID <- c(feature.ID, colnames(dat[[d]]$Y))
  }
  feature.ID <- sort(unique(feature.ID))
  for(d in study.ID){
    dat[[d]]$Y <- dat[[d]]$Y[,sort(intersect(feature.ID, colnames(dat[[d]]$Y)))]
  }
  K <- length(feature.ID)

  #=== Switch the reference to the last column ===#
  null.obj <- list()
  for(d in study.ID){
    Y.sub <- dat[[d]]$Y
    X.sub <- cbind(1, dat[[d]]$X)
    colnames(X.sub) <- c("Intercept", paste0("V_", as.character(1:(ncol(X.sub)-1))))
    rownames(X.sub) <- rownames(Y.sub)
    feature.set.tmp <- colSums(Y.sub != 0) > prev.filter
    feature.set.tmp[ncol(Y.sub)] <- TRUE
    Y.sub.tmp <- Y.sub[,feature.set.tmp]

    ## Check if any correlations are 1.
    cors <- cor(Y.sub.tmp)
    lo_tri <- lower.tri(cors, diag = TRUE)
    cors[lo_tri] <- 0
    # if(!all(abs(cors) != 1)){
    #   row_col <- which(cors == 1, arr.ind = TRUE)
    #   row_id <- unique(row_col[,"row"])
    #   col_id <- unique(row_col[,"col"])
    #   tax.rm <- colnames(Y.sub)[row_id]
    #   tax.kp <- colnames(Y.sub)[col_id]
    #   feature.set.tmp[tax.rm] <- FALSE
    #   warning("Some features have high correlation in study ", d, ", Rmove features:\n",
    #           paste0("    ",tax.rm,"\n"))
    #
    # }
    Y.sub <- Y.sub[,feature.set.tmp]
    X.idx <- ncol(X.sub)
    n = nrow(Y.sub)
    K = ncol(Y.sub)
    dm = ncol(X.sub)

    if(dm == 2){
      est.single <- matrix(NA, nrow = 2, ncol = K, dimnames = list(c(":(Intercept)", ":X"), colnames(Y.sub)))
    }else{
      est.single <- matrix(NA, nrow = d, ncol = K, dimnames = list(c(":(Intercept)", paste0(":XV_", 1:(dm-1))), colnames(Y.sub)))
    }
    Y_b <- matrix(0, nrow = nrow(Y.sub), ncol = ncol(Y.sub), dimnames = list(rownames(Y.sub), colnames(Y.sub)))

    #=== Score-test summary statistics ===#
    if(is.null(depth)){
      N <- rowSums(Y.sub)
    }else{
      N <- (depth[[d]])[rownames(Y.sub)]
    }

    suppressWarnings(
      for(k in 1:K){
        # Try brglmFit model
        input.data.tmp = list(Y = Y.sub[,k], X = X.sub[,-1])
        glm.out.tmp = glm(Y ~ X, data = input.data.tmp,
                          family = poisson(link = "log"),
                          offset = log(N),
                          method = brglm2::brglmFit,
                          type = "AS_mean")

        names(glm.out.tmp$coefficients) <- paste0(":", names(glm.out.tmp$coefficients))
        est.single[names(glm.out.tmp$coefficients),k] <- glm.out.tmp$coefficients

        ## compute the bias part
        Qmat <- qr.Q(glm.out.tmp$qr)
        Y_b[,k] <- 0.5 * rowSums(Qmat * Qmat)
      }
    )
    tmp <- list(est = est.single, n = nrow(Y.sub))


    #=== summary null model ===#
    est <- t(est.single)
    pp_mat <- NULL
    for(i in 1:length(N)){
      dd <- colSums(matrix(rep(X.sub[i,], nrow(est)), nrow = dm) * t(est), na.rm = TRUE)
      pp <- exp(dd)
      pp_mat <- rbind(pp_mat, pp)
    }
    rownames(pp_mat) <- rownames(Y.sub)
    null.obj.one <- list(est = est, p = pp_mat, Y_b = Y_b,
                         N = N, Y = Y.sub, X = X.sub, para.id = d)

    #=== output ===#
    null.obj[[d]] <- null.obj.one
  }

  #=== reorder output ===#
  for(ll in 1:length(null.obj)){
    null.obj[[ll]]$para.id <- NULL
  }

  for(d in study.ID){
    null.obj[[d]]$rm.sample.idx <- rm.sample.idx[[d]]
  }

  #=== Check input data ===#
  study.ID <- names(null.obj)
  for(d in names(null.obj)){
    if(is.null(colnames(covariate.interest[[d]]))){
      stop("covariate.interest is not a matrix with column names. Please check your data.")
    }else{
      covariate.interest[[d]] <- data.frame(covariate.interest[[d]])
    }
  }

  cov.names <- NULL
  for(d in study.ID){
    cov.names <- c(cov.names, colnames(covariate.interest[[d]]))
  }
  cov.names <- unique(cov.names)
  if(verbose){
    message.num <- 0
  }

  ### match rel.data, covariate.interest and covariate.adjust
  SUB.id <- list()
  for(d in study.ID){
    if(is.null(covariate.interest[[d]])){
      stop("Study IDs in rel.data and covariate.interest don't match, please check the input data.")
    }
    if(!is.data.frame(covariate.interest[[d]])){
      stop("covariate.interest is not a list of data frames.\n")
    }
    if(nrow(covariate.interest[[d]]) != (length(null.obj[[d]]$N) + length(null.obj[[d]]$rm.sample.idx))){
      stop("The sample size of covariate.interest is not correct, please check the input data.")
    }else{
      if(length(null.obj[[d]]$rm.sample.idx) > 0){
        covariate.interest[[d]] <- covariate.interest[[d]] %>% dplyr::slice(-null.obj[[d]]$rm.sample.idx)
      }
      rownames(covariate.interest[[d]]) <- names(null.obj[[d]]$N)
    }

    if(is.null(cluster[[d]])){
      cluster.nm <- names(null.obj[[d]]$N)
      names(cluster.nm) <- names(null.obj[[d]]$N)
      SUB.id[[d]] <- cluster.nm
    }else{
      if(!is.vector(cluster[[d]])){
        stop("cluster is not a list of vectors. \n")
      }
      if(length(cluster[[d]]) != (length(null.obj[[d]]$N) + length(null.obj[[d]]$rm.sample.idx))){
        stop("The sample size of cluster is not correct, please check the input data. \n")
      }else{
        if(length(null.obj[[d]]$rm.sample.idx) > 0){
          cluster[[d]] <- (cluster[[d]])[-null.obj[[d]]$rm.sample.idx]
        }
      }
      cluster.nm <- cluster[[d]]
      names(cluster.nm) <- names(null.obj[[d]]$N)
      SUB.id[[d]] <- cluster.nm
    }
  }

  #=== Sample combination in each study corresponded to covariate.interest
  Sample.info <- lapply(covariate.interest, function(d){
    if(ncol(d) > 1){
      tmp.cov <- apply(d,2,function(r){which(is.na(r))})
    }else{
      tmp.cov <- list(which(is.na(d)))
      names(tmp.cov) <- colnames(d)
    }
    if(length(tmp.cov) != 0){
      uni.tmp.cov <- unique(tmp.cov)
      rm.sample.idx <- list()
      rm.sample.cov <- list()
      for(l in 1:length(uni.tmp.cov)){
        rm.sample.idx[[l]] <- uni.tmp.cov[[l]]
        rm.sample.cov[[l]] <- names(tmp.cov)[tmp.cov %in% uni.tmp.cov[l]]
      }
      return(list(rm.sample.cov = rm.sample.cov, rm.sample.idx = rm.sample.idx))
    }else{
      return(list(rm.sample.cov = list(colnames(d)), rm.sample.idx = list(integer())))
    }
  })

  #=== Generate summary statistics ===#
  if(verbose){
    message("++ Construct summary statistics. ++")
  }

  #=== Get summary statistics with given covariate.interest and ridge regularization on covariate matrix ===#
  summary.stat.study <- list()
  for(d in study.ID){
    cov.int.lst <- Sample.info[[d]]$rm.sample.cov
    cov.int.id <- Sample.info[[d]]$rm.sample.idx
    cov.int.nm <- colnames(covariate.interest[[d]])
    feat.id <- colnames(null.obj[[d]]$p)

    est.mat <- matrix(NA, nrow = length(feature.ID), ncol = length(cov.int.nm),
                      dimnames = list(feature.ID, cov.int.nm))
    cov.mat <- matrix(NA, nrow = length(feature.ID), ncol = length(cov.int.nm),
                      dimnames = list(feature.ID, cov.int.nm))

    for(cov.name in cov.int.nm){
      if(verbose){
        message("++ Construct summary statistics for study ", d, " and covariate of interest ", cov.name, ". ++")
      }
      l <- which(unlist(lapply(cov.int.lst, function(r){cov.name %in% r})))
      rm.sample.id <- cov.int.id[[l]]
      if(length(rm.sample.id) == 0){
        est_mat <- null.obj[[d]]$est
        Y_b <- null.obj[[d]]$Y_b
        pp_mat <- null.obj[[d]]$p
        Y.sub <- null.obj[[d]]$Y
        X.sub <- null.obj[[d]]$X
        N <- null.obj[[d]]$N
        SUBid <- as.character(SUB.id[[d]])
      }else{
        est_mat <- null.obj[[d]]$est[-rm.sample.id,]
        Y_b <- null.obj[[d]]$Y_b[-rm.sample.id,]
        pp_mat <- null.obj[[d]]$p[-rm.sample.id,]
        Y.sub <- null.obj[[d]]$Y[-rm.sample.id,]
        X.sub <- null.obj[[d]]$X[-rm.sample.id,]
        N <- null.obj[[d]]$N[-rm.sample.id]
        SUBid <- as.character(SUB.id[[d]])[-rm.sample.id]
      }
      K <- ncol(Y.sub)
      uniq.SUBid <- unique(SUBid)
      X.name <- colnames(X.sub)
      if(!is.numeric(covariate.interest[[d]][,cov.name])){
        stop("Covariate.interest should be numeric, please check your input.")
      }
      X.sub[,ncol(X.sub)] <- covariate.interest[[d]][rownames(X.sub),cov.name]

      covs <- NULL
      for(k in 1:K){
        name.cov <- paste0(X.name, ":", as.character(rep(k,ncol(X.sub))))
        s.i.lst <- NULL
        I_mat <- 0
        for(l in 1:length(uniq.SUBid)){
          s.i.SUB <- 0
          for(i in which(uniq.SUBid[l] == SUBid)){
            mu_hat <- pp_mat[i,k] * N[i]
            s.i.SUB <- s.i.SUB + t((Y.sub[i,k] - mu_hat + Y_b[i,k]) * X.sub[i,])
            I_mat <- I_mat + mu_hat * X.sub[i,] %*% t(X.sub[i,])
          }
          s.i.lst <- rbind(s.i.lst, s.i.SUB)
        }
        colnames(I_mat) <- name.cov
        rownames(I_mat) <- name.cov
        colnames(s.i.lst) <- name.cov

        ## Get est and cov
        cov_R <- solve(I_mat)
        colnames(cov_R) <- name.cov
        rownames(cov_R) <- name.cov
        colnames(s.i.lst) <- name.cov

        #=== solve GEE equation ===#
        sand <- 0
        nn <- length(uniq.SUBid)
        for(ii in 1:nn){
          sand  <- sand + (cov_R %*% s.i.lst[ii,]) %*% t(cov_R %*% s.i.lst[ii,])
        }
        Sigma_lambda <- sand / nn
        cov.ridge <- Sigma_lambda * nn
        colnames(cov.ridge) <- paste0(X.name,":", rep(k, each = length(X.name)))
        rownames(cov.ridge) <- paste0(X.name,":", rep(k, each = length(X.name)))
        cov_mat <- cov.ridge[paste0(X.name[length(X.name)],":", k),
                             paste0(X.name[length(X.name)],":", k)]
        covs <- c(covs, cov_mat)
      }
      cov.mat[colnames(Y.sub), cov.name] <- covs

      if(ncol(est_mat) == 2){
        est.mat[rownames(est_mat), cov.name] <- est_mat[, ":X"]
      }else{
        est.mat[rownames(est_mat), cov.name] <- est_mat[, paste0(":XV_", ncol(est_mat)-1)]
      }
    }

    summary.stat.study.one <- list(est = est.mat, stderr = sqrt(cov.mat), n = length(N))

    #=== output ===#
    summary.stat.study[[d]] <- summary.stat.study.one
  }

  ## Correct summary statistics
  if(correct){
    covariate.interest <- unique(unlist(lapply(summary.stat.study, function(d){colnames(d$est)})))
    for(cov.int in covariate.interest){
      for(d in study.ID){
        if(cov.int %in% colnames(summary.stat.study[[d]]$est)){
          ## adjust estimates
          min.delta <- median( - summary.stat.study[[d]]$est[,cov.int], na.rm = TRUE)
          non.na <- !is.na(summary.stat.study[[d]]$est[,cov.int])
          summary.stat.study[[d]]$est[non.na,cov.int] <- summary.stat.study[[d]]$est[non.na,cov.int] + min.delta

          ## adjust variances
          adjust.part <- sum(non.na) / (2 * sum(1 / sqrt(2*pi) / summary.stat.study[[d]]$stderr[non.na,cov.int]))^2
          summary.stat.study[[d]]$stderr[non.na,cov.int] <- sqrt((summary.stat.study[[d]]$stderr[non.na,cov.int])^2 + adjust.part)
        }
      }
    }
  }


  return(summary.stat.study)
}

