################### formula for multiple-test ###################
palm.get.summary <- function(null.obj,
                             covariate.interest,
                             cluster = NULL,
                             parallel.core = NULL,
                             verbose = FALSE) {

  #=== Check input data ===#
  study.ID <- names(null.obj)
  feature.ID <- NULL
  cov.names <- NULL
  for(d in names(null.obj)){
    feature.ID <- c(feature.ID, colnames(null.obj[[d]]$Y_I))
    if(is.null(colnames(covariate.interest[[d]]))){
      stop("covariate.interest is not a matrix with column names. Please check your data.")
    }else{
      covariate.interest[[d]] <- data.frame(covariate.interest[[d]])
      cov.names <- c(cov.names, colnames(covariate.interest[[d]]))
    }
  }
  feature.ID <- unique(feature.ID)
  cov.names <- unique(cov.names)
  K <- length(feature.ID)
  if(verbose){
    message.num <- 0
  }

  ## Match rel.data, covariate.interest and covariate.adjust.
  SUB.id <- list()
  for(d in study.ID){
    if(is.null(covariate.interest[[d]])){
      stop("Study IDs in rel.data and covariate.interest don't match, please check the input data.")
    }
    if(!is.data.frame(covariate.interest[[d]])){
      stop("covariate.interest is not a list of data frames.\n")
    }
    if(nrow(covariate.interest[[d]]) != (nrow(null.obj[[d]]$Y_I) + length(null.obj[[d]]$rm.sample.idx))){
      stop("The sample size of covariate.interest is not correct, please check the input data.")
    }else{
      if(length(null.obj[[d]]$rm.sample.idx) > 0){
        covariate.interest[[d]] <- covariate.interest[[d]] %>% dplyr::slice(-null.obj[[d]]$rm.sample.idx)
      }
      rownames(covariate.interest[[d]]) <- rownames(null.obj[[d]]$Y_I)
    }
    if(is.null(cluster[[d]])){
      cluster.nm <- rownames(null.obj[[d]]$Y_I)
      names(cluster.nm) <- rownames(null.obj[[d]]$Y_I)
      SUB.id[[d]] <- cluster.nm
    }else{
      if(!is.vector(cluster[[d]])){
        stop("cluster is not a list of vectors. \n")
      }
      if(length(cluster[[d]]) != (nrow(null.obj[[d]]$Y_I) + length(null.obj[[d]]$rm.sample.idx))){
        stop("The sample size of cluster is not correct, please check the input data. \n")
      }else{
        if(length(null.obj[[d]]$rm.sample.idx) > 0){
          cluster[[d]] <- (cluster[[d]])[-null.obj[[d]]$rm.sample.idx]
        }
      }
      cluster.nm <- cluster[[d]]
      names(cluster.nm) <- rownames(null.obj[[d]]$Y_I)
      SUB.id[[d]] <- cluster.nm
    }
  }

  #=== Sample combination in each study corresponded to covariate.interest ===#
  Sample.info <- lapply(covariate.interest, function(d){
    if(ncol(d) > 1){
      tmp.cov <- apply(d, 2, function(r){which(is.na(r))})
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
    est.mat <- matrix(NA, nrow = K, ncol = length(cov.int.nm),
                      dimnames = list(feature.ID, cov.int.nm))
    cov.mat <- matrix(NA, nrow = K, ncol = length(cov.int.nm),
                      dimnames = list(feature.ID, cov.int.nm))

    for(cov.name in cov.int.nm){
      if(verbose){
        message("++ Construct summary statistics for study ", d, " and covariate of interest ", cov.name, ". ++")
      }
      l <- which(unlist(lapply(cov.int.lst, function(r){cov.name %in% r})))
      rm.sample.id <- cov.int.id[[l]]
      if(length(rm.sample.id) == 0){
        Y_R <- null.obj[[d]]$Y_R
        Y_I <- null.obj[[d]]$Y_I
        X <- null.obj[[d]]$Z
        SUBid <- as.character(SUB.id[[d]])
      }else{
        Y_R <- null.obj[[d]]$Y_R[-rm.sample.id,,drop=FALSE]
        Y_I <- null.obj[[d]]$Y_I[-rm.sample.id,,drop=FALSE]
        X <- null.obj[[d]]$Z[-rm.sample.id,,drop=FALSE]
        SUBid <- as.character(SUB.id[[d]])[-rm.sample.id]
      }
      uniq.SUBid <- unique(SUBid)
      X.names <- colnames(X)
      if(!is.numeric(covariate.interest[[d]][,cov.name])){
        stop("Covariate.interest should be numeric, please check your input data.\n")
      }
      X[,ncol(X)] <- covariate.interest[[d]][rownames(X), cov.name]

      ## loop for summary statistics
      ests <- NULL
      covs <- NULL
      for(k in 1:ncol(Y_I)){
        name.cov <- paste0(X.names, ":", rep(k, ncol(X)))
        s.i.lst <- NULL
        I_mat <- 0
        for(l in 1:length(uniq.SUBid)){
          s.i.SUB <- 0
          for(i in which(uniq.SUBid[l] == SUBid)){
            s.i.SUB <- s.i.SUB + t(Y_R[i,k] * X[i,])
            I_mat <- I_mat + Y_I[i,k] * X[i,] %*% t(X[i,])
          }
          s.i.lst <- rbind(s.i.lst, s.i.SUB)
        }
        colnames(I_mat) <- name.cov
        rownames(I_mat) <- name.cov
        colnames(s.i.lst) <- name.cov

        ## Get estimate and variance
        beta.name <- paste0(X.names[length(X.names)],":", k)
        gamma.name <- setdiff(name.cov, beta.name)
        I_beta <- I_mat[beta.name, beta.name] - I_mat[beta.name, gamma.name]%*%solve(I_mat[gamma.name, gamma.name])%*%I_mat[gamma.name,beta.name]
        ests <- c(ests, solve(I_beta) %*% colSums(s.i.lst)[beta.name])
        core.U <- 0
        for(i in 1:nrow(s.i.lst)){
          tmp.U <- s.i.lst[i,beta.name] - I_mat[beta.name, gamma.name]%*%solve(I_mat[gamma.name, gamma.name]) %*% s.i.lst[i,gamma.name]
          core.U <- core.U + tmp.U %*% t(tmp.U)
        }
        covs <- c(covs, solve(I_beta) %*% (core.U) %*% solve(I_beta))

      }
      est.mat[colnames(Y_I), cov.name] <- ests
      cov.mat[colnames(Y_I), cov.name] <- covs
    }
    summary.stat.study.one <- list(est = est.mat, var = cov.mat, n = nrow(Y_I))

    #=== output ===#
    summary.stat.study[[d]] <- summary.stat.study.one
  }
  return(summary.stat.study)
}
