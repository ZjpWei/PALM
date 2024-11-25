palm.test <- function(summary.stats = summary.stats,
                      data.type,
                      p.adjust.method = "fdr",
                      cor.method = "original"){
  palm.meta <- list()
  covariate.interest <- unique(unlist(lapply(summary.stats, function(d){colnames(d$est)})))
  for(cov.int in covariate.interest){
    palm_fits <- list()
    study.ID <- NULL
    feature.ID <- NULL
    for(d in names(summary.stats)){
      if(cov.int %in% colnames(summary.stats[[d]]$est)){
        study.ID <- c(study.ID, d)
        feature.ID <- c(feature.ID, rownames(summary.stats[[d]]$est))
      }
    }
    feature.ID <- unique(feature.ID)

    AA.est <- matrix(NA, nrow = length(feature.ID), ncol = length(study.ID),
                     dimnames = list(feature.ID, as.character(study.ID)))
    AA.var <- matrix(NA, nrow = length(feature.ID), ncol = length(study.ID),
                     dimnames = list(feature.ID, as.character(study.ID)))

    for(d in study.ID){
      if(data.type == "AA"){
        min.delta <- 0
      }else if(data.type == "RA"){
        min.delta <- median(- summary.stats[[d]]$est[,cov.int], na.rm = TRUE)
      }else{
        stop("The data type should only be `AA` or `RA`, please check your input data.\n")
      }
      non.na <- !is.na(summary.stats[[d]]$est[,cov.int])
      median.var <- 0
      beta.coef <- summary.stats[[d]]$est[non.na,cov.int] + min.delta

      if(cor.method == "original"){
        pval <- 1 - pchisq(beta.coef^2 / summary.stats[[d]]$var[non.na,cov.int], df = 1)
        qval <- p.adjust(p = pval, method = p.adjust.method)
        palm_fits[[d]] <- data.frame(feature = names(beta.coef),
                                     coef = beta.coef,
                                     stderr = sqrt(summary.stats[[d]]$var[non.na,cov.int]),
                                     pval = pval,
                                     qval = qval)
      }else{
        pval <- 1 - pchisq(beta.coef^2 / (summary.stats[[d]]$var[non.na,cov.int] + sum(summary.stats[[d]]$var[non.na,cov.int])/sum(non.na)^2), df = 1)
        qval <- p.adjust(p = pval, method = p.adjust.method)
        palm_fits[[d]] <- data.frame(feature = names(beta.coef),
                                     coef = beta.coef,
                                     stderr = sqrt(summary.stats[[d]]$var[non.na,cov.int] + sum(summary.stats[[d]]$var[non.na,cov.int])/sum(non.na)^2),
                                     pval = pval,
                                     qval = qval)
      }


      AA.est[names(beta.coef),d] <- beta.coef
      if(cor.method == "original"){
        AA.var[names(beta.coef),d] <- summary.stats[[d]]$var[non.na,cov.int]
      }else{
        AA.var[names(beta.coef),d] <- summary.stats[[d]]$var[non.na,cov.int] + sum(summary.stats[[d]]$var[non.na,cov.int])/sum(non.na)^2
      }
    }

    if(length(study.ID) > 1){
      ## Meta statistics
      est.statics <- rowSums(AA.est / AA.var, na.rm = TRUE)
      var.statics <- rowSums(1 / AA.var, na.rm = TRUE)
      meta.coef <- est.statics / var.statics
      meta.var <- 1 / var.statics
      q.coef <- (meta.coef)^2 / meta.var

      pval.sin <- 1 - pchisq(q.coef, df = 1)
      qval.sin <- p.adjust(pval.sin, method = p.adjust.method)

      meta_fits <- data.frame(feature = rownames(AA.est),
                              coef = meta.coef,
                              stderr = sqrt(var.statics),
                              pval = pval.sin,
                              qval = qval.sin)

      palm.meta[[cov.int]] <- list(meta_fits = meta_fits, palm_fits = palm_fits)
    }else{
      palm.meta[[cov.int]] <- palm_fits[[study.ID]]
    }
  }
  return(palm.meta)
}
