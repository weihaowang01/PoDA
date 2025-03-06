.imputation=function(Y,Z,formulaa){
  
  n=ncol(Y)
  m=nrow(Y)
  ###################### LinDA补值法 ####################################
  adaptive=TRUE
  corr.cut = 0.1
  allvars <- all.vars(as.formula(formulaa))
  names(Z) <- allvars
  Y<-as.data.frame(Y)
  Z<-as.data.frame(Z)
  if (any(Y == 0)) {
    N <- colSums(Y)
    
    if (adaptive) {
      logN <- log(N)
      tmp <- lm(as.formula(paste0("logN", formulaa)),
                Z)
      
      corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr.pval <= corr.cut)) {
        # cat("Imputation approach is used.\n")
        imputation <- TRUE
      }else {
        # cat("Pseudo-count approach is used.\n")
        imputation <- FALSE
      }
    }
    if (imputation) {
      N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
      N.mat[Y > 0] <- 0
      tmp <- N[max.col(N.mat)]
      Y <- Y + N.mat/tmp
    }else {
      Y <- Y + 0.5
    }
  }

  return(as.matrix(Y))
}