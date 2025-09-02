# the imputation method proposed by Zhou et al.
#' @title The imputation method proposed by Zhou et al.

#' @param Y  matrix or data frame; observed abundance matrix. Row: taxa; column: samples. NAs are not expected in observed abundance matrix.
#' @param Z  matrix or data frame of covariates; The rows of \code{Z} correspond to the columns of \code{Y}. NAs are allowed.
#' @param formula character. For example: \code{formula = '~x1*x2+x3'}.


#' @return The observed abundance matrix after value imputation.
imputation=function(Y,Z,formula){
  n=ncol(Y)
  m=nrow(Y)
  adaptive=TRUE
  corr.cut = 0.1
  allvars <- all.vars(as.formula(formula))
  names(Z) <- allvars
  Y<-as.data.frame(Y)
  Z<-as.data.frame(Z)
  if (any(Y == 0)) {
    N <- colSums(Y)

    if (adaptive) {
      logN <- log(N)
      tmp <- lm(as.formula(paste0("logN", formula)),
                Z)

      corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr.pval <= corr.cut)) {
        imputation <- TRUE
      }else {
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

