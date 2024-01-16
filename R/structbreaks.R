

#' @title Ordinary Least Squares Parameter Vector Estimation
#' 
#' @description
#' The function `olsqr` estimates the parameters of a linear regression model
#' using Ordinary Least Squares (OLS) method.
#' 
#' @details
#' This function is an adaptation of the code originally written by Yohei Yamamoto
#' and Pierre Perron for MATLAB. The original codes can be found on Pierre Perron's
#' website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @param y Numeric vector: the response variable.
#' @param x Numeric matrix: the matrix of explanatory variables.
#' 
#' @return
#' A numeric matrix representing the estimated parameters of the linear regression model.
#' 
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @examples
#' # Example usage:
#' # Set up a simple linear regression problem
#' set.seed(123)
#' y <- rnorm(100)
#' x <- cbind(1, rnorm(100))
#' 
#' # Estimate parameters using olsqr
#' parameters <- olsqr(y, x)
#' 
#' @keywords internal
#' 
#' @export
olsqr <- function(y, x){
  b <- solve(t(x)%*%x)%*%t(x)%*%y 
  return(as.matrix(b))
}



#' @title Estimate OLS with White Robust Standard Errors
#' 
#' @description
#' This function estimates an ordinary least squares model with White robust standard errors.
#' 
#' @details
#' Note: This code is an adaptation of the one originally written by Yohei Yamamoto
#' and Pierre Perron for MATLAB. Original codes can be found on Pierre Perron's
#' website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @param Y Numeric vector:  the dependent variable.
#' @param X Numeric matrix:  the regressors.
#' @param add_constant Logical: indicating whether to add a constant to X. Default is FALSE.
#' 
#' @return
#' A list containing the following components:
#' -`y`: Dependent variable (numeric vector).
#' -`w`: Regressors (numeric matrix).
#' -`coef`: Estimator corresponding to the k regressors (numeric vector).
#' -`SE`: Standard error of the estimator (numeric vector).
#' -`SE_robust`: White robust standard error of the estimator (numeric vector).
#' -`sigma2`: Estimated variance of disturbances (numeric).
#' -`resid`: Residuals series of the regression (numeric vector).
#' -`SSR`: Sum of squared residuals (numeric).
#' -`t_stat`: T-statistics for each coefficient (numeric vector).
#' -`p_val`: P-values for each coefficient (numeric vector).
#' -`R2`: R-squared value (numeric).
#' -`logLike`: Log-likelihood of the model (numeric).
#' -`DW`: Durbin-Watson statistic (numeric).
#' 
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @examples
#' # Example usage:
#' # Set up a simple linear regression problem
#' set.seed(123)
#' Y <- rnorm(100)
#' X <- cbind(1, rnorm(100))
#' 
#' # Estimate parameters using OLS with White robust standard errors
#' result <- OLS(Y, X, add_constant = TRUE)
#' 
#' @keywords internal
#' 
#' @export
OLS <- function(Y, X, add_constant){
  
  # Purpose: 
  # Ordinary least squares with White robust standard error
  # ---
  # Model:
  # Yi = Xi * Beta + ui , where ui ~ N(0,s^2)
  # ---
  # Algorithm: 
  # inv(X'*X)* (X'*Y)
  # ---
  # Usage:
  # Y = dependent variable (n * 1 vector)
  # X = regressors (n * k matrix)
  # add_constant = whether to add a constant to X (default = 0)
  # ---
  # Returns:
  # estimator = estimator corresponding to the k regressors
  # SE = standard error of the estimator
  # SE_robust = White robust standard error of the estimator
  # sigma2 = estimated variance of disturbances
  # resid = residuals series of the regression
  # In the absence of returning arguments, estimation results will be displayed on screen
  # 
  # Written by Hang Qian, Iowa State University
  # Contact me:  matlabist@gmail.com
  nrow_x <- nrow(X)
  ncol_x <- ncol(X)
  nrow_y <- nrow(Y)
  ncol_y <- ncol(Y)
  
  
  nobs <- nrow(X)
  nvar <- ncol(X)
  if (add_constant == 1){
    X <- cbind(matrix(1,nobs,1),X)  
    nvar <- nvar + 1
  }
  XX <- t(X)%*%X
  inv_XX <- solve(XX)
  estimator <- inv_XX%*%(t(X)%*%Y)
  resid <- (Y - X%*%estimator)
  RSS <- t(resid)%*%resid
  sigma2 <- RSS/(nobs-nvar-1)
  cov_mat <- inv_XX*c(sigma2)
  SE <- sqrt(diag(cov_mat))
  
  X_trans <- X*resid[,rep(1,nvar)]
  cov_mat_robust <- inv_XX%*%(t(X_trans)%*%X_trans)%*%inv_XX
  cov_mat_robust <- nobs/(nobs-nvar-1)*cov_mat_robust
  SE_robust <- sqrt(diag(cov_mat_robust))
  
  t_stat <- estimator/ SE_robust
  
  p_val <- (1 - pt(abs(t_stat),nobs-nvar-1)) * 2
  Y_demean <- Y - mean(Y)
  R2 <- 1-RSS/(t(Y_demean)%*%Y_demean)
  log_like <- -nobs/2 * (1 + log(2*pi) + log(RSS/nobs))
  DW <- sum(diff(resid)^2) / RSS
  
  # --- Output
  output <- list(y = Y, w = X, coef = estimator, SE = SE, 
                 SE_robust = SE_robust, sigma2 = sigma2, resid = resid, 
                 SSR = RSS, t_stat = t_stat, p_val = p_val, R2 = R2, 
                 logLike = log_like, DW = DW)
  return(output)
}


#' @title Diagonal Partition of Observations by Regimes
#' 
#' @description
#' Procedure to construct the diagonal partition of z with m breaks at date b.
#' 
#' @details
#' Note: This code is an adaptation of the one originally written by Yohei Yamamoto
#' and Pierre Perron for MATLAB. Original codes can be found on Pierre Perron's
#' website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @param zz Numeric matrix:  the observations.
#' @param m Integer: indicating the number of breaks.
#' @param bb Numeric matrix: the break dates.
#' 
#' @return
#' A numeric matrix representing the diagonal partition of z with m breaks at date b.
#' 
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @keywords internal
#' 
#' @export
pzbar <- function(zz, m, bb){
  # procedure to construct the diagonal partition of z with m break at date b.  
  nt <- nrow(zz)
  q1 <- ncol(zz)
  zb <- matrix(0, nt, (m+1)*q1)
  zb[(1:bb[1,1]),(1:q1)] <- zz[(1:bb[1,1]),]
  i <- 2
  while (i <= m){
    zb[((bb[i-1,1]+1):bb[i,1]),(((i-1)*q1+1):(i*q1))] <- zz[((bb[i-1,1]+1):bb[i,1]),]
    i <- i + 1
  }
  zb[((bb[m,1]+1):nt),((m*q1+1):((m+1)*q1))] <- zz[((bb[m,1]+1):nt),]
  return(zb)
}


#' @title Determine if Square Matrix is Positive Definite
#' 
#' @description
#' This function tests if a square matrix A is positive definite.
#' A must be a square matrix.
#' 
#' @param A Numeric matrix: matrix to be tested for positive definiteness.
#' 
#' @return
#' A logical value indicating whether the matrix is positive definite.
#' 
#' @references
#' Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020),
#' "Testing Jointly for Structural Changes in the Error Variance and Coefficients
#' of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @keywords internal
#' 
#' @export
isposdef <- function(A){
  # Tests if Matrix A is positive definited.
  # A must be square
  #
  #  Brett Gray 16-Jun-03 (originally written in MATLAB)
  nr = nrow(A)
  nc = ncol(A)
  if (nr != nc){
    stop('ERROR - A must be a square matrix')
  }
  EigenVals <- eigen(A)$values
  Result = TRUE
  n = 1
  while (n <= nr){
    if (EigenVals[n] <= 0){
      Result = FALSE
    }
    n = n + 1
  }
  return(Result)    
}



#' @title Inverse of Positive Definite Matrix 
#' 
#' @description Computes the inverse of a positive definite matrix. If the matrix is not
#' positive definite, it regularizes it and then computes the inverse.
#' 
#' @param x Numeric matrix: A square matrix.
#' 
#' @return A list containing:
#'   - `xinv`: The inverse of the matrix x.
#'   - `flag`: A flag indicating if the original matrix was positive definite 
#'         (flag = 0) or if regularization was applied (flag = 1).
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for 
#' Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" 
#' \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @keywords internal
#' 
#' @export
invpd <- function(x){
  if (isposdef(x)){
    xinv  <- solve(x)
    flag  <- 0
  }else{
    flag      <- 1
    n         <- nrow(x)
    svd_out   <- svd(x)
    xchk      <- svd_out$u%*%diag(svd_out$d)%*%t(svd_out$v)
    dd        <- svd_out$d + 1000*pracma::eps(x = 1.0)
    di        <- rep(1,n)/dd
    xinv      <- svd_out$u%*%diag(di)%*%t(svd_out$v)
  }
  return(list(xinv = xinv, flag = flag))
}



#' @title Number of Cases Given Breaks
#' 
#' @description Determine Number of Cases to Check Given the Number of Breaks in Mean or Variance
#' 
#' @details Replaces `brcvcase()` and `numcase()` in the original MATLAB code.
#' 
#' @param m Integer: Number of breaks in the mean.
#' @param n Integer: Number of breaks in the variance.
#' 
#' @return A list containing:
#' -`num`: The number of cases to consider.
#' -`cvbrind`: A list of matrices representing cases to check given the breaks 
#'         in mean (`m`) and breaks in variance (`n`).
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for 
#' Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" 
#' \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @keywords internal
#' 
#' @export
numcase2 <- function(m,n){
  
  # No breaks in mean or var
  if (m==0){
    tmpls <- list(rbind(matrix(0,1,n),matrix(1,1,n)))
  }else if (n==0){
    tmpls <- list(rbind(matrix(1,1,m),matrix(0,1,m)))
  }else{ # some breaks in either mean or var  
    # Cases with no overlap
    K_min <- max(m,n)
    M_tmp <- t(kronecker(diag(2), rep(1,K_min)))
    if (m!=n){
      if (m>n){
        M <- M_tmp[,-((ncol(M_tmp)-abs(m-n)+1):ncol(M_tmp))]   
      }else if (m<n){
        M <- M_tmp[,-(1:abs(m-n))]     
      }
    }else{
      M <- M_tmp
    }
    pcM <- combinat::permn(ncol(M))   # must find way to reduce here. There are many more combinations than we need. 
    expP <- expand.grid(1:length(pcM), 1:length(pcM))
    M_all <- Map(
      function(a,b) rbind( M[1, pcM[[a]]], M[2, pcM[[a]]] ),
      expP[,1],
      expP[,2]
    )
    tmpls <- unique(M_all)
    
    # Cases with some overlap
    overlap_count <- 1
    while (overlap_count<=min(m,n)){
      M_tmp <- t(kronecker(diag(2), rep(1,K_min)))
      M_tmp <- rbind(M_tmp[1,-((ncol(M_tmp)-overlap_count+1):ncol(M_tmp))],
                     M_tmp[2,-(1:overlap_count)])
      if (m!=n){
        if (m>n){
          M <- M_tmp[,-((ncol(M_tmp)-abs(m-n)+1):ncol(M_tmp))]   
        }else if (m<n){
          M <- M_tmp[,-(1:abs(m-n))]     
        }
      }else{
        M <- M_tmp
      }
      pcM <- combinat::permn(ncol(M))
      expP <- expand.grid(1:length(pcM), 1:length(pcM))
      M_all <- Map(
        function(a,b) rbind( M[1, pcM[[a]]], M[2, pcM[[a]]] ),
        expP[,1],
        expP[,2]
      )
      tmpls <- c(tmpls,unique(M_all))
      overlap_count <- overlap_count + 1
    }
  }
  # number of cases to consider
  num <- length(tmpls)
  return(list(num = num, cvbrind = tmpls))
}


#' @title Compute Quadratic Kernel
#' 
#' @description This function computes the quadratic kernel at a given value `x`.
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @param x Value at which to evaluate the quadratic kernel.
#' 
#' @return The value of the quadratic kernel at the given `x`.
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @keywords internal
#' 
#' @export
kern <- function(x){
  # procedure to evaluate the quadratic kernel at some value x.
  del <- 6*pi*x/5
  ker <- 3*(sin(del)/del-cos(del))/(del*del)
  return(ker)
}


#' @title Bartlett Kernel
#' 
#' @description This function computes the Bartlett kernel at a given value `x`.
#' 
#' @param x Value at which to evaluate the Bartlett kernel.
#' 
#' @return The value of the Bartlett kernel at the given `x`.
#' 
#' @keywords internal
#' 
#' @export
bartlett_kern <- function(x){
  if (abs(x) <= 1){
    k <- 1 - abs(x)
  }else if (abs(x)>1){
    k <- 0
  }
  return(k)
}

#' @title Main Kernel Function
#' 
#' @description This function calls other kernel functions depending on the specified "type" of kernel.
#' 
#' @param x Value at which to evaluate the kernel.
#' @param type String: Type of kernel to use. Options include "quadratic" and "bartlett".
#' 
#' @return The value of the specified kernel at the given `x`.
#' 
#' @keywords internal
#' 
#' @export
kernMain <- function(x, type = "quadratic"){
  if (type == "quadratic"){
    ker <- kern(x)
  }else if (type == "bartlett")
    ker <- bartlett_kern(x)
  return(ker)
}


#' @title Bandwidth Based on AR(1) Approximation
#' 
#' @description This function computes the automatic bandwidth based on AR(1)
#' approximation for each vector of the matrix `vhat`. Each vector is given equal weight (weight = 1).
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @param vhat Numeric matrix: Matrix where each column represents a vector.
#' 
#' @return The computed bandwidth based on AR(1) approximation.
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @keywords internal
#' 
#' @export
bandw <- function(vhat){
  # procedure that compute the automatic bandwidth based on AR(1)
  # approximation for each vector of the matrix vhat. Each are given equal
  # weight 1.
  nt <- nrow(vhat)
  d <- ncol(vhat)
  a2n <- 0
  a2d <- 0
  for (i in 1:d){
    b <- olsqr(vhat[(2:nt),i], vhat[(1:(nt-1)),i])
    sig <- (vhat[(2:nt),i] - b%*%vhat[(1:(nt-1)),i])%*%t(vhat[(2:nt),i]-b%*%vhat[(1:(nt-1)),i])
    sig <- sig/(nt-1)
    a2n <- a2n + 4*b*b*sig*sig/(1-b)^8
    a2d <- a2d + sig*sig/(1-b)^4
  }
  a2 <- a2n/a2d
  st <- 1.3221*(a2*nt)^0.2
  return(st)
}


#' @title Compute Robust Standard Errors for the `correct1()` Function
#' 
#' @description This function computes the long-run covariance matrix of `vmat` and applies a small sample correction.
#' 
#' @param vmat Numeric matrix: Matrix for which the long-run covariance matrix is computed.
#' @param vmata Numeric matrix: Matrix used for bandwidth selection.
#' @param kerntype String: Type of kernel function to be used in bandwidth selection.
#' 
#' @return The computed long-run covariance matrix with small sample correction.
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @keywords internal
#' 
#' @export
jhatpr1 <- function(vmat,vmata,kerntype){
  # procedure to compute the long-run covariance matrix of vmat.  
  nt <- nrow(vmat)
  d <- ncol(vmat)
  jhat <- matrix(0, d, d)
  # calling the automatic bandwidth selection
  st <- bandw(vmata)
  # lag 0 covariance
  jhat <- t(vmat)%*%vmat
  # forward sum
  for (j in 1:(nt-1)){
    jhat <- jhat + c(kernMain(j/st, kerntype))*t(matrix(vmat[((j+1):nt),], nt-j, ncol(vmat)))%*%matrix(vmat[(1:(nt-j)),], nt-j, ncol(vmat))
  }
  # backward sum
  for (j in 1:(nt-1)){
    jhat <- jhat + c(kernMain(j/st, kerntype))*t(matrix(vmat[(1:(nt-j)),], nt-j, ncol(vmat)))%*%matrix(vmat[((j+1):nt),], nt-j, ncol(vmat))
  }
  # small sample correction
  jhat <- jhat/(nt-d)
  return(jhat)
}




#' @title Compute Robust Standard Errors in Joint Test of Structural Changes (Perron et al., 2021)
#' 
#' @description This function computes robust standard errors for the joint test of structural changes.
#' 
#' @param vmat Numeric matrix: Matrix of residuals under the null hypothesis or alternative hypothesis.
#' @param vmata Numeric matrix: Matrix of residuals used for bandwidth selection.
#' @param prewhit Logical: indicating whether to prewhiten the residuals (default is FALSE).
#' @param typekb Integer: indicating the type of residuals to be used: 0 for residuals under H0, 1 for residuals under H1, and 2 for hybrid method.
#' @param kerntype String: Type of kernel function to be used in bandwidth selection.
#' 
#' @return Robust standard errors for the joint test of structural changes.
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @keywords internal
#' 
#' @export
correct1 <- function(vmat, vmata, prewhit, typekb, kerntype){
  # typekb=0 LRV constructed residuals under H0
  # typekb=1 LRV constructed residuals under H1
  # typekb=2 LRV constructed by hybrid method
  
  if (typekb==0){
    vmata    <- vmat
  }else if (typekb==1){
    vmat     <- vmata    
  }
  nt <- nrow(vmat) 
  d <- ncol(vmat)
  bmat <- matrix(0, d, d)
  
  bmata <- matrix(0, d, d)
  vstar <- matrix(0, nt-1, d)
  vstara <- matrix(0, nt-1, d)
  
  if (prewhit  == TRUE){
    for (i in 1:d){
      b           <- olsqr(vmat[2:nt,i], vmat[1:(nt-1),])
      bmat[i,]    <- t(b)
      vstar[,i]   <- vmat[2:nt,i] - vmat[1:(nt-1),]%*%b
      ba          <- olsqr(vmata[2:nt,i], vmata[1:(nt-1),])
      bmata[i,]   <- t(ba)
      vstara[,i]  <- vmata[2:nt,i] - vmata[1:(nt-1),]%*%ba
    }
    jh            <- jhatpr1(vstar, vstara, kerntype)
    hac           <- invpd(diag(d)-bmat)$xinv%*%jh%*%t(inv(diag(d)-bmat))  
  }else{
    hac           <- jhatpr1(vmat, vmata, kerntype)
  }
  return(hac)
}


#' @title Diagonal Matrix with Break Fractions
#' 
#' @description This function constructs a diagonal matrix of dimension m+1 with ith entry (T_i - T_i-1)/T.
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @param b Numeric matrix: Matrix of breakpoints.
#' @param m Integer: Number of breakpoints.
#' @param Tsize Integer: Size of the dataset.
#' 
#' @return Diagonal matrix with break fractions.
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @keywords internal
#' 
#' @export
plambda <- function(b, m, Tsize){
  lambda <- matrix(0, m+1, m+1)
  lambda[1,1] <- b[1,1]/Tsize
  k <- 2  
  while (k <= m){
    lambda[k,k] <-(b[k,1]-b[k-1,1])/Tsize
    k <- k + 1    
  }  
  lambda[m+1,m+1] <- (Tsize-b[m,1])/Tsize
  return(lambda)
}



#' @title Log Likelihood Given Residuals
#' 
#' @param res Numeric vector: Vector of residuals.
#' @param n Integer: Number of structural changes.
#' @param brv Numeric vector: Vector of breakpoints.
#' 
#' @return A list containing the log likelihood and tao values.
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @keywords internal
#' 
#' @export
ploglik <- function(res,n,brv){
  bigt <- nrow(res)
  seg  <- as.matrix(c(0,brv,bigt))
  tao  <- matrix(0,bigt,1)
  loglik <- 0
  for (k in 1:(n+1)){
    i <- seg[k,] + 1
    j <- seg[(k+1),]
    vvar <- (t(res[i:j,])%*%res[i:j,])/(j-i+1)
    tao[i:j,] <- ((res[i:j,]^2)/c(vvar)) - 1
    loglik <- loglik - (1/2)*(j-i+1)*(log(2*pi)+1+log(vvar))
  }
  return(list(loglik = loglik, tao = tao))
}


#' @title Compute global break dates in partial structural break model
#' 
#' @description
#' This function computes global break dates in a partial structural break model.
#' 
#' @param y Numeric vector or matrix: the response variable.
#' @param z Numeric matrix: the matrix of covariates related to the structural breaks.
#' @param x Numeric matrix: the matrix of covariates not related to the structural breaks.
#' @param m Integer: the maximum number of breaks to consider.
#' @param h Integer: bandwidth parameter for the kernel density estimator in dating_purescSSR.
#' @param thtol Numeric: convergence threshold for the iterative algorithm (default is 1e-6).
#' @param maxi Integer: maximum number of iterations for the iterative algorithm (default is 10000).
#' 
#' @details
#' Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @return
#' A list containing:
#' -`glob`: Numeric matrix, the global break dates.
#' -`datevec`: Numeric matrix, the dates corresponding to each break.
#' -`bigvec`: Numeric matrix, the result of the iterative algorithm.
#' -`convergence_msg`: Character, a message indicating whether the algorithm converged or reached the maximum iterations.
#' 
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @examples
#' # Example usage:
#' # dating_partscSSR(y, z, x, m, h)
#' 
#' @seealso
#' \code{\link{dating_purescSSR}}, \code{\link{pzbar}}, \code{\link{olsqr}}
#' 
#' 
#' @export
dating_partscSSR <- function(y, z, x, m, h, thtol = 1e-6, maxi = 10000){
  # ----- Set some values
  bigT      <- nrow(y)
  q         <- ncol(z)
  p         <- ncol(x)
  glb       <- matrix(0, m, 1)
  datevec   <- matrix(0, m, m)
  for (mi in 1:m){
    qq <- p+q
    zz <- cbind(z, x)
    puresc_out <- dating_purescSSR(y, zz, mi, h)
    globnl <- puresc_out$glob
    datenl <- puresc_out$datevec
    bigvec <- puresc_out$bigvec
    
    xbar <- pzbar(x, mi, as.matrix(datenl[(1:mi),mi]))
    zbar <- pzbar(z, mi, as.matrix(datenl[(1:mi),mi]))
    theta <- olsqr(y, cbind(zbar, xbar))
    delta1 <- as.matrix(theta[(1:(q*(mi+1))),1])
    beta1 <- olsqr(y - zbar%*%delta1, x)
    ssr1 <- t(y-zbar%*%delta1-x%*%beta1)%*%(y-zbar%*%delta1-x%*%beta1)
    # Starting the iterations
    length <- 99999999
    i <- 1
    while ((length > thtol) & (i <=maxi)){  # May be useful to writ write this in C++ if it takes many iterations
      nlsc_out  <- dating_purescSSR(y - x%*%beta1, z, mi, h)
      datenl    <- nlsc_out$datevec
      bigvec  <- nlsc_out$bigvec
      
      zbar      <- pzbar(z, mi, as.matrix(datenl[(1:mi),mi]))
      theta1    <- olsqr(y,cbind(zbar,x))
      delta1    <- as.matrix(theta[(1:(q*(mi+1))),1])
      beta1     <- as.matrix(theta1[((q*(mi+1)+1):(q*(mi+1)+p)),1])
      ssrn      <- t(y - cbind(zbar, x)%*%theta1)%*%(y - cbind(zbar, x)%*%theta1)
      
      # Calculate overall SRR and check if significantly smaller.
      length    <- abs(ssrn-ssr1)
      i     <- i + 1
      ssr1  <- ssrn
      glb[mi,1] <- ssrn
      datevec[(1:mi),mi] <- datenl[(1:mi), mi]
    }
    if ((length > thtol) & (i > maxi)){
      msg     <- 'The number of iterations has reached the upper limit'
    }else{
      msg   <- 'converged'
    }
  }
  return(list(glob = glb, datevec = datevec, bigvec = bigvec, convergence_msg = msg))
}

#' @title Get bigvec of residuals
#' 
#' @description
#' This function calculates the bigvec of residuals, replacing the residuals() function in the original MATLAB code.
#' 
#' @param y Numeric vector or matrix: the response variable.
#' @param z Numeric matrix: the matrix of covariates.
#' @param b Numeric matrix: the matrix of coefficients.
#' @param q Integer: the number of coefficients.
#' @param m Integer: the number of structural breaks.
#' 
#' @return
#' A matrix containing the bigvec of residuals.
#' 
#' @references
#' Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @examples
#' # Example usage:
#' # bigvec_residuals(y, z, b, q, m)
#' 
#' @keywords internal
#' 
#' @export
bigvec_residuals <- function(y,z,b,q,m){
  bigt <- nrow(y)
  bigvec <- matrix(0, bigt*(m+1), 1)
  for (i in 1:(m+1)){
    bigvec[((i-1)*bigt+1):(i*bigt),1] <- (y - z%*%b[((i-1)*q+1):(i*q),1])^2
  }
  return(bigvec)
}

#' @title Find segments
#' 
#' @description
#' This function finds segments based on specified coefficient and variance breaks.
#' 
#' @param K Integer: total number of breaks.
#' @param brk Numeric vector: vector of break points.
#' @param m Integer: number of coefficient breaks.
#' @param n Integer: number of variance breaks.
#' @param cbrind Numeric vector: indicator vector for coefficient breaks.
#' @param vbrind Numeric vector: indicator vector for variance breaks.
#' @param q Integer: number of coefficients.
#'
#' @return
#' A list containing:
#' -`brc`: Numeric matrix, matrix of coefficient breaks.
#' -`brv`: Numeric matrix, matrix of variance breaks.
#' -`R`: Numeric matrix, matrix used for testing joint structural changes.
#'
#' @references
#' Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @examples
#' # Example usage:
#' # segmake(K, brk, m, n, cbrind, vbrind, q)
#' 
#' @keywords internal
#' 
#' @export
segmake <- function(K,brk,m,n,cbrind,vbrind,q){
  if (sum(cbrind)!=m){
    print('number of coefficient breaks is misspecified.')   
  }
  if (sum(vbrind)!=n){
    print('number of variance breaks is misspecified.')
  }
  
  brc        <- matrix(0,0,0)
  brv        <- matrix(0,0,0)
  R          <- matrix(0, q*(K+1), q*(m+1))
  R[(1:q),(1:q)] <- diag(q)
  
  ri         <- 0
  for (i in 1:K){
    if (cbrind[i,1]==1){
      brc = c(brc,brk[i]) 
      ri  = ri+1
    }
    R[(q*i+1):(q*(i+1)),(q*ri+1):(q*(ri+1))] <- diag(q)
    if (vbrind[i,1]==1){
      brv = c(brv, brk[i])
    }
  }
  return(list(brc = as.matrix(brc), brv = as.matrix(brv), R = R))
}


#' @title Estimate FGLS
#' 
#' @description
#' Procedure to estimate coefficients and residuals by FGLS given break dates in coefficients and variance.
#' 
#' @param y Numeric vector or matrix: the response variable.
#' @param z Numeric matrix: the matrix of covariates.
#' @param q Integer: number of coefficients.
#' @param x Numeric matrix: the matrix of covariates not related to the structural breaks.
#' @param p Integer: the number of covariates not related to the structural breaks.
#' @param bigt Integer: the total number of observations.
#' @param K Integer: the total number of breaks.
#' @param brk Numeric vector: vector of break points.
#' @param R Numeric matrix: matrix used for testing joint structural changes.
#' @param n Integer: number of variance breaks.
#' @param brv Numeric vector: vector of variance breaks.
#' @param rest Integer: indicator for whether to restrict or not.
#'
#' @return
#' A list containing:
#' -`nbeta`: Numeric matrix, estimated coefficients.
#' -`res`: Numeric matrix, residuals.
#' 
#' @references
#' Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#'
#' @examples
#' # Example usage:
#' # estimbr(y, z, q, x, p, bigt, K, brk, R, n, brv, rest)
#' 
#' 
#' @export
estimbr <- function(y,z,q,x,p,bigt,K,brk,R,n,brv,rest){
  if (p>=1){
    y <- y - x%*%invpd(t(x)%*%x)$xinv%*%t(x)%*%y
    z <- z - x%*%invpd(t(x)%*%x)$xinv%*%t(x)%*%z
  }
  if (K==0){
    zbar <- z  
  }else{
    zbar <- pzbar(z, K, brk)
  }
  if (rest==1){
    zbar <- zbar%*%R
  }
  beta0 <- invpd(t(zbar)%*%zbar)$xinv%*%t(zbar)%*%y
  res0  <- y - zbar%*%beta0
  
  ibigv <- diag(bigt)
  vseg  <- as.matrix(c(0,brv,bigt))
  for (k in 1:(n+1)){
    i <- vseg[k,]+1
    j <- vseg[(k+1),]
    vvar  <- (t(res0[i:j,])%*%res0[i:j,])/(j-i+1)
    ibigv[i:j,i:j] <- diag(j-i+1)*c(invpd(vvar)$xinv)
  }
  beta <- invpd(t(zbar)%*%ibigv%*%zbar)$xinv%*%t(zbar)%*%ibigv%*%y
  
  bstar      <- beta+10
  itr        <- 1
  while (max(abs(bstar-beta))>=1e-6 & itr<=100){
    bstar    <- beta
    res      <- y - zbar%*%beta
    
    ibigv    <- diag(bigt)
    for (k in 1:(n+1)){
      i     <- vseg[k,]+1
      j     <- vseg[(k+1),]
      vvar  <- (t(res[i:j,1])%*%res[i:j,1])/(j-i+1)
      ibigv[i:j,i:j] <- diag(j-i+1)*c(invpd(vvar)$xinv)
    }
    beta    <- invpd(t(zbar)%*%ibigv%*%zbar)$xinv%*%t(zbar)%*%ibigv%*%y
    itr     <- itr+1
    if (itr==100){
      print('The iteration has reached the upper limit')
    }
  }
  if (rest==1){
    nbeta   <- R%*%beta  
  }else{
    nbeta   <- beta
  }
  res       <- y - zbar%*%beta
  return(list(nbeta = nbeta, res = res))
}


#' @title Joint estimation of mean and variance coefficients
#' 
#' @description
#' Procedure to jointly estimate the coefficient and variance break dates following section 5.1 of Qu and Perron (2007). The code follows est.m of Qu and Perron (2007).
#' 
#' @param y Numeric vector or matrix: the response variable.
#' @param z Numeric matrix: the matrix of covariates.
#' @param q Integer: number of coefficients.
#' @param x Numeric matrix: the matrix of covariates not related to the structural breaks.
#' @param p Integer: the number of covariates not related to the structural breaks.
#' @param K Integer: the total number of breaks.
#' @param bigt Integer: the total number of observations.
#' @param h Integer: bandwidth parameter for the kernel density estimator in dating_loglik.
#' @param m Integer: the number of coefficient breaks to consider.
#' @param n Integer: the number of variance breaks to consider.
#' @param cbrind Numeric vector: indicator vector for coefficient breaks.
#' @param vbrind Numeric vector: indicator vector for variance breaks.
#'
#' @return
#' A list containing:
#' -`brk`: Numeric matrix, estimated break dates.
#' -`beta`: Numeric matrix, estimated coefficients.
#' -`brc`: Numeric matrix, matrix of coefficient breaks.
#' -`brv`: Numeric matrix, matrix of variance breaks.
#' -`res`: Numeric matrix, residuals.
#'
#' @references
#' Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#'
#' @examples
#' # Example usage:
#' # estdate(y, z, q, x, p, K, bigt, h, m, n, cbrind, vbrind)
#' 
#' @export
estdate <- function(y,z,q,x,p,K,bigt,h,m,n,cbrind,vbrind){
  # step 1
  datevec <- dating_MLE(y,z,q,x,p,h,K,bigt)$datevec
  brk <- as.matrix(datevec[,K])
  
  # step 2
  segmake_out <- segmake(K,brk,m,n,cbrind,vbrind,q)
  estimbr_out <- estimbr(y,z,q,x,p,bigt,K,brk,segmake_out$R,n,segmake_out$brv,1)
  if (p>=1){
    y <- y - x%*%invpd(t(x)%*%x)$xinv%*%t(x)%*%y
    z <- z - x%*%invpd(t(x)%*%x)$xinv%*%t(x)%*%z  
  }
  
  # step 3-5
  diff <- 0
  maxiter <- 10
  while ((diff!=-1) & (diff<maxiter)){
    bigvec      <- bigvec_residuals(y,z,estimbr_out$nbeta,q,K)
    datevec     <- dating_loglik(bigvec,h,K,bigt)$datevec
    if (all(as.matrix(datevec[,K])==brk)){
      diff      <- -1  
    }else{
      brk         <- as.matrix(datevec[,K])
      segmake_out <- segmake(K,brk,m,n,cbrind,vbrind,q)
      estimbr_out <- estimbr(y,z,q,x,p,bigt,K,brk,segmake_out$R,n,segmake_out$brv,1)
      diff        <- diff + 1
    }
    if (diff==maxiter){
      print('cannot find the converging break dates') 
    }
  }
  return(list(brk = brk, beta = estimbr_out$nbeta, brc = segmake_out$brc, brv = segmake_out$brv, res = estimbr_out$res))
}


#' @title Compute long run variance
#' 
#' @description 
#' This function computes the long-run covariance matrix of a given matrix.
#' 
#' @details 
#' Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @param vmat Numeric matrix: the input matrix for which the long-run covariance matrix is to be computed.
#' 
#' @return
#' A numeric matrix representing the long-run covariance matrix of the input matrix.
#'
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#'
#' @keywords internal
#' 
#' @export
jhatpr <- function(vmat){
  # procedure to compute the long-run covariance matrix of vmat.  
  nt <- nrow(vmat)
  d <- ncol(vmat)
  jhat <- matrix(0, d, d)
  # calling the automatic bandwidth selection
  st <- bandw(vmat)
  # lag 0 covariance
  jhat <- t(vmat)%*%vmat
  # forward sum
  for (j in 1:(nt-1)){
    jhat <- jhat + c(kern(j/st))*t(matrix(vmat[((j+1):nt),], nt-j, ncol(vmat)))%*%matrix(vmat[(1:(nt-j)),], nt-j, ncol(vmat))
  }
  # backward sum
  for (j in 1:(nt-1)){
    jhat <- jhat + c(kern(j/st))*t(matrix(vmat[(1:(nt-j)),], nt-j, ncol(vmat)))%*%matrix(vmat[((j+1):nt),], nt-j, ncol(vmat))
  }
  # small sample correction
  jhat <- jhat/(nt-d)
  
  return(jhat)
}


#' @title Compute robust standard errors
#' 
#' @description
#' This function computes robust standard errors using a procedure that activates the computation of robust standard errors.
#' 
#' @details
#' Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @param reg Numeric matrix: the regression matrix.
#' @param res Numeric vector: the residuals vector.
#' @param prewhit Integer: indicator for prewhitening. If set to 1, prewhitening is applied; otherwise, it is skipped.
#'
#' @return
#' A numeric matrix representing the robust standard errors.
#'
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#'
#' @keywords internal
#' 
#' @export
correct <- function(reg, res, prewhit){
  # main procedures which activates the computation of robust standard errors  
  nt    <- nrow(reg)
  d     <- ncol(reg)
  b     <- matrix(0, d, 1)
  bmat  <- matrix(0, d, d)
  vstar <- matrix(0, nt-1, d)
  vmat  <- matrix(0, nt, d)
  
  # First construct the matrix z_t*u_t.  
  for (i in 1:d){
    vmat[,i] <- reg[,i]*res
  }
  # Procedure that applies prewhitenning to the matrix vmat by filtering with
  # a VAR(1). If prewhit=0, it is skipped.
  if (prewhit==1){
    for (i in 1:d){
      b = olsqr(vmat[(2:nt), i], vmat[(1:(nt-1)),])
      bmat[i,] <- b
      vstar[,i] <- vmat[(2:nt),i] - vmat[(1:(nt-1)),]%*%b  
    }
    # Call the kernel on the residuals
    jh <- jhatpr(vstar)
    # recolor
    hac <- solve(diag(d)-bmat)%*%jh%*%t(solve(diag(d)-bmat))
  }else{
    hac <- jhatpr(vmat)
  }
  return(hac)
}

#' @title Diagonal matrix with variance for regime i
#' 
#' @description
#' This function computes a diagonal matrix of dimension i+1 with the ith entry
#' being the estimate of the variance of the residuals for segment i.
#' 
#' @details
#' Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @param res Numeric vector: the residuals vector.
#' @param b Numeric matrix: the break dates matrix.
#' @param q Integer: the number of regressors.
#' @param i Integer: the segment index.
#' @param nt Integer: the total number of observations.
#' 
#' @return
#' A diagonal matrix of dimension i+1 with the ith entry being the estimate of the variance of the residuals for segment i.
#'
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @keywords internal
#' 
#' @export
psigmq <- function(res, b, q, i, nt){
  # procedure that computes a diagonal matrix of dimension i+1 with ith entry
  # the estimate of the variance of the residuals for segment i.  
  sigmat <- matrix(0, i+1,i+1)
  sigmat[1,1] <- t(res[(1:b[1,1]),1])%*%res[(1:b[1,1]),1]/b[1,1]
  kk <- 2
  while (kk <= i){
    sigmat[kk,kk] <- t(res[((b[(kk-1),1]+1):b[kk,1]),1])%*%res[((b[(kk-1),1]+1):b[kk,1]),1]/(b[kk,1]-b[kk-1,1])
    kk <- kk + 1
  }
  sigmat[i+1,i+1] <- t(res[((b[i,1]+1):nt),1])%*%res[((b[i,1]+1):nt),1]/(nt-b[i,1])
  return(sigmat)
}


#' @title Covariance matrix of estimates delta.
#' 
#' @description
#' This function computes the covariance matrix of the estimates delta.
#' 
#' @details
#' Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#' 
#' @param y Numeric vector: the dependent variable.
#' @param z Numeric matrix: the matrix of regressors.
#' @param i Integer: the number of structural breaks.
#' @param q Integer: the number of regressors.
#' @param Tsize Integer: the total number of observations.
#' @param b Numeric matrix: the break dates matrix.
#' @param prewhit Integer: set to 1 to apply AR(1) prewhitening prior to estimating 
#'                       the long-run covariance matrix.
#' @param robust Integer: set to 1 to allow for heterogeneity and autocorrelation 
#'                       in the residuals, 0 otherwise.
#' @param x Numeric matrix: additional exogenous regressors.
#' @param p Integer: the number of additional regressors.
#' @param withb Integer: estimate covar matrix of param estimates including the constant betas.
#' @param hetdat Integer: option for the construction of the F tests.
#'                       Set to 1 to allow different moment matrices of the regressors across segments.
#'                       If hetdat=0, the same moment matrices are assumed for each segment and estimated from the full sample.
#'                       It is recommended to set hetdat=1 if p>0.
#' @param hetvar Integer: option for heteroscedasticity in the covariance matrix.
#' 
#' @return
#' The covariance matrix of the estimates delta.
#'
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @keywords internal
#' 
#' @export
pvdel <- function(y, z, i, q, Tsize, b, prewhit, robust, 
                  x, p, withb, hetdat, hetvar){
  # procedure that compute the covariance matrix of the estimates delta.
  
  # prewhit: set to 1 if want to apply AR(1) prewhitening prior to estimating 
  #          the long run covariance matrix.
  
  # robust:  set to 1 if want to allow for heterogeneity and autocorrelation 
  #          in the residuals, 0 otherwise. The method used is Andrews(1991) 
  #          automatic bandwidth with AR(1) approximation and the quadratic 
  #          kernel. Note: Do not set to 1 if lagged dependent variables are 
  #          included as regressors.
  
  # hetdat: option for the construction of the F tests. Set to 1 if want to
  #         allow different moment matrices of the regressors across segments. 
  #         If hetdat=0, the same moment matrices are assumed for each segment 
  #         and estimated from the full sample. It is recommended to set 
  #         hetdat=1 if p>0.
  
  # withb: Estimate covar matrix of param estimates including the constant betas 
  #         (only used in model estimation after break dates are chosen)
  
  ev <- matrix(1, i+1,1)
  zbar <- pzbar(z, i, b)
  
  if (p==0){
    delv <- olsqr(y, zbar)
    res <- y - zbar%*%delv
    reg <- zbar
  }else{
    delv <- olsqr(y, cbind(zbar, x))
    res <- y - cbind(zbar, x)%*%delv
    if (withb==0){
      reg <- zbar - x%*%solve(t(x)%*%x)%*%t(x)%*%zbar
    }else{
      reg <- cbind(zbar, x)
    }
  }
  
  if (robust==0){
    # section on testing with no serial correlation in errors
    if (p==0){
      if (hetdat==1 & hetvar==0){
        sig <- c(t(res)%*%res/Tsize)
        vdel <- sig*solve(t(reg)%*%reg)
      }else if (hetdat==1 && hetvar==1){
        sig <- psigmq(res, b, q, i, Tsize)
        vdel <- kronecker(sig,diag(q))%*%solve(t(reg)%*%reg)
      }else if (hetdat==0 && hetvar==0){
        lambda <- plambda(b, i, Tsize)
        sig <- c(t(res)%*%res/Tsize)
        vdel <- sig*solve(kronecker(lambda,(t(z)%*%z)))
      }else if (hetdat==0 && hetvar==1){
        lambda <- plambda(b, i, Tsize)
        sig <- psigmq(res, b, q, i, Tsize)
        vdel <- kronecker(sig,diag(q))%*%solve(kronecker(lambda,(t(z)%*%z)))
      }
    }else{
      if (hetdat==0){
        warning("The case hetdat=0 is not allowed. 'vdel' is returned as zeros.")
        vdel <- matrix(0,q*(i+1),q*(i+1))
      }
      if (hetdat==1 & hetvar==0){
        sig <- c(t(res)%*%res/Tsize)
        vdel <- sig*solve(t(reg)%*%reg)
      }else if (hetdat==1 & hetvar==1){
        wbar <- pzbar(reg, i, b)
        ww <- t(wbar)%*%wbar
        sig <-psigmq(res, b, q, i, Tsize)
        gg <- matrix(0,(i+1)*q+p*withb, (i+1)*q+p*withb)
        ie <-1
        while (ie<=(i+1)){
          gg <- gg + sig[ie,ie]*ww[((ie-1)*((i+1)*q+p*withb)+1:ie*((i+1)*q+p*withb)),((ie-1)*((i+1)*q+p*withb)+1:ie*((i+1)*q+p*withb))]
          ie <- ie + 1
        }
        vdel <- solve(t(reg)%*%reg)%*%gg%*%solve(t(reg)%*%reg)
      }
    }
  }else{
    if (hetdat==0){
      warning("The case hetdat=0 is not allowed. 'vdel' is returned as zeros.")
      vdel <- matrix(0,q*(i+1),q*(i+1))
    }
    if (p==0){
      if (hetvar==1){
        hac <- matrix(0, (i+1)*q, (i+1)*q)
        vdel <- matrix(0, (i+1)*q, (i+1)*q)
        hac[(1:q),(1:q)] <- b[1,1]*correct(as.matrix(z[(1:b[1,1]),]), res[(1:b[1,1]),1], prewhit)
        if (i>=2){
          for (j in 2:i){
            hac[(((j-1)*q+1):(j*q)),(((j-1)*q+1):(j*q))] <- (b[j,1]-b[(j-1),1])*correct(as.matrix(z[((b[(j-1),1]+1):b[j,1]),]), res[((b[(j-1),1]+1):b[j,1]),1], prewhit)
          }  
        }
        hac[((i*q+1):((i+1)*q)),((i*q+1):((i+1)*q))] <- (Tsize-b[i,1])*correct(as.matrix(z[((b[i,1]+1):Tsize),]),
                                                                               res[((b[i,1]+1):Tsize),1], prewhit)
        vdel <- solve(t(reg)%*%reg)%*%hac%*%solve(t(reg)%*%reg)
      }else{
        hac <- correct(as.matrix(z), res, prewhit)
        lambda <-  plambda(b, i, Tsize)
        vdel <- Tsize*solve(t(reg)%*%reg)%*%kronecker(lambda,hac)%*%solve(t(reg)%*%reg)
      }
    }else{
      hac <- correct(as.matrix(reg), res, prewhit)
      vdel <- Tsize*solve(t(reg)%*%reg)%*%hac%*%solve(t(reg)%*%reg)
    }
  } 
  return(vdel)
}



#' @title P-value for critical values of break dates
#'
#' @details This function computes the p-value that is used to determine the critical values and hence the confidence interval for break dates.
#'
#' @param x Numeric: value at which to evaluate the function.
#' @param bet Numeric: Parameter for the function.
#' @param alph Numeric: Parameter for the function.
#' @param b Numeric: Parameter for the function.
#' @param deld Numeric: Parameter for the function.
#' @param gam Numeric: Parameter for the function.
#'
#' @return The computed p-value
#'
#' @keywords internal
#' 
#' @export
funcg <- function(x,bet,alph,b,deld,gam){
  if (x <= 0){
    xb <- bet*sqrt(abs(x))
    if (abs(xb) <= 30){
      g <- -sqrt(-x/(2*pi))*exp(x/8)-(bet/alph)*exp(-alph*x)*pnorm(-bet*sqrt(abs(x)))+((2*bet*bet/alph)-2-x/2)*pnorm(-sqrt(abs(x))/2)         
    }else{
      aa <- log(bet/alph)-alph*x-xb^2/2-log(sqrt(2*pi))-log(xb)
      g <- -sqrt(-x/(2*pi))*exp(x/8)-exp(aa)*pnorm(-sqrt(abs(x))/2)+((2*bet*bet/alph)-2-x/2)*pnorm(-sqrt(abs(x))/2)
    }
  }else{
    xb <- deld*sqrt(x)
    if (abs(xb) <= 30){
      g <- 1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)+(b*deld/gam)*exp(gam*x)*pnorm(-deld*sqrt(x))+(2-b*b*x/2-2*deld*deld/gam)*pnorm(-b*sqrt(x)/2)
    }else{
      aa <- log((b*deld/gam))+gam*x-xb^2/2-log(sqrt(2*pi))-log(xb)
      g <- 1+(b/sqrt(2*pi))*sqrt(x)*exp(-b*b*x/8)+exp(aa)+(2-b*b*x/2-2*deld*deld/gam)*pnorm(-b*sqrt(x)/2) 
    }
  }
  return(g)
}

#' @title Compute critical values for confidence intervals
#'
#' @description
#' The function `cvg` is used to determine critical values for the break dates, which are essential for
#' constructing confidence intervals around the estimated break dates.
#' 
#' @details
#' The critical values are computed using an iterative procedure to find the values at specific
#' significance levels (0.025, 0.05, 0.95, 0.975) for constructing confidence intervals around break dates.
#' The iterative process involves the function `funcg`.
#'
#' @param eta A numeric value representing the ratio of variances.
#' @param phi1s A numeric value representing a parameter.
#' @param phi2s A numeric value representing a parameter.
#'
#' @return A matrix containing critical values for confidence intervals at specified significance levels.
#'
#' @keywords internal
#' 
#' @export
cvg <- function(eta,phi1s,phi2s){
  cvec <- matrix(0, 4, 1)
  a <- phi1s/phi2s
  gam <- ((phi2s/phi1s)+1)*eta/2
  b <- sqrt(phi1s*eta/phi2s)
  deld <- sqrt(phi2s*eta/phi1s)+b/2
  alph <- a*(1+a)/2
  bet <- (1+2*a)/2
  sig <- c(0.025, 0.05, 0.95, 0.975)
  isig <- 1
  while (isig <=4){
    upb <- 2000
    lwb <- -2000
    crit <- 999999
    cct <- 1
    while (abs(crit) >= .000001){
      cct <- cct + 1
      if (cct > 1000){
        stop('the procedure to get critical values for the break dates has 
             reached the upper bound on the number of iterations. This occured
             in the procedure cvg. The resulting confidence interval for this
             break date is incorect')
      }else{
        xx <- lwb+(upb-lwb)/2
        pval <- funcg(xx,bet,alph,b,deld,gam)
        crit <- pval-sig[isig]
        if (crit <= 0){
          lwb <- xx  
        }else{
          upb <- xx
        }
      }
    }
    cvec[isig,1] <- xx
    isig <- isig + 1
  }
  return(cvec)
}


#' @title Compute confidence intervals for breaks
#'
#' @description
#' The function `interval` computes confidence intervals for the break dates based on the "shrinking shifts"
#' asymptotic framework. It uses an iterative procedure to find critical values at specific significance levels
#' for constructing confidence intervals.
#'
#' @details
#' The computation involves calculating critical values using the `cvg` function. The confidence intervals are then
#' determined based on these critical values, and adjustments are made for rounding.
#' 
#' @param y A matrix or vector representing the dependent variable.
#' @param z A matrix representing the regressor matrix.
#' @param zbar A matrix representing the pre-whitened regressor matrix.
#' @param b A matrix representing the break dates.
#' @param q An integer representing the number of regressors.
#' @param m An integer representing the number of breaks.
#' @param robust An indicator (0 or 1) for robust standard errors.
#' @param prewhit An indicator (0 or 1) for AR(1) pre-whitening.
#' @param hetomega An indicator (0 or 1) for heteroscedasticity in the omega matrix.
#' @param hetq An indicator (0 or 1) for heteroscedasticity in the Q matrix.
#' @param x A matrix representing additional regressors.
#' @param p An integer representing the number of additional regressors.
#'
#' @return
#' A matrix containing confidence intervals for each break date. Each row corresponds to a break date, and
#' the columns represent the lower and upper bounds of the confidence interval.
#'
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#'
#' @keywords internal
#' 
#' @export
interval <- function(y,z,zbar,b,q,m,robust,prewhit,hetomega,hetq,x,p){
  cvf <- matrix(0, 4, 1)
  nt <- nrow(y)
  
  if (p==0){
    delta <- olsqr(y, zbar)
    res <- y - zbar%*%delta
  }else{
    dbdel <- olsqr(y, cbind(zbar, x))
    res <- y-cbind(zbar,x)%*%dbdel
    delta <- as.matrix(dbdel[1:((m+1)*q),1])
  }
  bound <- matrix(0, m, 4)
  bf <- matrix(0, m+2, 1)
  bf[2:(m+1), 1] <- b[1:m,1]
  bf[m+2, 1] <- nt
  
  for (i in 1:m){
    delv <- as.matrix(delta[(i*q+1):((i+1)*q), 1] - delta[((i-1)*q+1):(i*q), 1])
    
    if (robust==0){
      if (hetq==1){
        qmat <- t(z[(bf[i,1]+1):bf[i+1,1],])%*%z[(bf[i,1]+1):bf[i+1,1],]/(bf[i+1,1]-bf[i,1])
        qmat1 <- t(z[(bf[i+1,1]+1):bf[i+2,1],])%*%z[(bf[i+1,1]+1):bf[i+2,1],]/(bf[i+2,1]-bf[i+1,1])
      }else{
        qmat <- t(z)%*%z/nt
        qmat1 <- qmat
      }
      
      if (hetomega==1){
        phi1s <- t(res[(bf[i,1]+1):bf[i+1, 1], 1])%*%res[(bf[i,1]+1):bf[i+1,1], 1]/(bf[i+1,1]-bf[i,1])
        phi2s <- t(res[(bf[i+1,1]+1):bf[i+2,1],1])%*%res[(bf[i+1,1]+1):bf[i+2,1], 1]/(bf[i+2,1]-bf[i+1, 1])
      }else{
        phi1s <- t(res)%*%res/nt
        phi2s <- phi1s
      }
      
      eta <- t(delv)%*%qmat1%*%delv/(t(delv)%*%qmat%*%delv)
      cvf <- cvg(eta,phi1s,phi2s)
      
      
      a <- (t(delv)%*%qmat%*%delv)/phi1s
      
      bound[i,1] <- b[i,1] - cvf[4,1]/a
      bound[i,2] <- b[i,1] - cvf[1,1]/a
      bound[i,3] <- b[i,1] - cvf[3,1]/a
      bound[i,4] <- b[i,1] - cvf[2,1]/a
      
    }else{
      if (hetq==1){
        qmat <- t(z[(bf[i,1]+1):bf[i+1,1],])%*%z[(bf[i,1]+1):bf[i+1,1],]/(bf[i+1,1]-bf[i,1])
        qmat1<- t(z[(bf[i+1,1]+1):bf[i+2,1],])%*%z[(bf[i+1,1]+1):bf[i+2,1],]/(bf[i+2,1]-bf[i+1,1])
      }else{
        qmat <- t(z)%*%z/nt
        qmat1 <- qmat
      }
      
      if (hetomega==1){
        omega <- correct(as.matrix(z[(bf[i,1]+1):bf[i+1,1],]), as.matrix(res[(bf[i,1]+1):bf[i+1,1],1]), prewhit)
        omega1 <- correct(as.matrix(z[(bf[i+1,1]+1):bf[i+2,1],]), as.matrix(res[(bf[i+1,1]+1):bf[i+2,1],1]), prewhit)
      }else{
        omega <- correct(z, res, prewhit)
        omega1 <- omega
      }
      phi1s <- t(delv)%*%omega%*%delv/(t(delv)%*%qmat%*%delv)
      phi2s <- t(delv)%*%omega1%*%delv/(t(delv)%*%qmat%*%delv)
      
      eta <- t(delv)%*%qmat1%*%delv/(t(delv)%*%qmat%*%delv)
      
      cvf <- cvg(eta, phi1s, phi2s)
      
      a <- (t(delv)%*%qmat%*%delv)^2/(t(delv)%*%omega%*%delv)
      
      bound[i,1] <- b[i,1] - cvf[4,1]/a
      bound[i,2] <- b[i,1] - cvf[1,1]/a
      bound[i,3] <- b[i,1] - cvf[3,1]/a
      bound[i,4] <- b[i,1] - cvf[2,1]/a
    }
    bound[,1] <- round(bound[,1])
    bound[,2] <- round(bound[,2])+1
    bound[,3] <- round(bound[,3])
    bound[,4] <- round(bound[,4])+1
    
  }
  return(bound)
}


#' @title Estimate model with breaks
#'
#' @description
#' This function estimates the model after breaks are found. It also computes and reports confidence intervals for
#' the break dates. The method used depends on the specification for robust.
#'
#' @details
#' The function uses the `OLS` and `pvdel` functions to estimate the model and calculate corrected standard errors,
#' respectively. Confidence intervals are determined using the `interval` function. Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: [Pierre Perron Codes](https://blogs.bu.edu/perron/codes/)
#'
#' @param y A matrix or vector representing the dependent variable.
#' @param z A matrix representing the regressor matrix.
#' @param x A matrix representing additional regressors.
#' @param m An integer representing the number of breaks.
#' @param b A matrix representing the break dates.
#' @param robust An indicator (0 or 1) for robust standard errors.
#' @param prewhit An indicator (0 or 1) for AR(1) pre-whitening.
#' @param hetomega An indicator (0 or 1) for heteroscedasticity in the omega matrix.
#' @param hetq An indicator (0 or 1) for heteroscedasticity in the Q matrix.
#' @param hetdat An option for the construction of the F tests (0 or 1).
#' @param hetvar An indicator (0 or 1) for heteroscedasticity in the residuals.
#'
#' @return
#' A list containing the estimated model output, corrected standard errors, and confidence intervals for each break date.
#'
#' @references
#' Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#'
#' @export
estim <- function(y,z,x,m,b,robust,prewhit,hetomega,hetq,hetdat,hetvar){
  bigt <- nrow(y)
  q <- ncol(z)
  p <- ncol(x)
  d <- (m+1)*q + p
  
  # Construct the Z_bar matrix. The diagonal partition of Z at the
  # estimated break dates.
  zbar <- pzbar(z, m, b)
  
  if (p==0){
    reg <- zbar 
  }else{
    reg <- cbind(zbar,x)
  }
  
  mdl_out <- OLS(y, reg, FALSE)
  
  # corrected standard errors
  vdel <- pvdel(y, z, m, q, bigt, b, prewhit, robust, x, p, 1, hetdat, hetvar)
  SE_correct <- sqrt(diag(vdel))
  
  # confidence intervals
  bound <- interval(y,z,zbar,b,q,m,robust,prewhit,hetomega,hetq,x,p)
  colnames(bound) <- c("95% low","90% low","90% upp", "95% upp")
  # output
  mdl_out$SE_correct <- SE_correct
  mdl_out$bound <- bound
  return(mdl_out)
}




#' @title Get critical values using equations from Bai Perron (2003)
#'
#' @description
#' This function calculates critical values for multiple structural change tests based on the equations presented by Bai and Perron (2003).
#'
#' @param alpha The significance level (0.1, 0.05, 0.025, or 0.01).
#' @param q The number of regressors.
#' @param k The number of breaks.
#' @param eps The number of parameters for the mean and trend.
#'
#' @return
#' A list containing critical values for the left-tail (cvL), right-tail (cvUD), and two-tail (cvWD) tests. Additionally, a sequence of critical values (cvSeq) is provided.
#'
#' @references
#' Bai, J. and Perron, P. (2003), Critical values for multiple structural change tests. \emph{The Econometrics Journal}, 6: 72-78. https://doi.org/10.1111/1368-423X.00102
#'
#' @keywords internal
#' 
#' @export
getcv <- function(alpha, q, k, eps){
  if (alpha == 0.1){
    cL  <- (7.551+1.718*q-0.041*q^2-0.610*k-15.846*eps+0.025*q/eps)*exp(0.338/k-0.014/(eps*k))/q
    cUD <- (6.917+2.930*q-9.275*eps)*exp(-0.028*q-0.406*eps)/q
    cWD <- (7.316+3.128*q-8.624*eps)*exp(-0.029*q-0.412*eps)/q  
    cS  <- (8.397+3.702*q-0.209*q^2+0.317*(k+1)-3.736*(1/(k+1))-11.596*eps)*exp(-0.027*q+0.005*q^2)
  }else if (alpha == 0.05){
    cL  <- (8.238+1.756*q-0.043*q^2-0.659*k-15.436*eps+0.025*q/eps)*exp(0.389/k-0.013/(eps*k))/q
    cUD <- (8.228+3.095*q-9.644*eps)*exp(-0.029*q-0.291*eps)/q
    cWD <- (9.039+3.318*q-9.969*eps)*exp(-0.030*q-0.327*eps)/q  
    cS  <- (9.879+4.086*q-0.232*q^2+0.322*(k+1)-3.687*(1/(k+1))-11.931*eps)*exp(-0.039*q+0.006*q^2)
  }else if (alpha == 0.025){
    cL  <- (8.968+1.788*q-0.045*q^2-0.715*k-15.255*eps+0.025*q/eps)*exp(0.426/k-0.013/(eps*k))/q
    cUD <- (9.436+3.304*q-9.301*eps)*exp(-0.030*q-0.259*eps)/q
    cWD <- (10.703+3.465*q-11.119*eps)*exp(-0.031*q-0.250*eps)/q
    cS  <- (11.424+4.435*q-0.261*q^2+0.300*(k+1)-3.700*(1/(k+1))-12.292*eps)*exp(-0.047*q+0.006*q^2)
  }else if (alpha == 0.01){
    cL  <- (9.879+1.771*q-0.042*q^2-0.777*k-14.551*eps+0.025*q/eps)*exp(0.471/k-0.012/(eps*k))/q
    cUD <- (11.211+3.366*q-7.279*eps)*exp(-0.027*q-0.268*eps)/q
    cWD <- (13.189+3.346*q-11.870*eps)*exp(-0.026*q-0.173*eps)/q
    cS  <- (13.073+4.954*q-0.317*q^2+0.285*(k+1)-3.419*(1/(k+1))-12.452*eps)*exp(-0.058*q+0.008*q^2)
  }
  output <- list(cvL = cL, cvUD = cUD, cvWD = cWD, cvSeq = cS)
  return(output)
}


#' @title Get critical values for LRT_4 
#' 
#' @description
#' This function retrieves critical values for the LRT_4 test from the paper by Perron, Yamamoto, and Zhou (2020).
#'
#' @param alpha Numeric: The significance level (e.g., 0.05 for a 5% significance level).
#' @param q Integer: The number of coefficients for the regression model.
#' @param trm Numeric: The value of the trimming parameter
#'
#' @return
#' A matrix of critical values for the LRT_4 test. Rows correspond to the number of coefficient breaks, and columns correspond to the number of variance breaks.
#'
#' @references
#' Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#'
#' @keywords internal
#' 
#' @export
getcv4 <- function(alpha, q, trm){
  
  # input of the critical value of the supLR4 test
  # A matrix is provided for each set of (q,signif). Its row 
  # corresponds m (# of coefficient breaks) and its column
  # corresponds to n (# of variance breaks).
  
  siglev <- as.matrix(c(0.10,0.05,0.025,0.01))
  signif <- which(siglev==alpha)
  
  if (trm == .05){
    
    cv <- matrix(0,9,9)
    
    if (q==1){
      if (signif == 1){
        cv<- rbind(
          c(7.08,  7.11, 6.58, 6.27, 5.88, 5.60, 5.28, 5.00, 4.79),
          c(7.13,  7.11, 6.76, 6.38, 6.03, 5.73, 5.47, 5.20, NaN),
          c(6.64,  6.70, 6.43, 6.16, 5.87, 5.63, 5.42, NaN,  NaN),
          c(6.30,  6.38, 6.17, 5.97, 5.76, 5.51, NaN,  NaN,  NaN),
          c(5.92,  6.05, 5.89, 5.77, 5.57, NaN,  NaN,  NaN,  NaN),
          c(5.61,  5.78, 5.63, 5.53, NaN,  NaN,  NaN,  NaN,  NaN),
          c(5.30,  5.48, 5.39, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(4.98,  5.22, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(4.77,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(8.03,  7.82, 7.26, 6.80, 6.37, 6.04, 5.73, 5.39, 5.14),
          c(7.83,  7.71, 7.27, 6.85, 6.46, 6.11, 5.84, 5.56, NaN),
          c(7.24,  7.22, 6.88, 6.58, 6.30, 6.02, 5.75, NaN,  NaN),
          c(6.83,  6.84, 6.63, 6.35, 6.10, 5.90, NaN,  NaN,  NaN),
          c(6.40,  6.49, 6.31, 6.14, 5.93, NaN,  NaN,  NaN,  NaN),
          c(6.02,  6.19, 5.99, 5.92, NaN,  NaN,  NaN,  NaN,  NaN),
          c(5.74,  5.83, 5.74, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(5.38,  5.58, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(5.11,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(8.88,  8.42, 7.87, 7.26, 6.81, 6.41, 6.06, 5.73, 5.44),
          c(8.54,  8.25, 7.76, 7.30, 6.84, 6.46, 6.16, 5.84, NaN),
          c(7.82,  7.68, 7.33, 7.00, 6.63, 6.33, 6.06, NaN,  NaN),
          c(7.29,  7.29, 7.05, 6.71, 6.44, 6.22, NaN,  NaN,  NaN),
          c(6.87,  6.88, 6.64, 6.49, 6.27, NaN,  NaN,  NaN,  NaN),
          c(6.40,  6.56, 6.35, 6.28, NaN,  NaN,  NaN,  NaN,  NaN),
          c(6.10,  6.13, 6.04, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(5.71,  5.89, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(5.42,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(10.08,  9.29, 8.56, 7.91, 7.37, 6.83, 6.56, 6.15, 5.83),
          c(9.39,  8.94, 8.45, 7.88, 7.33, 6.87, 6.58, 6.27, NaN),
          c(8.39,  8.29, 7.88, 7.53, 6.98, 6.72, 6.47, NaN,  NaN),
          c(7.84,  7.84, 7.54, 7.19, 6.84, 6.60, NaN,  NaN,  NaN),
          c(7.47,  7.38, 7.08, 6.85, 6.64, NaN,  NaN,  NaN,  NaN),
          c(6.85,  6.98, 6.76, 6.67, NaN,  NaN,  NaN,  NaN,  NaN),
          c(6.60,  6.51, 6.39, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(6.09,  6.25, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(5.77,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
    }
    if (q==2){
      if (signif == 1){
        cv<- rbind(
          c(8.41,  8.01, 7.30, 6.85, 6.35, 5.97, 5.58, 5.32, 5.05),
          c(8.74,  8.36, 7.69, 7.20, 6.73, 6.32, 6.01, 5.69, NaN),
          c(8.47,  8.17, 7.64, 7.23, 6.81, 6.44, 6.11, NaN,  NaN),
          c(8.27,  7.97, 7.54, 7.19, 6.81, 6.52, NaN,  NaN,  NaN),
          c(7.88,  7.74, 7.40, 7.04, 6.74, NaN,  NaN,  NaN,  NaN),
          c(7.56,  7.49, 7.19, 6.92, NaN,  NaN,  NaN,  NaN,  NaN),
          c(7.29,  7.21, 6.99, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(6.98,  6.94, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(6.74,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(9.44,  8.77, 7.94, 7.39, 6.85, 6.43, 6.00, 5.70, 5.41),
          c(9.48,  9.02, 8.22, 7.70, 7.18, 6.73, 6.38, 6.03, NaN),
          c(9.16,  8.73, 8.18, 7.70, 7.23, 6.83, 6.47, NaN,  NaN),
          c(8.88,  8.51, 7.99, 7.59, 7.20, 6.87, NaN,  NaN,  NaN),
          c(8.40,  8.22, 7.80, 7.44, 7.12, NaN,  NaN,  NaN,  NaN),
          c(8.07,  7.88, 7.63, 7.31, NaN,  NaN,  NaN,  NaN,  NaN),
          c(7.75,  7.62, 7.37, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(7.42,  7.34, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(7.13,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(10.33,  9.44, 8.60, 7.85, 7.32, 6.84, 6.36, 6.04, 5.71),
          c(10.14,  9.61, 8.75, 8.16, 7.63, 7.08, 6.71, 6.36, NaN),
          c(9.79,  9.24, 8.65, 8.11, 7.61, 7.16, 6.82, NaN,  NaN),
          c(9.40,  9.02, 8.43, 7.95, 7.56, 7.20, NaN,  NaN,  NaN),
          c(8.85,  8.67, 8.21, 7.83, 7.45, NaN,  NaN,  NaN,  NaN),
          c(8.47,  8.27, 8.00, 7.65, NaN,  NaN,  NaN,  NaN,  NaN),
          c(8.14,  7.97, 7.70, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(7.79,  7.70, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(7.48,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(11.68, 10.26, 9.24, 8.44, 7.82, 7.41, 6.80, 6.47, 6.10),
          c(10.97, 10.42, 9.35, 8.76, 8.13, 7.51, 7.10, 6.75, NaN),
          c(10.54,  9.85, 9.22, 8.62, 8.05, 7.59, 7.18, NaN,  NaN),
          c(10.05,  9.63, 8.97, 8.42, 8.01, 7.61, NaN,  NaN,  NaN),
          c(9.44,  9.21, 8.71, 8.30, 7.89, NaN,  NaN,  NaN,  NaN),
          c(9.05,  8.76, 8.34, 8.08, NaN,  NaN,  NaN,  NaN,  NaN),
          c(8.56,  8.46, 8.10, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(8.26,  8.09, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(7.92,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
    }
    if (q==3){
      if (signif == 1){
        cv<- rbind(
          c(9.57,  8.72, 7.82, 7.24, 6.70, 6.26, 5.90, 5.53, 5.23),
          c(10.23,  9.38, 8.48, 7.90, 7.33, 6.87, 6.45, 6.09, NaN),
          c(10.05,  9.46, 8.68, 8.08, 7.61, 7.13, 6.75, NaN,  NaN),
          c(9.86,  9.36, 8.67, 8.15, 7.72, 7.26, NaN,  NaN,  NaN),
          c(9.49,  9.15, 8.58, 8.14, 7.73, NaN,  NaN,  NaN,  NaN),
          c(9.27,  8.96, 8.47, 8.08, NaN,  NaN,  NaN,  NaN,  NaN),
          c(8.95,  8.72, 8.29, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(8.70,  8.50, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(8.39,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(10.63,  9.51, 8.48, 7.79, 7.20, 6.69, 6.31, 5.90, 5.59),
          c(11.04, 10.03, 9.04, 8.41, 7.78, 7.28, 6.79, 6.47, NaN),
          c(10.79, 10.01, 9.18, 8.55, 8.06, 7.56, 7.13, NaN,  NaN),
          c(10.47,  9.93, 9.16, 8.60, 8.12, 7.64, NaN,  NaN,  NaN),
          c(10.09,  9.68, 9.04, 8.59, 8.10, NaN,  NaN,  NaN,  NaN),
          c(9.77,  9.43, 8.90, 8.49, NaN,  NaN,  NaN,  NaN,  NaN),
          c(9.45,  9.15, 8.75, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(9.15,  8.89, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(8.80,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(11.63, 10.25, 9.08, 8.26, 7.62, 7.12, 6.70, 6.25, 5.87),
          c(11.79, 10.67, 9.56, 8.88, 8.19, 7.66, 7.16, 6.80, NaN),
          c(11.45, 10.55, 9.71, 8.97, 8.43, 7.96, 7.48, NaN,  NaN),
          c(11.05, 10.42, 9.64, 9.03, 8.47, 7.98, NaN,  NaN,  NaN),
          c(10.61, 10.19, 9.45, 8.95, 8.44, NaN,  NaN,  NaN,  NaN),
          c(10.26,  9.81, 9.31, 8.80, NaN,  NaN,  NaN,  NaN,  NaN),
          c(9.88,  9.56, 9.13, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(9.52,  9.26, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(9.18,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(12.91, 11.08,  9.83, 8.93, 8.18, 7.71, 7.14, 6.65, 6.31),
          c(12.68, 11.43, 10.21, 9.46, 8.71, 8.08, 7.60, 7.18, NaN),
          c(12.28, 11.25, 10.35, 9.52, 8.94, 8.35, 7.88, NaN,  NaN),
          c(11.75, 11.05, 10.16, 9.50, 8.94, 8.45, NaN,  NaN,  NaN),
          c(11.17, 10.71,  9.96, 9.37, 8.84, NaN,  NaN,  NaN,  NaN),
          c(10.77, 10.34,  9.81, 9.26, NaN,  NaN,  NaN,  NaN,  NaN),
          c(10.46, 10.03,  9.59, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(10.06,  9.65,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(9.65,  NaN,   NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
    }
    if (q==4){
      if (signif == 1){
        cv<- rbind(
          c(10.61,  9.43, 8.29, 7.60, 7.03, 6.53, 6.14, 5.77, 5.41),
          c(11.54, 10.41, 9.31, 8.49, 7.87, 7.38, 6.88, 6.46, NaN),
          c(11.50, 10.55, 9.65, 8.88, 8.25, 7.75, 7.28, NaN,  NaN),
          c(11.30, 10.57, 9.71, 9.09, 8.54, 8.03, NaN,  NaN,  NaN),
          c(11.11, 10.44, 9.74, 9.14, 8.62, NaN,  NaN,  NaN,  NaN),
          c(10.77, 10.23, 9.69, 9.12, NaN,  NaN,  NaN,  NaN,  NaN),
          c(10.48, 10.08, 9.54, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(10.24,  9.86, NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(9.94,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(11.66, 10.20,  8.96, 8.18, 7.56, 6.97, 6.55, 6.14, 5.75),
          c(12.42, 11.07,  9.89, 9.01, 8.35, 7.81, 7.31, 6.81, NaN),
          c(12.26, 11.16, 10.20, 9.33, 8.70, 8.13, 7.66, NaN,  NaN),
          c(11.99, 11.14, 10.21, 9.51, 8.99, 8.42, NaN,  NaN,  NaN),
          c(11.70, 11.00, 10.23, 9.60, 9.02, NaN,  NaN,  NaN,  NaN),
          c(11.31, 10.75, 10.18, 9.52, NaN,  NaN,  NaN,  NaN,  NaN),
          c(10.99, 10.54,  9.99, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(10.71, 10.31,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(10.40,  NaN,   NaN,  NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(12.72, 10.90,  9.57, 8.68, 8.04, 7.37, 6.92, 6.48, 6.12),
          c(13.21, 11.69, 10.48, 9.48, 8.79, 8.16, 7.69, 7.13, NaN),
          c(13.00, 11.73, 10.72, 9.80, 9.09, 8.52, 8.00, NaN,  NaN),
          c(12.60, 11.65, 10.66, 9.93, 9.39, 8.78, NaN,  NaN,  NaN),
          c(12.24, 11.51, 10.62, 9.98, 9.38, NaN,  NaN,  NaN,  NaN),
          c(11.80, 11.27, 10.62, 9.92, NaN,  NaN,  NaN,  NaN,  NaN),
          c(11.52, 10.93, 10.38, NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(11.18, 10.67, NaN,   NaN,  NaN,  NaN,  NaN,  NaN,  NaN),
          c(10.77,  NaN,  NaN,   NaN,  NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(13.82, 11.79, 10.31,  9.30, 8.69, 7.87, 7.37, 6.88, 6.53),
          c(14.26, 12.40, 11.10, 10.15, 9.31, 8.66, 8.19, 7.58, NaN),
          c(13.79, 12.43, 11.26, 10.38, 9.58, 8.94, 8.44, NaN,  NaN),
          c(13.24, 12.28, 11.20, 10.41, 9.90, 9.25, NaN,  NaN,  NaN),
          c(12.91, 12.06, 11.19, 10.56, 9.80, NaN,  NaN,  NaN,  NaN),
          c(12.37, 11.85, 11.17, 10.36, NaN,  NaN,  NaN,  NaN,  NaN),
          c(12.11, 11.46, 10.87, NaN,   NaN,  NaN,  NaN,  NaN,  NaN),
          c(11.81, 11.28, NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  NaN),
          c(11.29,  NaN,  NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  NaN))
      }
    }
    if (q==5){
      if (signif == 1){
        cv<- rbind(
          c(11.51, 10.05,  8.80,  7.97, 7.32, 6.83, 6.38, 5.96, 5.62),
          c(12.73, 11.36, 10.03,  9.11, 8.37, 7.84, 7.27, 6.81, NaN),
          c(12.81, 11.62, 10.46,  9.65, 8.91, 8.30, 7.80, NaN,  NaN),
          c(12.72, 11.77, 10.73,  9.99, 9.28, 8.68, NaN,  NaN,  NaN),
          c(12.47, 11.71, 10.78, 10.09, 9.45, NaN,  NaN,  NaN,  NaN),
          c(12.30, 11.54, 10.79, 10.15, NaN,  NaN,  NaN,  NaN,  NaN),
          c(11.96, 11.37, 10.70, NaN,   NaN,  NaN,  NaN,  NaN,  NaN),
          c(11.68, 11.18, NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  NaN),
          c(11.38, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(12.67, 10.84,  9.46,  8.51, 7.83, 7.30, 6.78, 6.35, 5.96),
          c(13.73, 12.12, 10.64,  9.68, 8.85, 8.26, 7.67, 7.20, NaN),
          c(13.63, 12.27, 10.99, 10.18, 9.38, 8.76, 8.18, NaN,  NaN),
          c(13.43, 12.35, 11.27, 10.46, 9.73, 9.07, NaN,  NaN,  NaN),
          c(13.12, 12.33, 11.25, 10.58, 9.85, NaN,  NaN,  NaN,  NaN),
          c(12.89, 12.09, 11.24, 10.60, NaN,  NaN,  NaN,  NaN,  NaN),
          c(12.50, 11.89, 11.16, NaN,   NaN,  NaN,  NaN,  NaN,  NaN),
          c(12.16, 11.65, NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  NaN),
          c(11.85,  NaN,  NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(13.73, 11.61, 10.08,  9.08,  8.30, 7.70, 7.18, 6.71, 6.30),
          c(14.56, 12.80, 11.21, 10.19,  9.26, 8.67, 8.04, 7.52, NaN),
          c(14.36, 12.80, 11.48, 10.63,  9.78, 9.12, 8.61, NaN,  NaN),
          c(14.09, 12.92, 11.79, 10.90, 10.14, 9.41, NaN,  NaN,  NaN),
          c(13.69, 12.76, 11.72, 10.99, 10.23, NaN,  NaN,  NaN,  NaN),
          c(13.40, 12.55, 11.70, 11.00, NaN,   NaN,  NaN,  NaN,  NaN),
          c(13.02, 12.31, 11.62, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(12.59, 12.08, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(12.28, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(14.99, 12.66, 10.76,  9.72,  8.85, 8.22, 7.64, 7.17, 6.70),
          c(15.59, 13.55, 11.81, 10.79,  9.75, 9.18, 8.48, 7.91, NaN),
          c(15.31, 13.58, 12.09, 11.14, 10.24, 9.56, 9.01, NaN,  NaN),
          c(14.89, 13.60, 12.33, 11.39, 10.60, 9.88, NaN,  NaN,  NaN),
          c(14.43, 13.36, 12.32, 11.47, 10.71, NaN,  NaN,  NaN,  NaN),
          c(14.09, 13.10, 12.19, 11.43, NaN,   NaN,  NaN,  NaN,  NaN),
          c(13.63, 12.84, 12.10, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(13.13, 12.52, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(12.76, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
    }
    if (q==6){
      if (signif == 1){
        cv<- rbind(
          c(12.41,  10.66, 9.18,  8.33,  7.61, 7.05, 6.58, 6.15, 5.81),
          c(13.88,  12.19, 10.69, 9.71,  8.86, 8.19, 7.65, 7.18, NaN),
          c(14.08,  12.62, 11.35, 10.37, 9.61, 8.91, 8.31, NaN,  NaN),
          c(14.12,  12.79, 11.70, 10.75,10.02, 9.32, NaN,  NaN,  NaN),
          c(13.91,  12.84, 11.82, 11.02,10.28, NaN,  NaN,  NaN,  NaN),
          c(13.64,  12.77, 11.86, 11.13, NaN,  NaN,  NaN,  NaN,  NaN),
          c(13.35,  12.58, 11.80, NaN,   NaN,  NaN,  NaN,  NaN,  NaN),
          c(13.07,  12.44, NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  NaN),
          c(12.76,  NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(13.59,  11.49, 9.86,  8.90,  8.13,  7.51, 7.00, 6.55, 6.17),
          c(14.80,  12.92, 11.32, 10.26, 9.36,  8.65, 8.04, 7.54, NaN),
          c(14.94,  13.31, 11.87, 10.90, 10.08, 9.33, 8.73, NaN,  NaN),
          c(14.80,  13.42, 12.28, 11.25, 10.44, 9.76, NaN,  NaN,  NaN),
          c(14.58,  13.42, 12.31, 11.48, 10.69, NaN,  NaN,  NaN,  NaN),
          c(14.28,  13.32, 12.38, 11.60, NaN,   NaN,  NaN,  NaN,  NaN),
          c(13.91,  13.07, 12.28, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(13.62,  12.95, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(13.25,  NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(14.57,  12.25, 10.54, 9.44,  8.62,  7.90, 7.41, 6.94, 6.49),
          c(15.61,  13.55, 11.86, 10.71, 9.83,  9.06, 8.42, 7.90, NaN),
          c(15.66,  13.96, 12.42, 11.35, 10.48, 9.68, 9.11, NaN,  NaN),
          c(15.43,  13.98, 12.79, 11.68, 10.85, 10.12,NaN,  NaN,  NaN),
          c(15.16,  13.94, 12.83, 11.91, 11.11, NaN,  NaN,  NaN,  NaN),
          c(14.82,  13.78, 12.85, 12.03, NaN,   NaN,  NaN,  NaN,  NaN),
          c(14.37,  13.58, 12.70, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(14.11,  13.38, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(13.77,  NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(15.61,  13.25, 11.42, 10.04, 9.16,  8.42,  7.85, 7.41, 6.89),
          c(16.71,  14.41, 12.58, 11.36, 10.35, 9.57,  8.81, 8.38, NaN),
          c(16.54,  14.71, 13.07, 11.89, 10.91, 10.18, 9.60, NaN,  NaN),
          c(16.25,  14.74, 13.45, 12.22, 11.42, 10.58, NaN,  NaN,  NaN),
          c(15.80,  14.50, 13.39, 12.43, 11.57, NaN,   NaN,  NaN,  NaN),
          c(15.43,  14.37, 13.40, 12.60, NaN,   NaN,   NaN,  NaN,  NaN),
          c(14.98,  14.13, 13.20, NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(14.70,  13.88, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(14.36,  NaN,   NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
    }
    if (q==7){
      if (signif == 1){
        cv<- rbind(
          c(13.28, 11.26, 9.68,  8.65,  7.89,  7.29, 6.77, 6.33, 5.91),
          c(15.03, 12.97, 11.35, 10.20, 9.31,  8.61, 7.98, 7.45, NaN),
          c(15.31, 13.65, 12.22, 11.10, 10.16, 9.41, 8.84, NaN,  NaN),
          c(15.34, 13.94, 12.63, 11.56, 10.74, 10.01,NaN,  NaN,  NaN),
          c(15.21, 14.04, 12.82, 11.91, 11.09, NaN,  NaN,  NaN,  NaN),
          c(14.97, 13.92, 12.89, 12.02, NaN,   NaN,  NaN,  NaN,  NaN),
          c(14.69, 13.82, 12.89, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(14.42, 13.65, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(14.17, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(14.55, 12.10, 10.35, 9.25,  8.43,  7.75, 7.20, 6.70, 6.31),
          c(15.99, 13.75, 12.01, 10.78, 9.81,  9.07, 8.36, 7.86, NaN),
          c(16.14, 14.37, 12.80, 11.62, 10.62, 9.88, 9.26, NaN,  NaN),
          c(16.14, 14.59, 13.22, 12.07, 11.19, 10.45,NaN,  NaN,  NaN),
          c(15.88, 14.66, 13.39, 12.40, 11.53, NaN,  NaN,  NaN,  NaN),
          c(15.60, 14.46, 13.37, 12.51, NaN,   NaN,  NaN,  NaN,  NaN),
          c(15.31, 14.39, 13.37, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(14.99, 14.15, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(14.71, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(15.63, 12.85, 11.06, 9.78,  8.89,  8.22, 7.58, 7.08, 6.66),
          c(16.86, 14.44, 12.59, 11.29, 10.25, 9.46, 8.72, 8.17, NaN),
          c(16.87, 14.98, 13.32, 12.07, 11.06, 10.28,9.57, NaN,  NaN),
          c(16.89, 15.17, 13.75, 12.56, 11.65, 10.82,NaN,  NaN,  NaN),
          c(16.47, 15.21, 13.88, 12.86, 11.94, NaN,  NaN,  NaN,  NaN),
          c(16.13, 15.00, 13.83, 12.87, NaN,   NaN,  NaN,  NaN,  NaN),
          c(15.84, 14.90, 13.80, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(15.53, 14.59, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(15.16, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(17.11, 13.80, 11.99, 10.40, 9.44,  8.70,  8.10, 7.55, 7.01),
          c(18.15, 15.19, 13.30, 11.88, 10.79, 10.01, 9.20, 8.54, NaN),
          c(17.82, 15.78, 14.00, 12.62, 11.50, 10.86, 10.09,NaN,  NaN),
          c(17.53, 15.86, 14.46, 13.09, 12.16, 11.37, NaN,  NaN,  NaN),
          c(17.11, 15.90, 14.47, 13.37, 12.50, NaN,   NaN,  NaN,  NaN),
          c(16.76, 15.76, 14.32, 13.33, NaN,   NaN,   NaN,  NaN,  NaN),
          c(16.54, 15.50, 14.44, NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(16.18, 15.03, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(15.75,  NaN,  NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
    }
    if (q==8){
      if (signif == 1){
        cv<- rbind(
          c(14.16, 11.82, 10.02, 8.96,  8.20,  7.54, 6.95, 6.52, 6.11),
          c(16.14, 13.77, 12.00, 10.76, 9.76,  8.98, 8.34, 7.83, NaN),
          c(16.60, 14.60, 12.94, 11.75, 10.75, 9.99, 9.26, NaN,  NaN),
          c(16.60, 15.00, 13.50, 12.34, 11.42, 10.66,NaN,  NaN,  NaN),
          c(16.49, 15.10, 13.79, 12.74, 11.82, NaN,  NaN,  NaN,  NaN),
          c(16.29, 15.09, 13.93, 12.98, NaN,   NaN,  NaN,  NaN,  NaN),
          c(16.05, 14.99, 13.95, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(15.73, 14.81, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(15.50,  NaN,  NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(15.38, 12.64, 10.80, 9.59,  8.70,  7.99, 7.40, 6.95, 6.47),
          c(17.11, 14.55, 12.64, 11.31, 10.27, 9.47, 8.79, 8.24, NaN),
          c(17.51, 15.32, 13.60, 12.28, 11.26, 10.46,9.69, NaN,  NaN),
          c(17.39, 15.68, 14.04, 12.93, 11.92, 11.10,NaN,  NaN,  NaN),
          c(17.18, 15.75, 14.35, 13.22, 12.28, NaN,  NaN,  NaN,  NaN),
          c(16.95, 15.71, 14.45, 13.47, NaN,   NaN,  NaN,  NaN,  NaN),
          c(16.66, 15.53, 14.47, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(16.36, 15.33, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(16.02, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(16.62, 13.47, 11.50, 10.20, 9.22,  8.43, 7.77, 7.27, 6.80),
          c(18.15, 15.30, 13.24, 11.83, 10.75, 9.88, 9.18, 8.57, NaN),
          c(18.27, 16.00, 14.10, 12.81, 11.67, 10.89,10.09,NaN,  NaN),
          c(18.07, 16.29, 14.58, 13.40, 12.38, 11.46,NaN,  NaN,  NaN),
          c(17.82, 16.29, 14.89, 13.64, 12.68, NaN,  NaN,  NaN,  NaN),
          c(17.53, 16.25, 14.92, 13.90, NaN,   NaN,  NaN,  NaN,  NaN),
          c(17.18, 16.00, 14.97, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(16.88, 15.77, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(16.48,  NaN,  NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(18.21, 14.40, 12.33, 10.77, 9.86,  8.93, 8.25, 7.75, 7.27),
          c(19.26, 16.11, 13.91, 12.41, 11.32, 10.35,9.53, 8.95, NaN),
          c(19.22, 16.75, 14.77, 13.46, 12.20, 11.32,10.57,NaN,  NaN),
          c(19.07, 17.06, 15.24, 13.97, 12.89, 12.01,NaN,  NaN,  NaN),
          c(18.61, 17.05, 15.54, 14.15, 13.22, NaN,  NaN,  NaN,  NaN),
          c(18.20, 16.90, 15.60, 14.41, NaN,   NaN,  NaN,  NaN,  NaN),
          c(17.87, 16.55, 15.54, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(17.43, 16.25, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(17.13,  NaN,  NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
    }
    if (q==9){
      if (signif == 1){
        cv<- rbind(
          c(14.91, 12.31, 10.45, 9.34,  8.42,  7.75, 7.17, 6.71, 6.26),
          c(17.21, 14.56, 12.66, 11.28, 10.22, 9.42, 8.73, 8.12, NaN),
          c(17.74, 15.52, 13.71, 12.40, 11.35, 10.46,9.75, NaN,  NaN),
          c(17.82, 16.01, 14.40, 13.11, 12.06, 11.25,NaN,  NaN,  NaN),
          c(17.72, 16.19, 14.79, 13.57, 12.62, NaN,  NaN,  NaN,  NaN),
          c(17.58, 16.16, 14.92, 13.83, NaN,   NaN,  NaN,  NaN,  NaN),
          c(17.29, 16.17, 14.99, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(17.05, 15.96, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(16.72,  NaN,  NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(16.19, 13.25, 11.20, 9.93,  8.96,  8.24,  7.64, 7.11, 6.61),
          c(18.19, 15.36, 13.35, 11.85, 10.74, 9.91,  9.13, 8.51, NaN),
          c(18.60, 16.27, 14.33, 12.94, 11.82, 10.93, 10.18,NaN,  NaN),
          c(18.65, 16.72, 15.01, 13.65, 12.57, 11.69, NaN,  NaN,  NaN),
          c(18.48, 16.84, 15.33, 14.08, 13.09, NaN,   NaN,  NaN,  NaN),
          c(18.21, 16.78, 15.47, 14.31, NaN,   NaN,   NaN,  NaN,  NaN),
          c(17.92, 16.73, 15.53, NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(17.66, 16.53, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(17.26, NaN,   NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(17.36, 14.01, 11.93, 10.55, 9.43,  8.69,  8.01, 7.48, 6.96),
          c(19.06, 16.05, 13.95, 12.38, 11.18, 10.32, 9.54, 8.87, NaN),
          c(19.41, 16.97, 14.92, 13.47, 12.33, 11.36, 10.49,NaN,  NaN),
          c(19.28, 17.35, 15.55, 14.16, 13.03, 12.08, NaN,  NaN,  NaN),
          c(19.17, 17.38, 15.84, 14.60, 13.52, NaN,   NaN,  NaN,  NaN),
          c(18.81, 17.29, 16.00, 14.77, NaN,   NaN,   NaN,  NaN,  NaN),
          c(18.46, 17.28, 16.06, NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(18.12, 17.00, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(17.78,  NaN,  NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(18.89,  15.11, 12.82, 11.12, 10.02, 9.23,  8.50, 7.97, 7.36),
          c(20.09,  16.97, 14.82, 13.04, 11.74, 10.89, 9.99, 9.31, NaN),
          c(20.46,  17.79, 15.56, 14.06, 12.95, 11.76, 10.98,NaN,  NaN),
          c(20.14,  18.10, 16.17, 14.77, 13.56, 12.54, NaN,  NaN,  NaN),
          c(19.96,  18.17, 16.42, 15.12, 13.98, NaN,   NaN,  NaN,  NaN),
          c(19.55,  17.89, 16.62, 15.26, NaN,   NaN,   NaN,  NaN,  NaN),
          c(19.11,  18.00, 16.65, NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(18.80,  17.58, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(18.52,  NaN,   NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
    }
    if (q==10){
      if (signif == 1){
        cv<- rbind(
          c(15.92, 12.83, 10.84, 9.67,  8.72,  7.99, 7.38, 6.89, 6.44),
          c(18.15, 15.36, 13.22, 11.79, 10.70, 9.80, 9.09, 8.41, NaN),
          c(18.83, 16.42, 14.52, 13.04, 11.95, 10.98,10.20,NaN,  NaN),
          c(19.12, 16.98, 15.27, 13.84, 12.75, 11.86,NaN,  NaN,  NaN),
          c(18.95, 17.26, 15.63, 14.41, 13.34, NaN,  NaN,  NaN,  NaN),
          c(18.80, 17.33, 15.88, 14.72, NaN,   NaN,  NaN,  NaN,  NaN),
          c(18.56, 17.26, 15.98, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(18.28, 17.15, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(18.01, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(17.18, 13.77, 11.62, 10.29, 9.31,  8.51, 7.86, 7.28, 6.82),
          c(19.23, 16.18, 13.92, 12.41, 11.28, 10.26,9.54, 8.82, NaN),
          c(19.69, 17.19, 15.14, 13.63, 12.46, 11.47,10.62,NaN,  NaN),
          c(19.95, 17.69, 15.87, 14.38, 13.27, 12.32,NaN,  NaN,  NaN),
          c(19.70, 17.86, 16.23, 14.94, 13.87, NaN,  NaN,  NaN,  NaN),
          c(19.48, 17.94, 16.45, 15.23, NaN,   NaN,  NaN,  NaN,  NaN),
          c(19.21, 17.87, 16.50, NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(18.89, 17.74, NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN),
          c(18.60, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(18.28, 14.60, 12.32, 10.91, 9.78,  8.98,  8.26, 7.63, 7.14),
          c(20.14, 16.93, 14.50, 12.91, 11.73, 10.77, 9.89, 9.19, NaN),
          c(20.56, 17.88, 15.75, 14.23, 12.93, 11.83, 11.02,NaN,  NaN),
          c(20.71, 18.38, 16.44, 14.90, 13.71, 12.69, NaN,  NaN,  NaN),
          c(20.40, 18.48, 16.77, 15.42, 14.30, NaN,   NaN,  NaN,  NaN),
          c(20.04, 18.49, 17.01, 15.68, NaN,   NaN,   NaN,  NaN,  NaN),
          c(19.86, 18.40, 16.98, NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(19.47, 18.24, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(19.10,  NaN,  NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(19.87, 15.60, 13.14, 11.55, 10.38, 9.47,  8.77, 8.07, 7.53),
          c(21.12, 17.83, 15.29, 13.53, 12.37, 11.27, 10.35,9.57, NaN),
          c(21.61, 18.79, 16.49, 14.88, 13.47, 12.44, 11.47,NaN,  NaN),
          c(21.59, 19.10, 17.04, 15.44, 14.24, 13.23, NaN,  NaN,  NaN),
          c(21.25, 19.10, 17.39, 15.90, 14.82, NaN,   NaN,  NaN,  NaN),
          c(20.80, 19.16, 17.66, 16.29, NaN,   NaN,   NaN,  NaN,  NaN),
          c(20.59, 19.06, 17.63, NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(20.08, 18.83, NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN),
          c(19.74, NaN,   NaN,   NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      } 
    }
  }
  
  if (trm == .10){
    cv <- matrix(0,7,7)
    if (q==1){
      if (signif == 1){
        cv<- rbind(
          c(6.59,   6.32,   5.77,   5.24,   4.78,   4.34,   3.93), 
          c(6.34,   6.20,   5.74,   5.38,   4.91,   4.49,    NaN),
          c(5.79,   5.74,   5.43,   5.11,   4.75,    NaN,    NaN), 
          c(5.24,   5.32,   5.11,   4.86,    NaN,    NaN,    NaN),  
          c(4.76,   4.96,   4.79,    NaN,    NaN,    NaN,    NaN),   
          c(4.35,   4.49,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(3.91,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(7.63,   7.10,   6.43,   5.76,   5.26,   4.73,  4.32), 
          c(7.12,   6.83,   6.26,   5.86,   5.35,   4.90,   NaN),
          c(6.39,   6.27,   5.90,   5.54,   5.12,    NaN,   NaN), 
          c(5.78,   5.81,   5.52,   5.25,    NaN,    NaN,   NaN),  
          c(5.19,   5.33,   5.16,    NaN,    NaN,    NaN,   NaN),   
          c(4.81,   4.84,    NaN,    NaN,    NaN,    NaN,   NaN),   
          c(4.28,    NaN,    NaN,    NaN,    NaN,    NaN,   NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(8.54,   7.75,   6.93,   6.26,   5.69,   5.10,  4.73), 
          c(7.78,   7.44,   6.77,   6.24,   5.76,   5.27,   NaN),
          c(6.99,   6.72,   6.31,   5.86,   5.45,    NaN,   NaN), 
          c(6.28,   6.25,   5.92,   5.61,    NaN,    NaN,   NaN),  
          c(5.61,   5.72,   5.52,    NaN,    NaN,    NaN,   NaN),   
          c(5.19,   5.18,    NaN,    NaN,    NaN,    NaN,   NaN),   
          c(4.59,    NaN,    NaN,    NaN,    NaN,    NaN,   NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(9.79,   8.70,   7.77,   6.86,   6.16,   5.48,  5.16), 
          c(8.73,   8.17,   7.36,   6.76,   6.18,   5.72,   NaN),
          c(7.60,   7.25,   6.83,   6.30,   5.91,    NaN,   NaN), 
          c(6.87,   6.81,   6.32,   6.07,    NaN,    NaN,   NaN),  
          c(6.09,   6.17,   5.92,    NaN,    NaN,    NaN,   NaN),   
          c(5.62,   5.65,    NaN,    NaN,    NaN,    NaN,   NaN),   
          c(4.97,    NaN,    NaN,    NaN,    NaN,    NaN,   NaN))
      }
    }  
    if (q==2){
      if (signif == 1){
        cv<- rbind(
          c(7.88,   7.18,   6.40,   5.78,   5.25,   4.72,   4.25), 
          c(7.96,   7.41,   6.70,   6.10,   5.53,   5.08,    NaN),
          c(7.52,   7.12,   6.59,   6.06,   5.63,    NaN,    NaN), 
          c(7.03,   6.82,   6.35,   5.93,    NaN,    NaN,    NaN),  
          c(6.64,   6.42,   6.08,    NaN,    NaN,    NaN,    NaN),   
          c(6.11,   6.03,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(5.65,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(8.87,   7.94,   6.96,   6.34,   5.73,   5.16,   4.63), 
          c(8.78,   8.03,   7.23,   6.58,   5.98,   5.48,    NaN),
          c(8.15,   7.68,   7.08,   6.52,   6.04,    NaN,    NaN), 
          c(7.63,   7.32,   6.81,   6.32,    NaN,    NaN,    NaN),  
          c(7.14,   6.88,   6.51,    NaN,    NaN,    NaN,    NaN),   
          c(6.58,   6.49,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(6.04,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(9.85,   8.69,   7.54,   6.80,   6.16,   5.58,   4.94), 
          c(9.52,   8.69,   7.77,   7.02,   6.36,   5.85,    NaN),
          c(8.73,   8.16,   7.48,   6.96,   6.44,    NaN,    NaN), 
          c(8.20,   7.78,   7.23,   6.69,    NaN,    NaN,    NaN),  
          c(7.58,   7.26,   6.88,    NaN,    NaN,    NaN,    NaN),   
          c(7.03,   6.85,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(6.44,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(11.12,   9.52,   8.37,   7.43,   6.77,   6.04,   5.35), 
          c(10.55,   9.52,   8.37,   7.57,   6.89,   6.27,    NaN),
          c(9.68,   8.80,   8.01,   7.54,   6.94,    NaN,    NaN), 
          c(8.82,   8.35,   7.81,   7.27,    NaN,    NaN,    NaN),  
          c(8.25,   7.82,   7.47,    NaN,    NaN,    NaN,    NaN),   
          c(7.60,   7.39,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(6.98,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
    }  
    if (q==3){
      if (signif == 1){
        cv<- rbind(
          c(8.98,   7.93,   6.94,   6.20,   5.54,   5.04,   4.49), 
          c(9.34,   8.44,   7.52,   6.76,   6.18,   5.61,    NaN),
          c(9.04,   8.33,   7.56,   6.94,   6.35,    NaN,    NaN), 
          c(8.62,   8.05,   7.47,   6.87,    NaN,    NaN,    NaN),  
          c(8.12,   7.75,   7.27,    NaN,    NaN,    NaN,    NaN),   
          c(7.66,   7.37,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(7.08,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(10.06,   8.72,   7.58,   6.76,   6.03,   5.46,   4.89), 
          c(10.23,   9.11,   8.09,   7.24,   6.63,   6.02,    NaN),
          c(9.73,   8.91,   8.09,   7.39,   6.77,    NaN,    NaN), 
          c(9.29,   8.61,   7.99,   7.31,    NaN,    NaN,    NaN),  
          c(8.69,   8.26,   7.73,    NaN,    NaN,    NaN,    NaN),   
          c(8.22,   7.84,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(7.54,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(11.08,   9.43,   8.18,   7.30,   6.47,   5.82,   5.28), 
          c(10.98,   9.75,   8.60,   7.68,   7.04,   6.40,    NaN),
          c(10.36,   9.46,   8.57,   7.77,   7.18,    NaN,    NaN), 
          c(9.87,   9.06,   8.45,   7.71,    NaN,    NaN,    NaN),  
          c(9.20,   8.74,   8.10,    NaN,    NaN,    NaN,    NaN),   
          c(8.73,   8.27,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(7.97,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(12.43,  10.33,   8.91,   7.88,   7.01,   6.34,   5.66), 
          c(12.01,  10.53,   9.20,   8.27,   7.50,   6.83,    NaN),
          c(11.18,  10.09,   9.15,   8.31,   7.68,    NaN,    NaN), 
          c(10.58,   9.80,   8.98,   8.29,    NaN,    NaN,    NaN),  
          c(9.80,   9.34,   8.66,    NaN,    NaN,    NaN,    NaN),   
          c(9.23,   8.81,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(8.48,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
    }  
    if (q==4){
      if (signif == 1){
        cv<- rbind(
          c(9.96,   8.54,   7.41,   6.54,   5.93,   5.27,   4.76), 
          c(10.60,   9.32,   8.25,   7.37,   6.63,   6.01,    NaN),
          c(10.39,   9.43,   8.48,   7.73,   7.00,    NaN,    NaN), 
          c(10.01,   9.27,   8.49,   7.80,    NaN,    NaN,    NaN),  
          c(9.59,   8.94,   8.27,    NaN,    NaN,    NaN,    NaN),   
          c(9.02,   8.60,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(8.45,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(11.10,   9.38,   8.08,   7.09,   6.42,   5.68,   5.14), 
          c(11.51,  10.05,   8.85,   7.88,   7.10,   6.41,    NaN),
          c(11.14,  10.06,   9.02,   8.29,   7.45,    NaN,    NaN), 
          c(10.72,   9.87,   9.04,   8.28,    NaN,    NaN,    NaN),  
          c(10.20,   9.46,   8.74,    NaN,    NaN,    NaN,    NaN),   
          c(9.60,   9.07,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(8.93,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(12.17,  10.13,   8.72,   7.64,   6.89,   6.10,   5.49), 
          c(12.30,  10.72,   9.39,   8.38,   7.55,   6.78,    NaN),
          c(11.79,  10.56,   9.52,   8.76,   7.85,    NaN,    NaN), 
          c(11.33,  10.43,   9.47,   8.72,    NaN,    NaN,    NaN),  
          c(10.79,   9.97,   9.22,    NaN,    NaN,    NaN,    NaN),   
          c(10.13,   9.55,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(9.37,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(13.50,  11.07,   9.44,   8.33,   7.41,   6.64,   5.98), 
          c(13.36,  11.59,  10.02,   8.99,   8.03,   7.30,    NaN),
          c(12.62,  11.30,  10.09,   9.31,   8.29,    NaN,    NaN), 
          c(12.18,  11.10,  10.08,   9.25,    NaN,    NaN,    NaN),  
          c(11.49,  10.61,   9.69,    NaN,    NaN,    NaN,    NaN),   
          c(10.62,  10.17,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(9.86,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN)) 
      }
    }  
    if (q==5){
      if (signif == 1){
        cv<- rbind(
          c(10.94,   9.19,   7.88,   6.94,   6.19,   5.52,   4.95), 
          c(11.81,  10.21,   8.93,   7.92,   7.17,   6.46,    NaN),
          c(11.74,  10.44,   9.29,   8.46,   7.67,    NaN,    NaN), 
          c(11.38,  10.36,   9.38,   8.56,    NaN,    NaN,    NaN),  
          c(10.86,  10.08,   9.29,    NaN,    NaN,    NaN,    NaN),   
          c(10.36,   9.77,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(9.76,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(12.14,  10.00,   8.59,   7.46,   6.66,   5.93,   5.31), 
          c(12.76,  10.99,   9.55,   8.44,   7.61,   6.87,    NaN),
          c(12.52,  11.10,   9.88,   8.97,   8.12,    NaN,    NaN), 
          c(12.06,  10.99,   9.93,   9.04,    NaN,    NaN,    NaN),  
          c(11.52,  10.63,   9.84,    NaN,    NaN,    NaN,    NaN),   
          c(10.95,  10.29,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(10.26,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(13.22,  10.74,   9.18,   7.97,   7.10,   6.33,   5.70), 
          c(13.68,  11.63,  10.08,   8.87,   8.05,   7.30,    NaN),
          c(13.13,  11.66,  10.36,   9.43,   8.51,    NaN,    NaN), 
          c(12.76,  11.59,  10.41,   9.47,    NaN,    NaN,    NaN),  
          c(12.20,  11.15,  10.32,    NaN,    NaN,    NaN,    NaN),   
          c(11.52,  10.75,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(10.80,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(14.47,  11.77,   9.83,   8.77,   7.62,   6.82,   6.13), 
          c(14.66,  12.50,  10.70,   9.51,   8.55,   7.76,    NaN),
          c(14.00,  12.31,  10.91,  10.05,   9.00,    NaN,    NaN), 
          c(13.52,  12.19,  11.10,  10.05,    NaN,    NaN,    NaN),  
          c(12.77,  11.80,  10.97,    NaN,    NaN,    NaN,    NaN),   
          c(12.20,  11.24,    NaN,    NaN,    NaN,    NaN,    NaN),   
          c(11.38,    NaN,    NaN,    NaN,    NaN,    NaN,    NaN))
      }
    }  
    if (q==6){
      if (signif == 1){
        cv<- rbind(
          c(11.76,  9.76,   8.25,  7.27,  6.45, 5.79, 5.16), 
          c(12.91,  11.08,  9.59,  8.53,  7.63, 6.90, NaN),
          c(12.94,  11.42,  10.12, 9.12,  8.25, NaN,  NaN), 
          c(12.60,  11.39,  10.32, 9.41,  NaN,  NaN,  NaN),  
          c(12.21,  11.24,  10.29, NaN,   NaN,  NaN,  NaN),   
          c(11.64,  10.89,  NaN,   NaN,   NaN,  NaN,  NaN),   
          c(11.03,  NaN,    NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(13.00,  10.56, 8.96,  7.85, 7.00, 6.25, 5.55), 
          c(13.93,  11.88, 10.23, 9.09, 8.15, 7.33, NaN),
          c(13.83,  12.11, 10.69, 9.65, 8.71, NaN,  NaN), 
          c(13.33,  12.02, 10.85, 9.93, NaN,  NaN,  NaN),  
          c(12.89,  11.85, 10.84, NaN,  NaN,  NaN,  NaN),   
          c(12.23,  11.49, NaN,   NaN,  NaN,  NaN,  NaN),   
          c(11.59,  NaN,   NaN,   NaN,  NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(14.09,  11.36, 9.62,  8.44,  7.48, 6.70, 5.91), 
          c(14.79,  12.56, 10.83, 9.59,  8.57, 7.70, NaN),
          c(14.60,  12.73, 11.24, 10.17, 9.17, NaN,  NaN), 
          c(14.04,  12.58, 11.37, 10.42, NaN,  NaN,  NaN),  
          c(13.55,  12.35, 11.34, NaN,   NaN,  NaN,  NaN),   
          c(12.81,  11.99, NaN,   NaN,   NaN,  NaN,  NaN),   
          c(12.09,  NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(15.24,  12.25,  10.41, 8.99,  7.99, 7.22, 6.37), 
          c(15.80,  13.40,  11.60, 10.17, 9.13, 8.24, NaN),
          c(15.60,  13.68,  11.79, 10.80, 9.66, NaN,  NaN), 
          c(14.75,  13.35,  11.94, 10.97, NaN,  NaN,  NaN),  
          c(14.23,  12.92,  11.87, NaN,   NaN,  NaN,  NaN),   
          c(13.46,  12.54,  NaN,   NaN,   NaN,  NaN,  NaN),   
          c(12.73,  NaN,    NaN,   NaN,   NaN,  NaN,  NaN))
      }
    }  
    if (q==7){
      if (signif == 1){
        cv<- rbind(
          c(12.63, 10.38, 8.73,  7.57,  6.76, 6.02, 5.35), 
          c(14.03, 11.88, 10.24, 9.08,  8.09, 7.29, NaN),
          c(14.10, 12.38, 10.92, 9.75,  8.82, NaN,  NaN), 
          c(13.88, 12.45, 11.19, 10.20, NaN,  NaN,  NaN),  
          c(13.41, 12.34, 11.17, NaN,   NaN,  NaN,  NaN),   
          c(12.95, 12.00, NaN,   NaN,   NaN,  NaN,  NaN),   
          c(12.24, NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(13.89, 11.20, 9.39,  8.15,  7.30, 6.49, 5.78), 
          c(14.97, 12.69, 10.88, 9.64,  8.62, 7.71, NaN),
          c(14.96, 13.08, 11.54, 10.31, 9.31, NaN,  NaN), 
          c(14.67, 13.11, 11.81, 10.72, NaN,  NaN,  NaN),  
          c(14.06, 12.96, 11.74, NaN,   NaN,  NaN,  NaN),   
          c(13.60, 12.55, NaN,   NaN,   NaN,  NaN,  NaN),   
          c(12.84, NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(15.08, 11.92, 10.10, 8.72,  7.73, 6.91, 6.19), 
          c(15.88, 13.44, 11.43, 10.11, 9.06, 8.09, NaN),
          c(15.71, 13.76, 12.12, 10.82, 9.74, NaN,  NaN), 
          c(15.40, 13.70, 12.27, 11.15, NaN,  NaN,  NaN),  
          c(14.64, 13.51, 12.24, NaN,   NaN,  NaN,  NaN),   
          c(14.19, 13.14, NaN,   NaN,   NaN,  NaN,  NaN),   
          c(13.28, NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(16.37, 12.87, 10.96, 9.36,  8.27,  7.42, 6.65), 
          c(16.96, 14.08, 12.22, 10.70, 9.61,  8.57, NaN),
          c(16.52, 14.57, 12.77, 11.41, 10.28, NaN,  NaN), 
          c(16.39, 14.42, 12.88, 11.67, NaN,   NaN,  NaN),  
          c(15.39, 14.12, 12.85, NaN,   NaN,   NaN,  NaN),   
          c(14.91, 13.76, NaN,   NaN,   NaN,   NaN,  NaN),   
          c(13.96, NaN,   NaN,   NaN,   NaN,   NaN,  NaN))
      }
    }  
    if (q==8){
      if (signif == 1){
        cv<- rbind(
          c(13.52, 10.89, 9.07,  7.91,  7.03, 6.26, 5.56), 
          c(15.14, 12.70, 10.88, 9.60,  8.55, 7.67, NaN),
          c(15.28, 13.27, 11.65, 10.48, 9.42, NaN,  NaN), 
          c(15.11, 13.43, 12.04, 10.90, NaN,  NaN,  NaN),  
          c(14.71, 13.37, 12.12, NaN,   NaN,  NaN,  NaN),   
          c(14.18, 13.06, NaN,   NaN,   NaN,  NaN,  NaN),   
          c(13.48, NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(14.73, 11.78, 9.82,  8.52,  7.57, 6.71, 5.99), 
          c(16.18, 13.53, 11.52, 10.18, 9.07, 8.12, NaN),
          c(16.22, 13.99, 12.25, 11.06, 9.91, NaN,  NaN), 
          c(15.86, 14.12, 12.63, 11.45, NaN,  NaN,  NaN),  
          c(15.45, 14.03, 12.67, NaN,   NaN,  NaN,  NaN),   
          c(14.86, 13.68, NaN,   NaN,   NaN,  NaN,  NaN),   
          c(14.15, NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(15.98, 12.54, 10.49, 9.10,  8.06,  7.13, 6.34), 
          c(17.10, 14.21, 12.20, 10.72, 9.52,  8.56, NaN),
          c(17.02, 14.61, 12.82, 11.55, 10.32, NaN,  NaN), 
          c(16.66, 14.77, 13.17, 11.92, NaN,   NaN,  NaN),  
          c(16.06, 14.61, 13.18, NaN,   NaN,   NaN,  NaN),   
          c(15.51, 14.24, NaN,   NaN,   NaN,   NaN,  NaN),   
          c(14.63, NaN,   NaN,   NaN,   NaN,   NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(17.47, 13.41, 11.29, 9.72,  8.63,  7.64, 6.80), 
          c(18.12, 15.17, 12.93, 11.37, 10.02, 9.03, NaN),
          c(18.06, 15.42, 13.48, 12.19, 10.92, NaN,  NaN), 
          c(17.53, 15.46, 13.84, 12.54, NaN,   NaN,  NaN),  
          c(16.81, 15.41, 13.78, NaN,   NaN,   NaN,  NaN),   
          c(16.30, 14.81, NaN,   NaN,   NaN,   NaN,  NaN),   
          c(15.31, NaN,   NaN,   NaN,   NaN,   NaN,  NaN))
      }
    }  
    if (q==9){
      if (signif == 1){
        cv<- rbind(
          c(14.30, 11.43, 9.53,  8.27,  7.24, 6.49, 5.75), 
          c(16.22, 13.46, 11.52, 10.10, 8.99, 8.06, NaN),
          c(16.36, 14.22, 12.44, 11.10, 9.98, NaN,  NaN), 
          c(16.24, 14.42, 12.85, 11.65, NaN,  NaN,  NaN),  
          c(15.85, 14.39, 12.98, NaN,   NaN,  NaN,  NaN),   
          c(15.36, 14.10, NaN,   NaN,   NaN,  NaN,  NaN),   
          c(14.67, NaN,   NaN,   NaN,   NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(15.53, 12.27, 10.25, 8.87,  7.75,  6.91, 6.19), 
          c(17.29, 14.31, 12.22, 10.72, 9.53,  8.54, NaN),
          c(17.24, 14.96, 13.07, 11.69, 10.51, NaN,  NaN), 
          c(17.08, 15.10, 13.49, 12.22, NaN,   NaN,  NaN),  
          c(16.58, 15.03, 13.62, NaN,   NaN,   NaN,  NaN),   
          c(16.08, 14.70, NaN,   NaN,   NaN,   NaN,  NaN),   
          c(15.31, NaN,   NaN,   NaN,   NaN,   NaN,  NaN)) 
      }
      if (signif == 3){
        cv<- rbind(
          c(16.72, 13.13, 10.99, 9.44,  8.27,  7.38, 6.53), 
          c(18.19, 15.06, 12.84, 11.24, 10.04, 8.94, NaN),
          c(18.10, 15.59, 13.65, 12.23, 10.98, NaN,  NaN), 
          c(17.81, 15.74, 14.06, 12.74, NaN,   NaN,  NaN),  
          c(17.27, 15.56, 14.13, NaN,   NaN,   NaN,  NaN),   
          c(16.72, 15.31, NaN,   NaN,   NaN,   NaN,  NaN),   
          c(15.86, NaN,   NaN,   NaN,   NaN,   NaN,  NaN))  
      }
      if (signif == 4){
        cv<- rbind(
          c(18.29,  14.09,  11.85, 10.15, 8.85,  7.95, 7.00), 
          c(19.32,  15.94,  13.64, 11.92, 10.60, 9.41, NaN),
          c(19.00,  16.40,  14.37, 12.86, 11.50, NaN,  NaN), 
          c(18.77,  16.40,  14.71, 13.36, NaN,   NaN,  NaN),  
          c(18.17,  16.23,  14.69, NaN,   NaN,   NaN,  NaN),   
          c(17.52,  16.02,  NaN,   NaN,   NaN,   NaN,  NaN),   
          c(16.65,  NaN,    NaN,   NaN,   NaN,   NaN,  NaN))  
      }
    }  
    if (q==10){
      if (signif == 1){
        cv<- rbind(
          c(15.20, 11.94,  9.87,  8.56,  7.55,  6.73, 5.97), 
          c(17.19, 14.21,  12.03, 10.60, 9.37,  8.43, NaN),
          c(17.61, 15.11,  13.17, 11.71, 10.51, NaN,  NaN), 
          c(17.41, 15.41,  13.69, 12.41, NaN,   NaN,  NaN),  
          c(17.02, 15.44,  13.93, NaN,   NaN,   NaN,  NaN),   
          c(16.51, 15.09,  NaN,   NaN,   NaN,   NaN,  NaN),   
          c(15.80, NaN,    NaN,   NaN,   NaN,   NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(16.55, 12.83, 10.62, 9.19,  8.11,  7.20, 6.39), 
          c(18.21, 15.07, 12.79, 11.20, 9.90,  8.87, NaN),
          c(18.51, 15.88, 13.89, 12.27, 11.02, NaN,  NaN), 
          c(18.28, 16.18, 14.37, 12.98, NaN,   NaN,  NaN),  
          c(17.76, 16.11, 14.51, NaN,   NaN,   NaN,  NaN),   
          c(17.24, 15.76, NaN,   NaN,   NaN,   NaN,  NaN),   
          c(16.48,   NaN, NaN,   NaN,   NaN,   NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(17.71, 13.69, 11.37, 9.76,  8.60,  7.65, 6.79), 
          c(19.26, 15.85, 13.32, 11.79, 10.41, 9.31, NaN),
          c(19.39, 16.66, 14.52, 12.81, 11.48, NaN,  NaN), 
          c(19.12, 16.90, 15.05, 13.47, NaN,   NaN,  NaN),  
          c(18.48, 16.72, 15.11, NaN,   NaN,   NaN,  NaN),   
          c(17.94, 16.27, NaN,   NaN,   NaN,   NaN,  NaN),   
          c(17.05, NaN,   NaN,   NaN,   NaN,   NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(19.36, 14.74, 12.09, 10.56, 9.26,  8.18, 7.26), 
          c(20.29, 16.76, 14.05, 12.38, 10.97, 9.83, NaN),
          c(20.36, 17.47, 15.17, 13.48, 12.05, NaN,  NaN), 
          c(20.01, 17.66, 15.73, 14.04, NaN,   NaN,  NaN),  
          c(19.25, 17.47, 15.86, NaN,   NaN,   NaN,  NaN),   
          c(18.77, 16.98, NaN,   NaN,   NaN,   NaN,  NaN),   
          c(17.70, NaN,   NaN,   NaN,   NaN,   NaN,  NaN))  
      }
    }  
  }
  
  if (trm == .15){
    cv <- matrix(0,4,4)
    if (q==1){
      if (signif == 1){
        cv<- rbind(
          c(6.21,  5.72,  5.09,  4.39), 
          c(5.75,  5.46,  5.00,  NaN),   
          c(5.06,  4.97,  NaN,   NaN),   
          c(4.38,  NaN,   NaN,   NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(7.18,  6.46,  5.68,  4.86), 
          c(6.49,  6.13,  5.50,  NaN),   
          c(5.69,  5.48,  NaN,   NaN),   
          c(4.91,  NaN,   NaN,   NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(8.12,  7.23,  6.16,  5.33), 
          c(7.17,  6.71,  5.95,  NaN),   
          c(6.17,  5.94,  NaN,   NaN),   
          c(5.34,  NaN,   NaN,   NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(9.24,  8.00,  6.79,  5.91), 
          c(7.98,  7.45,  6.48,  NaN),   
          c(6.76,  6.57,  NaN,   NaN),   
          c(5.90,  NaN,   NaN,   NaN))
      }
    }
    if (q==2){
      if (signif == 1){
        cv<- rbind(
          c(7.45, 6.54, 5.68, 4.91), 
          c(7.31, 6.66, 5.87, NaN),   
          c(6.74, 6.29, NaN,  NaN),   
          c(6.01, NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(8.45, 7.36, 6.27, 5.42), 
          c(8.12, 7.33, 6.45, NaN),   
          c(7.46, 6.86, NaN,  NaN),   
          c(6.53, NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(9.45, 8.02, 6.84, 5.90), 
          c(8.91, 7.88, 6.94, NaN),   
          c(8.15, 7.41, NaN,  NaN),   
          c(7.06, NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(10.73, 8.93, 7.57, 6.49), 
          c(9.90, 8.73, 7.58, NaN),   
          c(8.84, 8.04, NaN,  NaN),   
          c(7.64, NaN,  NaN,  NaN))
      }
    }   
    if (q==3){
      if (signif == 1){
        cv<- rbind(
          c(8.53, 7.30, 6.20, 5.32), 
          c(8.63, 7.63, 6.66, NaN),   
          c(8.12, 7.41, NaN,  NaN),   
          c(7.44, NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(9.52, 8.07, 6.84, 5.88), 
          c(9.51, 8.31, 7.23, NaN),   
          c(8.86, 8.05, NaN,  NaN),   
          c(8.05, NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(10.61, 8.80, 7.39, 6.39), 
          c(10.30, 8.98, 7.74, NaN),   
          c(9.48, 8.63, NaN,  NaN),   
          c(8.61, NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(11.87, 9.67, 8.16, 6.89), 
          c(11.30, 9.80, 8.38, NaN),   
          c(10.25, 9.33, NaN,  NaN),   
          c(9.30, NaN,  NaN,  NaN))
      }
    }   
    if (q==4){
      if (signif == 1){
        cv<- rbind(
          c(9.51, 7.87, 6.67, 5.66), 
          c(9.90, 8.56, 7.38, NaN),   
          c(9.45, 8.40, NaN,  NaN),   
          c(8.86, NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(10.54, 8.73, 7.32, 6.17), 
          c(10.83, 9.30, 8.02, NaN),   
          c(10.27, 9.03, NaN,  NaN),   
          c(9.53, NaN,  NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(11.61, 9.47, 7.92, 6.73), 
          c(11.62, 9.98, 8.60, NaN),   
          c(11.01, 9.58, NaN,  NaN),   
          c(10.10, NaN,  NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(13.08, 10.42, 8.71, 7.26), 
          c(12.62, 10.73, 9.28, NaN),   
          c(11.84, 10.29, NaN,  NaN),   
          c(10.86, NaN,   NaN,  NaN))   
      }
    }   
    if (q==5){
      if (signif == 1){
        cv<- rbind(
          c(10.45, 8.53, 7.18, 6.03), 
          c(11.03, 9.41, 8.10, NaN),   
          c(10.69, 9.43, NaN,  NaN),   
          c(10.05, NaN,  NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(11.66, 9.33,  7.82, 6.56), 
          c(12.01, 10.13, 8.70, NaN),   
          c(11.49, 10.09, NaN,  NaN),   
          c(10.78, NaN,   NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(12.72, 10.09, 8.43, 7.10), 
          c(12.89, 10.82, 9.23, NaN),   
          c(12.28, 10.72, NaN,  NaN),   
          c(11.34, NaN,   NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(14.06, 11.15, 9.12, 7.70), 
          c(14.13, 11.67, 9.87, NaN),   
          c(13.19, 11.38, NaN,  NaN),   
          c(12.07, NaN,   NaN,  NaN))
      }
    }   
    if (q==6){
      if (signif == 1){
        cv<- rbind(
          c(11.27, 9.07,  7.52, 6.38), 
          c(12.08, 10.21, 8.73, NaN),   
          c(11.85, 10.34, NaN,  NaN),   
          c(11.22, NaN,   NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(12.50, 9.90,  8.20, 6.92), 
          c(13.02, 11.01, 9.36, NaN),   
          c(12.72, 11.11, NaN,  NaN),   
          c(12.02, NaN,   NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(13.62, 10.61, 8.81, 7.49), 
          c(13.90, 11.72, 9.98, NaN),   
          c(13.49, 11.76, NaN,  NaN),   
          c(12.62, NaN,   NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(14.90, 11.53, 9.66,  8.21), 
          c(14.99, 12.53, 10.63, NaN),   
          c(14.35, 12.59, NaN,   NaN),   
          c(13.42, NaN,   NaN,   NaN))
      }
    }   
    if (q==7){
      if (signif == 1){
        cv<- rbind(
          c(12.17, 9.69,  7.98, 6.66), 
          c(13.22, 11.03, 9.33, NaN),   
          c(13.01, 11.30, NaN,  NaN),   
          c(12.46, NaN,   NaN,  NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(13.42, 10.55, 8.67, 7.23), 
          c(14.24, 11.76, 9.99, NaN),   
          c(13.93, 11.99, NaN,  NaN),   
          c(12.46, NaN,   NaN,  NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(14.63, 11.28, 9.34,  7.80), 
          c(15.16, 12.48, 10.65, NaN),   
          c(14.76, 12.63, NaN,   NaN),   
          c(13.82, NaN,   NaN,   NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(16.06, 12.10, 10.20, 8.40), 
          c(16.28, 13.33, 11.36, NaN),   
          c(15.74, 13.34, NaN,   NaN),   
          c(14.65, NaN,   NaN,   NaN))
      }
    }   
    if (q==8){
      if (signif == 1){
        cv<- rbind(
          c(13.02, 10.19, 8.36,  7.03), 
          c(14.23, 11.77, 10.01, NaN),   
          c(14.23, 12.19, NaN,   NaN),   
          c(13.53, NaN,   NaN,   NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(14.38, 11.12, 9.04,  7.60), 
          c(15.30, 12.65, 10.71, NaN),   
          c(15.10, 12.93, NaN,   NaN),   
          c(14.34, NaN,   NaN,   NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(15.55, 11.90, 9.73,  8.20), 
          c(16.28, 13.39, 11.32, NaN),   
          c(15.92, 13.67, NaN,   NaN),   
          c(15.03, NaN,   NaN,   NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(17.10, 12.93, 10.67, 8.86), 
          c(17.28, 14.22, 12.06, NaN),   
          c(16.78, 14.48, NaN,   NaN),   
          c(15.92, NaN,   NaN,   NaN))
      }
    }   
    if (q==9){
      if (signif == 1){
        cv<- rbind(
          c(13.81, 10.74, 8.77,  7.35), 
          c(15.27, 12.48, 10.61, NaN),   
          c(15.29, 13.08, NaN,   NaN),   
          c(14.70, NaN,   NaN,   NaN))
      }
      if (signif == 2){
        cv<- rbind( 
          c(15.05, 11.61, 9.53,  7.94), 
          c(16.35, 13.35, 11.30, NaN),   
          c(16.19, 13.83, NaN,   NaN),   
          c(15.48, NaN,   NaN,   NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(16.12, 12.45, 10.23, 8.47), 
          c(17.33, 14.13, 11.91, NaN),   
          c(17.05, 14.63, NaN,   NaN),   
          c(16.14, NaN,   NaN,   NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(17.61, 13.42, 11.04, 9.18), 
          c(18.49, 15.14, 12.75, NaN),   
          c(18.05, 15.44, NaN,   NaN),   
          c(17.06, NaN,   NaN,   NaN))
      }
    }   
    if (q==10){
      if (signif == 1){
        cv<- rbind(
          c(14.62, 11.23, 9.14,  7.64), 
          c(16.26, 13.36, 11.20, NaN),   
          c(16.33, 13.93, NaN,   NaN),   
          c(15.77, NaN,   NaN,   NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(16.03, 12.11,  9.86,  8.23), 
          c(17.35, 14.19,  11.86, NaN),   
          c(17.37, 14.72,  NaN,   NaN),   
          c(16.57,  NaN,   NaN,   NaN))
      }
      if (signif == 3){
        cv<- rbind(
          c(17.15, 12.94, 10.53, 8.78), 
          c(18.37, 14.99, 12.53, NaN),   
          c(18.24, 15.49,  NaN,  NaN),   
          c(17.30,  NaN,   NaN,  NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(18.73, 14.04, 11.40, 9.52), 
          c(19.57, 15.95, 13.26, NaN),   
          c(19.10, 16.38, NaN,   NaN),   
          c(18.23, NaN,   NaN,   NaN))
      }
    }   
  }
  
  if (trm == .20){
    cv <- matrix(0,2,2)
    if (q==1){
      if (signif == 1){
        cv<- rbind(
          c(5.83,5.18), 
          c(5.19,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(6.79,5.89), 
          c(5.93,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(7.70,6.70), 
          c(6.56,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(8.83,7.52), 
          c(7.42,NaN))
      }
    }  
    if (q==2){
      if (signif == 1){
        cv<- rbind(
          c(7.10,6.01), 
          c(6.72,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(8.12,6.77), 
          c(7.52,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(9.08,7.50), 
          c(8.34,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(10.27,8.33), 
          c(9.31,NaN))
      }
    }  
    if (q==3){   
      if (signif == 1){
        cv<- rbind(
          c(8.09,6.70), 
          c(7.94,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(9.11,7.50), 
          c(8.77,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(10.18,8.25), 
          c(9.59,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(11.50,9.09), 
          c(10.50,NaN))
      }
    }  
    if (q==4){   
      if (signif == 1){
        cv<- rbind(
          c(9.09,7.31), 
          c(9.17,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(10.14,8.14), 
          c(10.01,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(11.17,8.91), 
          c(10.89,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(12.67,9.76), 
          c(11.90,NaN))
      }
    }  
    if (q==5){   
      if (signif == 1){
        cv<- rbind(
          c(9.99,7.94), 
          c(10.36,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(11.20,8.75), 
          c(11.33,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(12.28,9.54), 
          c(12.22,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(13.56,10.52), 
          c(13.29,NaN))
      }
    }  
    if (q==6){   
      if (signif == 1){
        cv<- rbind(
          c(10.80,8.50), 
          c(11.32,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(12.02,9.30), 
          c(12.30,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(13.16,10.08), 
          c(13.14,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(14.49,10.94), 
          c(14.19,NaN))
      }
    }  
    if (q==7){   
      if (signif == 1){
        cv<- rbind(
          c(11.70,9.09), 
          c(12.56,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(12.89,9.97), 
          c(13.63,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(14.08,10.68), 
          c(14.57,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(15.46,11.58), 
          c(15.67,NaN))
      }
    }  
    if (q==8){   
      if (signif == 1){
        cv<- rbind(
          c(12.56,9.52), 
          c(13.46,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(13.84,10.46), 
          c(14.55,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(14.87,11.26), 
          c(15.45,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(16.39,12.43), 
          c(16.60,NaN))
      }
    }  
    if (q==9){   
      if (signif == 1){
        cv<- rbind(
          c(13.33,10.15), 
          c(14.37,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(14.59,11.08), 
          c(15.46,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(15.67,11.78), 
          c(16.51,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(17.03,12.85), 
          c(17.69,NaN))
      }
    }  
    if (q==10){   
      if (signif == 1){
        cv<- rbind(
          c(14.09,10.61), 
          c(15.37,NaN))
      }
      if (signif == 2){
        cv<- rbind(
          c(15.53,11.48), 
          c(16.51,NaN))
      }   
      if (signif == 3){
        cv<- rbind(
          c(16.73,12.39), 
          c(17.59,NaN))
      }
      if (signif == 4){
        cv<- rbind(
          c(18.15,13.47), 
          c(18.75,NaN))
      }
    }  
  }
  
  if (trm == .25){
    if (q==1){
      if (signif == 1){
        cv <- 5.48    
      }
      if (signif == 2){
        cv <- 6.43
      }
      if (signif == 3){
        cv <- 7.42
      }
      if (signif == 4){
        cv <- 8.56
      }
    }
    if (q==2){
      if (signif == 1){
        cv <- 6.70
      }
      if (signif == 2){
        cv <- 7.72
      }
      if (signif == 3){
        cv <- 8.69
      }
      if (signif == 4){
        cv <-9.94
      }
    }
    if (q==3){
      if (signif == 1){
        cv <- 7.67
      }
      if (signif == 2){
        cv <- 8.75
      }
      if (signif == 3){
        cv <- 9.73
      }
      if (signif == 4){
        cv <- 10.89
      }
    }
    if (q==4){
      if (signif == 1){
        cv <- 8.66
      }
      if (signif == 2){
        cv <- 9.73
      }
      if (signif == 3){
        cv <- 10.87
      }
      if (signif == 4){
        cv <- 12.33
      }
    }
    if (q==5){
      if (signif == 1){
        cv <- 9.56
      }
      if (signif == 2){
        cv <- 10.73
      }
      if (signif == 3){
        cv <- 11.93
      }
      if (signif == 4){
        cv <-13.23
      }
    }
    if (q==6){
      if (signif == 1){
        cv <- 10.37
      }
      if (signif == 2){
        cv <- 11.56
      }
      if (signif == 3){
        cv <- 12.62
      }
      if (signif == 4){
        cv <- 14.03
      }
    }
    if (q==7){
      if (signif == 1){
        cv <- 11.23
      }
      if (signif == 2){
        cv <- 12.43
      }
      if (signif == 3){
        cv <- 13.50
      }
      if (signif == 4){
        cv <- 14.87
      }
    }
    if (q==8){
      if (signif == 1){
        cv <- 12.07
      }
      if (signif == 2){
        cv<- 13.32
      }
      if (signif == 3){
        cv <- 14.44
      }
      if (signif == 4){
        cv <- 16.04
      }
    }
    if (q==9){
      if (signif == 1){
        cv <- 12.83
      }
      if (signif == 2){
        cv <- 14.11
      }
      if (signif == 3){
        cv <- 15.16
      }
      if (signif == 4){
        cv <- 16.62
      }
    }
    if (q==10){
      if (signif == 1){
        cv <- 13.61
      }
      if (signif == 2){
        cv <- 14.84
      }
      if (signif == 3){
        cv <- 16.28
      }
      if (signif == 4){
        cv <- 17.71
      }
    }
  }
  
  return(cv)
}


#' @title Get critical values for UDmaxLRT_4 
#' 
#' @description
#' This function retrieves critical values for the UDmaxLRT_4 test from the paper by Perron, Yamamoto, and Zhou (2020).
#'
#' @param alpha Numeric: The significance level (e.g., 0.05 for a 5% significance level).
#' @param q Integer: The number of coefficients for the regression model.
#' @param trm Numeric: The threshold multiplier.
#' @param M Integer: The number of coefficient breaks.
#' @param N Integer: The number of variance breaks.
#'
#' @return
#' A numeric value representing the critical value for the UDmaxLRT_4 test.
#'
#' @references
#' Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#'
#' @keywords internal
#' 
#' @export
getdmax4 <- function(alpha, q, trm, M, N){
  
  siglev <- as.matrix(c(0.10,0.05,0.025,0.01))
  signif <- which(siglev==alpha)
  
  if (M==2 & N==2){
    indmn <- 1
  }else if (M==3 & N==2){
    indmn <- 2
  }else if (M==2 & N==3){
    indmn <- 3
  }else if (M==3 & N==3){
    indmn <- 4
  }else{
    disp('The critical values are not available.')
  }
  
  # input of the critical values of the UDmax4 tests
  # For each matrix, the row is q (number of regressors) and the column
  # is number of maximum breaks (no column if unavailable).
  #  1st column: M=2, N=2
  #  2nd column: M=3, N=2
  #  3rd column: M=2, N=3
  #  4th column: M=3, N=3
  
  if (trm == .05){
    if (signif == 1){    
      cv<- rbind(
        c(7.93, 	7.94,  7.94,  7.96),
        c(9.28,   9.35,  9.29,  9.36),
        c(10.57,  10.68, 10.57, 10.68),
        c(11.74,  11.96, 11.74, 11.96),
        c(12.91,  13.22, 12.92, 13.22),
        c(14.01,  14.40, 14.01, 14.40),
        c(15.06,  15.55, 15.06, 15.55),
        c(16.13,  16.69, 16.13, 16.69),
        c(17.24,  17.88, 17.24, 17.88),
        c(18.30,  19.06, 18.30, 19.06))
    }
    if (signif == 2){
      cv<- rbind(
        c(8.69, 	8.72,  8.70,  8.74),
        c(10.13,  10.18, 10.14, 10.18),
        c(11.46,  11.57, 11.46, 11.57),
        c(12.68,  12.81, 12.68, 12.81),
        c(13.86,  14.06, 13.86, 14.06),
        c(14.89,  15.24, 14.89, 15.24),
        c(16.00,  16.42, 16.00, 16.42),
        c(17.13,  17.63, 17.13, 17.63),
        c(18.21,  18.80, 18.21, 18.80),
        c(19.30,  20.01, 19.30, 20.01))
    }
    if (signif == 3){
      cv<- rbind(
        c(9.41,   9.43,  9.43,  9.44),
        c(10.94,  10.95, 10.95, 10.96),
        c(12.25,  12.29, 12.25, 12.29),
        c(13.57,  13.66, 13.57, 13.66),
        c(14.65,  14.84, 14.65, 14.84),
        c(15.74,  16.03, 15.74, 16.03),
        c(16.98,  17.30, 16.98, 17.30),
        c(18.17,  18.52, 18.17, 18.52),
        c(19.12,  19.65, 19.12, 19.65),
        c(20.25,  20.88, 20.25, 20.88))
    }
    if (signif == 4){
      cv<- rbind(
        c(10.55, 	10.55, 10.56, 10.56),
        c(12.08,  12.10, 12.08, 12.10),
        c(13.39,  13.40, 13.39, 13.40),
        c(14.42,  14.58, 14.42, 14.58),
        c(15.67,  15.84, 15.67, 15.84),
        c(16.88,  17.04, 16.88, 17.04),
        c(18.42,  18.65, 18.42, 18.65),
        c(19.48,  19.63, 19.48, 19.63),
        c(20.31,  20.71, 20.31, 20.71),
        c(21.44,  21.83, 21.44, 21.83))
    }
  }
  
  if (trm == .10){
    if (signif == 1){
      cv<- rbind(
        c(7.18, 	7.18,  7.18,  7.18),
        c(8.47,	  8.51,  8.48,  8.52),
        c(9.73,	  9.85,  9.74,  9.85),
        c(10.88,	10.99, 10.88, 10.99),
        c(12.07,	12.28, 12.07, 12.28),
        c(13.11,  13.37, 13.11, 13.37),
        c(14.15,  14.50, 14.15, 14.50),
        c(15.20,  15.63, 15.20, 15.63),
        c(16.15,  16.63, 16.15, 16.63),
        c(17.28,  17.84, 17.28, 17.84))
    }
    if (signif == 2){
      cv<- rbind(
        c(8.03,	8.04,  8.04,  8.05),
        c(9.37,	9.37,  9.37,  9.37),
        c(10.66,	10.71, 10.66, 10.71),
        c(11.85,  11.94, 11.81, 11.94),
        c(13.06,	13.19, 13.06, 13.19),
        c(14.14,	14.32, 14.14, 14.32),
        c(15.12,  15.40, 15.12, 15.40),
        c(16.22,  16.53, 16.22, 16.53),
        c(17.26,  17.59, 17.26, 17.59),
        c(18.38,  18.80, 18.38, 18.80))
    }
    if (signif == 3){
      cv<- rbind(
        c(8.81,	  8.94,  8.94,  8.94),
        c(10.32,	10.32, 10.32, 10.32),
        c(11.48,	11.58, 11.54, 11.58),
        c(12.81,	12.82, 12.81, 12.82),
        c(13.99,	14.08, 13.99, 14.08),
        c(15.00,  15.15, 15.00, 15.15),
        c(16.12,  16.33, 16.12, 16.33),
        c(17.15,  17.46, 17.15, 17.46),
        c(18.25,  18.59, 18.25, 18.59),
        c(19.27,  19.77, 19.27, 19.77))
    }
    if (signif == 4){
      cv<- rbind(
        c(10.00,	10.00, 10.00, 10.00),
        c(11.47,	11.47, 11.47, 11.47),
        c(12.66,	12.72, 12.69, 12.72),
        c(13.99,	13.99, 13.99, 13.99),
        c(15.16,	15.16, 15.16, 15.16),
        c(16.05,  16.22, 16.05, 16.22),
        c(17.37,  17.51, 17.37, 17.51),
        c(18.37,  18.67, 18.37, 18.67),
        c(19.52,  19.74, 19.52, 19.74),
        c(20.45,  20.74, 20.45, 20.74))
    }
  }
  
  if (trm == .15){
    if (signif == 1){
      cv<- rbind(
        c(6.61,   6.61,  6.62),
        c(7.93,   7.93,  7.93),
        c(9.09,   9.15,  9.09),
        c(10.24,  10.33, 10.25),
        c(11.33,  11.49, 11.35),
        c(12.30,  12.49, 12.30),
        c(13.37,  13.56, 13.37),
        c(14.35,  14.63, 14.35),
        c(15.29,  15.65, 15.29),
        c(16.34,  16.78, 16.34))
    }
    if (signif == 2){
      cv<- rbind(
        c(7.51,  7.51,  7.51),
        c(8.88,  8.88,  8.88),
        c(10.08, 10.11, 10.08),
        c(11.19, 11.28, 11.21),
        c(12.38, 12.51, 12.45),
        c(13.34, 13.47, 13.34),
        c(14.40, 14.62, 14.40),
        c(15.49, 15.68, 15.49),
        c(16.35, 16.62, 16.35),
        c(17.51, 17.81, 17.51))
    }
    if (signif == 3){
      cv<- rbind(
        c(8.32,  8.32,  8.32),
        c(9.77,  9.77,  9.77),
        c(10.93, 10.99, 10.95),
        c(12.20, 12.20, 12.20),
        c(13.38, 13.38, 13.38),
        c(14.36, 14.47, 14.36),
        c(15.38, 15.53, 15.38),
        c(16.51, 16.66, 16.51),
        c(17.32, 17.56, 17.32),
        c(18.55, 18.78, 18.55))
    }
    if (signif == 4){
      cv<- rbind(
        c(9.42,  9.42,  9.42),
        c(10.96, 10.96, 10.96),
        c(12.19, 12.19, 12.19),
        c(13.39, 13.43, 13.41),
        c(14.50, 14.50, 14.50),
        c(15.53, 15.55, 15.53),
        c(16.67, 16.78, 16.67),
        c(17.82, 17.91, 17.82),
        c(18.77, 18.90, 18.77),
        c(19.85, 20.00, 19.85))
    }
  }
  
  
  if (trm == .20){
    if (signif == 1){
      cv<- rbind(
        c(6.15 ),
        c(7.39 ),
        c(8.55 ),
        c(9.64 ),
        c(10.70),
        c(11.73),
        c(12.72),
        c(13.71),
        c(14.61),
        c(15.60))
    }
    if (signif == 2){
      cv<- rbind(
        c(7.05),
        c(8.42),
        c(9.48),
        c(10.66),
        c(11.84),
        c(12.79),
        c(13.82),
        c(14.92),
        c(15.76),
        c(16.75))
    }
    if (signif == 3){
      cv<- rbind(
        c(7.87),
        c(9.40),
        c(10.41),
        c(11.53),
        c(12.86),
        c(13.80),
        c(14.76),
        c(15.99),
        c(16.90),
        c(17.78))
    }
    if (signif == 4){
      cv<- rbind(
        c(8.95),
        c(10.54),
        c(11.64),
        c(12.84),
        c(13.95),
        c(15.02),
        c(16.13),
        c(17.23),
        c(18.27),
        c(19.26))
    }
  }
  
  if (ncol(cv)<indmn){
    cvmn <- NaN
  }else{
    cvmn <- cv[q,indmn]
  }
  
  return(cvmn)
}
