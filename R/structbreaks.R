
#' @title ols parameter vector estimation
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @export
olsqr <- function(y, x){
  b <- solve(t(x)%*%x)%*%t(x)%*%y 
  return(as.matrix(b))
}

#' @title Estimate OLS with White robust standard errors
#' 
#' @description This function estimates an ordinary least squares model with White robust 
#' standard error.
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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

#' @title Compute quadratic kernel
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @export
kern <- function(x){
  # procedure to evaluate the quadratic kernel at some value x.
  del <- 6*pi*x/5
  ker <- 3*(sin(del)/del-cos(del))/(del*del)
  return(ker)
}


#' @title Barlette Kernel
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



#' @title Diagonal partition of observations by regimes
#' 
#' @description procedure to construct the diagonal partition of z with m break at date b.
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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

#' @title Create matrix if regime dates
#' 
#' @export
regim_dt <- function(bvec){
  mr <- length(bvec) - 1
  regime_dates <- matrix(0, mr, 2)
  for (xr in 1:mr){
    regime_dates[xr,1] <- bvec[xr]
    regime_dates[xr,2] <- bvec[xr+1]
  }
  return(regime_dates)
}


#' @title Covariance matrix of estimates delta.
#' 
#' @description procedure that compute the covariance matrix of the estimates delta.
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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

#' @title diagonal matrix with variance for regime i
#' 
#' @description procedure that computes a diagonal matrix of dimension i+1 with ith entry
#' the estimate of the variance of the residuals for segment i.  
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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

#' @title diagonal matrix with break fractions
#' 
#' @description procedure that construct a diagonal matrix of dimension m+1 with ith
#' entry (T_i-T_i-1)/T.  
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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


#' @title compute long run variance
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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

#' @title bandwidth based on AR(1) approximation
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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
      if (cct > 100){
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




#' @title Confidence intervals for breaks
#' 
#' @description This procedure that computes confidence intervals 
#' for the break dates based on the "shrinking shifts" asymptotic framework.
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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
#' @description This function estimates model after breaks are found. It also 
#' computes and reports confidence intervals for the break dates. The method used 
#' depends on the specification for robust
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
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
#' @references Bai, J. and Perron, P. (2003), Critical values for multiple structural change tests. \emph{The Econometrics Journal}, 6: 72-78. https://doi.org/10.1111/1368-423X.00102
#' @export
getcv <- function(alpha, q, k, eps = 0.1){
  if (alpha == 0.05){
    cL  <- (8.238+1.756*q-0.043*q^2-0.659*k-15.436*eps+0.025*q/eps)*exp(0.389/k-0.013/(eps*k))/q
    cUD <- (8.228+3.095*q-9.644*eps)*exp(-0.029*q-0.291*eps)/q
    cWD <- (9.039+3.318*q-9.969*eps)*exp(-0.030*q-0.327*eps)/q  
  }else if (alpha == 0.1){
    cL  <- (7.551+1.718*q-0.041*q^2-0.610*k-15.846*eps+0.025*q/eps)*exp(0.338/k-0.014/(eps*k))/q
    cUD <- (6.917+2.930*q-9.275*eps)*exp(-0.028*q-0.406*eps)/q
    cWD <- (7.316+3.128*q-8.624*eps)*exp(-0.029*q-0.412*eps)/q  
  }
  output <- list(cvL = cL, cvUD = cUD, cvWD = cWD)
  return(output)
}


#' @title Determine if square matrix is positive definite
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
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



#' @title Inverse of positive definite matrix 
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
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



#' @title Compute robust standard errors for correct1() function
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
jhatpr1 <- function(vmat,vmata){
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



#' @title Compute robust standard errors in joint test of sc (i.e., Perron et al. 2021)
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
correct1 <- function(vmat, vmata, prewhit, typekb){
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
    jh            <- jhatpr1(vstar, vstara)
    hac           <- invpd(diag(d)-bmat)$xinv%*%jh%*%t(inv(eye(d)-bmat))  
  }else{
    hac           <- jhatpr1(vmat,vmata)
  }
  return(hac)
}


#' @title Find segments 
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
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
#' @description procedure to estimate coefficients and residuals by FGLS given break dates in coefficients and variance.
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
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





#' @title Log likelihood given residuakls
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
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






#' @title check Determine number of cases to check given number of breaks in mean or var
#' 
#' @description replaces brcvcase() and numcase() in original MATLAB code.
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
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


#' @title Get bigvec of residuals
#' 
#' @description replaces residuals() in original MATLAB code
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
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


#' @title Compute global break dates for pure structural change model
#' 
#' @description This is the main procedure which calculates the break points that globally
#'  minimizes the SSR. It returns optimal dates and associated SSR for all numbers of breaks less than or equal to m.
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @param y A (\code{T x 1}) vector with endogenous variable.
#' @param z A (\code{T x q}) matrix with explanatory variables subject to change.
#' @param m An integer determining the number of breaks to find.
#' @param h An integer determining the minimum length of a regime. 
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @export
dating_purescSSR <- function(y, z, m, h){
  # ----- Set some values
  bigT      <- nrow(y)
  q         <- ncol(z)
  datevec   <- matrix(0,m,m)
  optdat    <- matrix(0,bigT,m)
  optssr    <- matrix(0,bigT,m)
  glb       <- matrix(0,m,1)
  dvec      <- matrix(0,bigT,1)
  bigvec    <- ssrbigvec(y, z, h)
  # Determine global optimum break dates
  if (m==1){
    ssrmin_datx <- parti(1, h, bigT-h, bigT, bigvec, bigT)
    datevec[1,1] <- ssrmin_datx$dx
    glb[1,1] <- ssrmin_datx$ssrmin
  }else{
    for (j1 in (2*h):bigT){
      ssrmin_datx <- parti(1, h, j1-h, j1, bigvec, bigT)
      optdat[j1,1] <- ssrmin_datx$dx
      optssr[j1,1] <- ssrmin_datx$ssrmin
    }
    glb[1,1] <- optssr[bigT,1]
    datevec[1,1] <- optdat[bigT,1]
    
    for (ib in 2:m){
      if (ib==m){
        jlast <- bigT
        for (jb in (ib*h):(jlast-h)){
          dvec[jb,1] <- optssr[jb,(ib-1)] + bigvec[(jb+1)*bigT-jb*(jb+1)/2,1]
        }
        optssr[jlast,ib] <- t(min(dvec[(ib*h):(jlast-h),1]))
        minindcdvec <-  which.min(dvec[(ib*h):(jlast-h),1])
        optdat[jlast,ib] <- (ib*h-1) + t(minindcdvec)
      }else{
        for (jlast in ((ib+1)*h):bigT){
          for (jb  in (ib*h):(jlast-h)){
            dvec[jb,1] <- optssr[jb,(ib-1)] + bigvec[jb*bigT-jb*(jb-1)/2+jlast-jb,1]
          }
          optssr[jlast,ib] <-  min(dvec[(ib*h):(jlast-h),1])
          minindcdvec <- which.min(dvec[(ib*h):(jlast-h),1])
          optdat[jlast,ib] <- (ib*h-1) + t(minindcdvec)
        }
      }
      datevec[ib,ib] <- optdat[bigT,ib]
      for (i in 1:(ib-1)){
        xx  <- ib-i
        datevec[xx,ib] <- optdat[datevec[(xx+1),ib],xx]  
      }
      glb[ib,1] <- optssr[bigT,ib]
    }
  }
  return(list(glob = glb, datevec = datevec, bigvec = bigvec))
}



#' @title Compute global break dates in partial structural break model
#' 
#' @details Note: This code is an adaptation of the one originally written by Yohei 
#' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
#' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
#' 
#' @export
dating_partscSSR <- function(y, z, x, m, h, thtol = 1e-6, maxi = 10000){
  # ----- Set some values
  bigT      <- nrow(y)
  q         <- ncol(z)
  p         <- ncol(x)
  glb       <- matrix(0, m, 1)
  datevec   <- matrix(0, m, m)
  #mi        <- 1
  for (mi in 1:m){
    qq <- p+q
    zz <- cbind(z, c(x))
    puresc_out <- dating_purescSSR(y, zz, mi, h, con)
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
    if ((length > thtol) & (i > con$maxi)){
      msg     <- 'The number of iterations has reached the upper limit'
    }else{
      msg   <- 'converged'
    }
  }
  return(list(glob = glb, datevec = datevec, bigvec = bigvec, convergence_msg = msg))
}

#' @title Compute global break dates using log-likelihood 
#' 
#' @description This is the main procedure which calculates the break points that globally maximize the loglikelihood function. It returns optimal dates and associated log likelihood for all numbers of breaks less than or equal to m. 
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
dating_loglik <- function(bigvec,h,m,bigt){
  datevec   <- matrix(0,m,m)
  optdat    <- matrix(0,bigt,m) 
  optlr     <- matrix(0,bigt,m)                    
  dvec      <- matrix(0,bigt,1)
  glob      <- matrix(0,m,1)
  
  if (m==1){
    dating_out <- parti_loglik(1,h,bigt-h,bigt,bigvec,bigt)
    datevec[1,1]    <- dating_out$dx
    glob[1,1]       <- dating_out$lrmax
  }else{
    for (j1 in (2*h):bigt){
      dating_out    <-  parti_loglik(1,h,j1-h,j1,as.matrix(bigvec[1:(2*bigt),]),bigt)    
      optlr[j1,1]   <- dating_out$lrmax
      optdat[j1,1]  <- dating_out$dx
    }                                                                      
    glob[1,1]       <- optlr[bigt,1]
    datevec[1,1]    <- optdat[bigt,1]
    for (ib in 2:m){
      if (ib==m){
        jlast <- bigt
        for (jb in (ib*h):(jlast-h)){
          dvec[jb,1] <- optlr[jb,(ib-1)]-0.5*(bigt-jb+1)*((log(2*pi)+1)+log(sum(bigvec[(m*bigt+jb+1):(bigt*(m+1)),])/(bigt-jb)))
        }
        optlr[jlast,ib]  <- max(dvec[(ib*h):(jlast-h),1])
        optdat[jlast,ib] <- (ib*h-1) + which.max(dvec[(ib*h):(jlast-h),1])
      }else{
        for (jlast in ((ib+1)*h):bigt){
          for (jb in (ib*h):(jlast-h)){
            dvec[jb,1] <- optlr[jb,(ib-1)]-0.5*(jlast-jb+1)*((log(2*pi)+1)+log(sum(bigvec[(ib*bigt+jb+1):(ib*bigt+jlast),])/(jlast-jb)))  
          }
          optlr[jlast,ib] <- max(dvec[(ib*h):(jlast-h),1])
          optdat[jlast,ib] <- (ib*h-1)+which.max(dvec[(ib*h):(jlast-h),1])
        }
      }
      datevec[ib,ib] <- optdat[bigt,ib]
      for (i in 1:(ib-1)){
        xx <- ib-i
        datevec[xx,ib] <- optdat[datevec[(xx+1),ib],xx]
      }
      glob[ib,1] <- optlr[bigt,ib]
    }
    
  }
  return(list(glob = glob, datevec = datevec))
}



#' @title Compute global break dates using log-likelihood 
#' 
#' @description This is the main procedure which calculates the break points that globally maximize the loglikelihood function. It returns optimal dates and associated log likelihood for all numbers of breaks less than or equal to m. 
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
dating_MLE <- function(y,z,q,x,p,h,m,bigt){
  datevec <- matrix(0, m, m) 
  optdat  <- matrix(0, bigt, m)
  optmle  <- matrix(0, bigt, m)      
  dvec    <- matrix(0, bigt, 1)                        
  glob    <- matrix(0, m, 1)
  bigvec  <- mlebigvec(y, z, q, x, p, h, bigt)
  # ----- up to here
  if (m == 1){
    parti_out <- parti(1,h,bigt-h,bigt,bigvec,bigt)
    datevec[1,1]  <- parti_out$dx 
    glob[1,1]     <- parti_out$ssrmin
  }else{
    for (j1 in (2*h):bigt){
      parti_out <- parti(1,h,j1-h,j1,bigvec,bigt)           
      optmle[j1,1]  <- parti_out$ssrmin
      optdat[j1,1]  <- parti_out$dx 
    }
    glob[1,1]      <- optmle[bigt,1]
    datevec[1,1]   <- optdat[bigt,1]
    for (ib in 2:m){
      if (ib == m){
        jlast <- bigt
        for (jb in (ib*h):(jlast-h)){
          dvec[jb,1]  <- optmle[jb,ib-1] + bigvec[(jb+1)*bigt-jb*(jb+1)/2,1]
        }
        optmle[jlast,ib]      <- min(dvec[(ib*h):(jlast-h),1])
        minindcdvec           <- which.min(dvec[(ib*h):(jlast-h),1])
        optdat[jlast,ib]      <- (ib*h-1) + t(minindcdvec)
      }else{
        for (jlast in ((ib+1)*h):bigt){
          for (jb in (ib*h):(jlast-h)){
            dvec[jb,1]   <- optmle[jb,ib-1] + bigvec[jb*bigt-jb*(jb-1)/2+jlast-jb,1] 
          }
          optmle[jlast,ib] <- min(dvec[(ib*h):(jlast-h),1])
          minindcdvec  <- which.min(dvec[(ib*h):(jlast-h),1])
          optdat[jlast,ib] <- (ib*h-1) + t(minindcdvec)
        }
      }
      datevec[ib,ib] <- optdat[bigt,ib]
      for (i in 1:(ib-1)){
        xx <- ib-i
        datevec[xx,ib] <- optdat[datevec[xx+1,ib],xx]  
      }
      glob[ib,1] <- optmle[bigt,ib]
    }
  }
  glob <- -glob
  return(list(glob = glob,datevec = datevec,bigvec = bigvec))
}




#' @title Joint estimation of mean and variance coefficients
#' 
#' @description procedure to jointly estimates the coefficient and variance break dates following section 5.1 of Qu and Perron (2007). The code follows est.m of Qu and Perron (2007)
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
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

