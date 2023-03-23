
#' @title Determine number of breaks & dates
#' 
#' @description This function determines the number of break and their dates using sequential and UDmax procedures. 
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
determineBreaks <- function(y, M, N, z, x, con){
  # start with UDmax test up to M breaks
  control_tmp9 <- list(robust = con$robust,
                       prewhit = con$prewhit,
                       typek = con$typekbc,
                       kerntype = con$kerntype,
                       alpha = con$alpha)
  control_tmp10 <- list(vrobust = con$vrobust,
                        prewhit = con$prewhit,
                        typek = con$typekbv,
                        kerntype = con$kerntype,
                        alpha = con$alpha)
  control_tmpm  <- list(robust = con$robust,
                       prewhit = con$prewhit,
                       typekc = con$typekbc,
                       kerntype = con$kerntype,
                       alpha = con$alpha)
  control_tmpn  <- list(vrobust = con$vrobust,
                       prewhit = con$prewhit,
                       typekbv = con$typekbv,
                       kerntype = con$kerntype,
                       alpha = con$alpha)
  control_tmpmn   <- list(robust = con$robust,
                          vrobust = con$vrobust,
                          prewhit = con$prewhit,
                          typekbv = con$typekbv,
                          typekbc = con$typekbc,
                          kerntype = con$kerntype,
                          alpha = con$alpha)
  
  if ((M>0) & (N==0)){
    out <- pslr00(y, M, con$trm, z, x, control_tmpm)
    # check if null hypothesis is rejected
    if (is.null(con$alpha)){
      UDmax_null_check <- out$UDmaxLRT>out$cvUDmax[1,1]  
    }else{
      UDmax_null_check <- out$UDmaxLRT>=out$cvUDmax[1,((100-as.numeric(gsub("%","",colnames(out$cvUDmax))))/100)==con$alpha]  
    }
    # if rejected, use Sequential procedure to find m starting with m=1 under null.
    if (UDmax_null_check){
      mi <- 1
      seq_null_check <- TRUE
      while ((mi < M) & (seq_null_check)){
        out <- pslr9(y, mi, N, con$trm, z, x, control_tmp9)
        if (is.null(con$alpha)){
          seq_null_check <- out$supSeq>=out$cv[1,1]
        }else{
          seq_null_check <- out$supSeq>=out$cv[1,((100-as.numeric(gsub("%","",colnames(out$cv))))/100)==con$alpha]  
        }
        if (seq_null_check){
          mi <- mi + 1  
        }
      }
      out <- pslr0(y, mi, con$trm, z, x, control_tmpm)
      brcdt <- out$brcstar # global break date
      brvdt <- NULL
      m <- mi
      n <- N
    }else{
      stop("No breaks in coefficients found.")
    }
  }else if ((M==0) & (N>0)){
    # start with UDmax test up to M breaks
    out <- pslr5(y, N, con$trm, z, x, control_tmpn)
    # check if null hypothesis is rejected
    if (is.null(con$alpha)){
      UDmax_null_check <- out$UDmaxLRT>out$cvUDmax[1,1]  
    }else{
      UDmax_null_check <- out$UDmaxLRT>=out$cvUDmax[1,((100-as.numeric(gsub("%","",colnames(out$cvUDmax))))/100)==con$alpha]  
    }
    # if rejected, use Sequential procedure to find m starting with m=1 under null.
    if (UDmax_null_check){
      ni <- 1
      seq_null_check <- TRUE
      while ((ni < N) & (seq_null_check)){
        out <- pslr10(y, M, ni, con$trm, z, x, control_tmp10)
        if (is.null(con$alpha)){
          seq_null_check <- out$supSeq>=out$cv[1,1]
        }else{
          seq_null_check <- out$supSeq>=out$cv[1,((100-as.numeric(gsub("%","",colnames(out$cv))))/100)==con$alpha]  
        }
        if (seq_null_check){
          ni <- ni + 1  
        }
      }
      out <- pslr1(y, ni, con$trm, z, x, control_tmpn)
      brcdt <- NULL
      brvdt <- out$brvstar # global break date
      m <- M
      n <- ni
    }else{
      stop("No breaks in variance found.")
    }
  }else if ((M>0) & (N>0)){
    mi <- 0
    seq_null_check <- TRUE
    while ((mi < M) & (seq_null_check)){
      out <- pslr9(y, mi, N, con$trm, z, x, control_tmp9)
      if (is.null(con$alpha)){
        seq_null_check <- out$supSeq>=out$cv[1,1]
      }else{
        seq_null_check <- out$supSeq>=out$cv[1,((100-as.numeric(gsub("%","",colnames(out$cv))))/100)==con$alpha]  
      }
      if (seq_null_check){
        mi <- mi + 1  
      }
    }
    ni <- 0
    seq_null_check <- TRUE
    while ((ni < N) & (seq_null_check)){
      out <- pslr10(y, M, ni, con$trm, z, x, control_tmp10)
      if (is.null(con$alpha)){
        seq_null_check <- out$supSeq>=out$cv[1,1]
      }else{
        seq_null_check <- out$supSeq>=out$cv[1,((100-as.numeric(gsub("%","",colnames(out$cv))))/100)==con$alpha]  
      }
      if (seq_null_check){
        ni <- ni + 1  
      }
    }
    m <- mi 
    n <- ni
    if ((m>0) & (n>0)){
      out <- pslr4(y, m, n, con$trm, z, x, control_tmpmn)  
      brcdt <- out$brcstar # global break date
      brvdt <- out$brvstar # global break date
    }else if ((m>0) & (n==0)){
      out <- pslr0(y, m, con$trm, z, x, control_tmpm)  
      brcdt <- out$brcstar # global break date
      brvdt <- NULL
    }else if ((m==0) & (n>0)){
      out <- pslr1(y, n, con$trm, z, x, control_tmpn)  
      brcdt <- NULL
      brvdt <- out$brvstar # global break date
    }
  }
  return(list(m = m, n = n, brcdt = brcdt, brvdt = brvdt, supLRT = out$suplr, cv = out$cv))
}





#' @title Estimate model 
#' 
#' @description Model estimation with known or unknown breaks
#' 
#' @param y A (\code{T x 1}) vector with endogenous variable.
#' @param z A (\code{T x q}) matrix with explanatory variables subject to change.
#' @param m An integer determining the number of breaks in coefficients.
#' @param n An integer determining the number of breaks in variance.
#' @param control List with options which include
#' \itemize{
#'   \item{\code{robust}: }{set to 1 if want to allow for heterogeneity and autocorrelation in the residuals, 0 otherwise. The method used is Andrews(1991) automatic bandwidth with AR(1) approximation and the quadratic kernel. Note: Do not set to 1 if lagged dependent variables are included as regressors. Default is \code{TRUE}.}
#'   \item{\code{prewhit}: }{set to 1 if want to apply AR(1) prewhitening prior to estimating the long run covariance matrix. Default is \code{TRUE}.}
#'   \item{\code{hetdat}: }{option for the construction of the F tests. Set to 1 if want to allow different moment matrices of the regressors across segments. If hetdat=0, the same moment matrices are assumed for each segment and estimated from the ful sample. It is recommended to set hetdat=1 if p>0. Default is \code{TRUE}.}
#'   \item{\code{hetvar}: }{option for the construction of the F tests.Set to 1 if want to allow for the variance of the residuals to be different across segments. If hetvar=0, the variance of the residuals is assumed constant across segments and constructed from the ful sample. This option is not available when robust=1. Default is \code{TRUE}.}
#'   \item{\code{hetomega}: }{used in the construction of the confidence intervals for the break dates. If hetomega=0, the long run covariance matrix of zu is assumed identical across segments (the variance of the errors u if robust=0). Default is \code{TRUE}.}
#'   \item{\code{hetq}: }{used in the construction of the confidence intervals for the break dates. If hetq=0, the moment matrix of the data is assumed identical across segments. Default is \code{TRUE}.}
#'   \item{\code{doglobal}: }{set to 1 if want to cal the procedure to obtain global minimizers. Default is \code{TRUE}.}
#'   \item{\code{hetomega}: }{used in the construction of the confidence intervals for the break dates. If hetomega=0, the long run covariance matrix of zu is assumed identical across segments (the variance of the errors u if robust=0Default is \code{TRUE}.}
#'   \item{\code{hetq}: }{used in the construction of the confidence intervals for the break dates. If hetq=0, the moment matrix of the data is assumed identical across segments. Default is \code{TRUE}.}
#'   \item{\code{}: }{Default is \code{TRUE}.}
#' }
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
estimdl <- function(y, m, n, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(robust = TRUE,
              vrobust = TRUE,
              prewhit = FALSE,
              typekbv = 2,
              typekbc = 2,
              brcdt = NULL,
              brvdt = NULL,
              trm = 0.1,
              hetdat    = TRUE,
              hetvar    = TRUE,
              hetomega  = TRUE,
              hetq      = TRUE,
              kerntype = "quadratic",
              alpha = NULL)
  # ----- Perform some checks for controls
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])){
    warning("unknown names in control: ", paste(noNms,collapse=", ")) 
  }
  # ----- Perform checks for inputs 
  if (is.matrix(z)){
    q <- ncol(z)
  }else{
    stop("z must be a (T x q) matrix.") 
  }
  if (length(x)!=0){
    if (is.matrix(x)){
      p <- ncol(x)
    }else{
      stop("x must be a (T x p) matrix.") 
    }  
  }else{
    p <-0
  }
  if (is.matrix(y)){
    bigt <- nrow(y)
  }else{
    stop("y must be a (T x 1) matrix.") 
  }
  # ----- Obtain h from trm
  h <- round(con$trm*bigt)
  mdlout <- list()
  mdlout$m <- m
  mdlout$n <- n
  # Find break dates if not given
  if ((is.null(con$brcdt)) & (is.null(con$brvdt))){
    brkout <- determineBreaks(y, m, n, z, x, con)
    con$brcdt <- brkout$brcdt
    con$brvdt <- brkout$brvdt
    mdlout$m <- brkout$m
    mdlout$n <- brkout$n
    mdlout$supLRT <- brkout$supLRT
    mdlout$cv <- brkout$cv
  }
  # ----- esimate model
  if ((is.null(con$brcdt)==FALSE) & (is.null(con$brvdt))){
    estim_out <- estim(y,z,x,mdlout$m,con$brcdt,con$robust,con$prewhit,con$hetomega,con$hetq,con$hetdat,con$hetvar)
    beta <- estim_out$coef
    stdev <- sqrt(estim_out$sigma2)
    brc <- con$brcdt
    brv <- matrix(0,0,0)
    res <-  estim_out$resid
  }else{
    # step 1
    if ((is.null(con$brcdt)) & (is.null(con$brvdt)==FALSE)){
      brc <- matrix(0,0,0)
      brk <- as.matrix(sort(unique(c(brcdt,con$brvdt))))
      K <- length(brk)
      cbrind <- as.matrix(as.numeric(brk %in% con$brcdt))
      vbrind <- as.matrix(as.numeric(brk %in% con$brvdt))
    }else if ((is.null(con$brcdt)==FALSE) & (is.null(con$brvdt)==FALSE)){
      brk <- as.matrix(sort(unique(c(con$brcdt,con$brvdt))))
      K <- length(brk)
      cbrind <- as.matrix(as.numeric(brk %in% con$brcdt))
      vbrind <- as.matrix(as.numeric(brk %in% con$brvdt))
    }
    # step 2
    segmake_out <- segmake(K,brk,mdlout$m,mdlout$n,cbrind,vbrind,q)
    estimbr_out <- estimbr(y,z,q,x,p,bigt,K,brk,segmake_out$R,mdlout$n,segmake_out$brv,1)
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
        segmake_out <- segmake(K,brk,mdlout$m,mdlout$n,cbrind,vbrind,q)
        estimbr_out <- estimbr(y,z,q,x,p,bigt,K,brk,segmake_out$R,mdlout$n,segmake_out$brv,1)
        diff        <- diff + 1
      }
      if (diff==maxiter){
        print('cannot find the converging break dates') 
      }
    }
    vseg  <- as.matrix(c(0,segmake_out$brv,bigt))
    nvar <- matrix(0,mdlout$n+1,1)
    for (k in 1:(mdlout$n+1)){
      i <- vseg[k,]+1
      j <- vseg[(k+1),]
      nvar[k,]  <- sqrt((t(estimbr_out$res[i:j,])%*%estimbr_out$res[i:j,])/(j-i+1))
    }
    beta <- estimbr_out$nbeta
    stdev <- nvar
    brc <- segmake_out$brc
    brv <- segmake_out$brv
    res <-  estimbr_out$res
  }
  # ----- Organize output
  mdlout$y = y 
  mdlout$z = z 
  mdlout$x = x
  mdlout$beta = beta
  mdlout$stdev = stdev
  mdlout$brk = brk
  mdlout$brc = brc
  mdlout$brv = brv
  mdlout$res = res
  mdlout$control = con
  # *** NOTE: CI are computed for m>0 and n==0 but left out for now after other 
  #           cases do not have these ready yet. 
  class(mdlout) <- "mdl"
  return(mdlout)
}