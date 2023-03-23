


#' @title Test stat for  0 vs m breaks in mean given n=0 breaks in var
#' 
#' @description Computes the supLRT test statistic for m coefficient changes given no variance changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr0 <- function(y, m, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(robust = TRUE,
              prewhit = FALSE,
              typekc = 2,
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
  h <- round(trm*bigt)
  # ----- 1) log-likelihood function under the null
  if (p==0){
    reg0 <- z
  }else{
    reg0 <- cbind(z, x)
  }
  res0  <-  y - reg0%*%invpd(t(reg0)%*%reg0)$xinv%*%t(reg0)%*%y
  vvar0 <- t(res0)%*%res0/bigt
  lr0   <- -(bigt/2)*(log(2*pi)+1+log(vvar0))
  # ----- 2) log-likelihood function under the alternative
  if (p==0){
    datevec <- dating_purescSSR(y, z, m, h)$datevec 
  }else{
    datevec <- dating_partscSSR(y, z, x, m, h)$datevec 
  }
  brc <- as.matrix(datevec[,m])
  zbar <- pzbar(z,m,brc)
  if (p==0){
    reg1 <- zbar
  }else{
    reg1 <- cbind(zbar, x)
  }
  beta1 <- invpd(t(reg1)%*%reg1)$xinv%*%t(reg1)%*%y
  res1  <- y - reg1%*%beta1
  vvar1 <- t(res1)%*%res1/bigt
  lr1   <- -(bigt/2)*(log(2*pi)+1+log(vvar1))
  
  if (con$robust==FALSE){
    suplr <- 2*(lr1-lr0)
  }else if(con$robust==TRUE){
    if (p==0){
      vmat0 <- matrix(0, bigt, q)
      vmat1 <- matrix(0, bigt, q)
      for (i in 1:q){
        vmat0[, i] <- reg0[, i]*res0
        vmat1[, i] <- reg0[, i]*res1
      }
      hac              <- correct1(vmat0, vmat1, con$prewhit, con$typekc,con$kerntype)
      lambda           <- plambda(brc, m, bigt)
      vdel             <- bigt*invpd(t(reg1)%*%reg1)$xinv%*%kronecker(lambda,hac)%*%invpd(t(reg1)%*%reg1)$xinv
      delta1           <- beta1
    }else{
      regm             <- zbar - x%*%invpd(t(x)%*%x)$xinv%*%t(x)%*%zbar
      vmat0            <- matrix(0, bigt, q*(m+1))
      vmat1            <- matrix(0, bigt, q*(m+1))
      for (i in 1:(q*(m+1))){
        vmat0[,i] <- regm[,i]*res0
        vmat1[,i] <- regm[,i]*res1
      }
      hac              <- correct1(vmat0, vmat1, con$prewhit, con$typekc,con$kerntype)
      vdel             <- bigt*invpd(t(regm)%*%regm)$xinv%*%hac%*%invpd(t(regm)%*%regm)$xinv
      delta1           <- beta1[1:((m+1)*q),1]
    }
    rsub <- matrix(0, m, m+1)
    for (j in 1:m){
      rsub[j,j] <- -1
      rsub[j,j+1] <- 1 
    }
    rmat <-  kronecker(rsub,diag(q))
    
    fstar <- t(delta1)%*%t(rmat)%*%invpd(rmat%*%vdel%*%t(rmat))$xinv%*%rmat%*%delta1
    suplr <- (bigt - (m+1)*q-p)*fstar/bigt   
  }
  
  brcstar <- brc
  suplr <- as.matrix(suplr/m)
  
  colnames(suplr) <- "supLRT_0"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,q,m,trm)$cvL)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,q,m,trm)$cvL,getcv(0.05,q,m,trm)$cvL, 
                        getcv(0.025,q,m,trm)$cvL, getcv(0.01,q,m,trm)$cvL)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr0_out <- list(suplr = suplr, brcstar = brcstar, cv = cv)
  class(pslr0_out) <- "test"
  return(pslr0_out)
}






#' @title Test stat for  0 vs n breaks in variance given m=0 breaks in mean
#' 
#' @description Computes the sup LR_1T test statistic for n variance changes given no coefficient changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @return suplr - corrected test statistic (see eq. 9 of PYZ 2020)
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr1 <- function(y, n, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(vrobust = TRUE,
              prewhit = FALSE,
              typekbv = 2,
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
  h <- round(trm*bigt)
  # ----- 1) log-likelihood function under the null
  if (p==0){
    reg0 <- z
  }else{
    reg0 <- cbind(z, x)
  }
  res0  <-  y - reg0%*%invpd(t(reg0)%*%reg0)$xinv%*%t(reg0)%*%y
  vvar0 <- t(res0)%*%res0/bigt
  tao0  <- ((res0^2)/c(vvar0))-1
  lr0   <- -(bigt/2)*(log(2*pi)+1+log(vvar0))
  # ----- 2) log-likelihood function under the alternative
  bigvec <- matrix(0, bigt*(n+1), 1)
  for (i in 1:(n+1)){
    bigvec[((i-1)*bigt+1):(i*bigt),] <- res0^2
  }
  datevec <- dating_loglik(bigvec,h,n,bigt)$datevec
  brv <- as.matrix(datevec[,n])
  seg_out   <- segmake(n,brv,0,n,matrix(0,n,1),matrix(1,n,1),q)
  estim_out <- estimbr(y,z,q,x,p,bigt,n,seg_out$brv,seg_out$R,n,seg_out$brv,1)
  plog_out  <- ploglik(estim_out$res, n, seg_out$brv)
  if (con$vrobust==0){
    phi <- (t(tao0)%*%tao0)/(bigt-1)
  }else if (con$vrobust==1){
    phi <- correct1(tao0,plog_out$tao,con$prewhit,con$typekbv,con$kerntype) 
  }
  lr1 <- plog_out$loglik
  suplr   <- (2/phi)*2*(lr1-lr0)
  brvstar <- brv
  suplr   <- as.matrix(suplr/n)
  
  colnames(suplr) <- "supLRT_1"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,1,n,trm)$cvL)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,1,n,trm)$cvL,getcv(0.05,1,n,trm)$cvL, 
                        getcv(0.025,1,n,trm)$cvL, getcv(0.01,1,n,trm)$cvL)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr1_out <- list(suplr = suplr, brvstar = brvstar, cv = cv)
  class(pslr1_out) <- "test"
  return(pslr1_out)
}


#' @title Test stat for  0 vs n breaks in variance given M=m breaks in mean
#' 
#' @description Computes the sup LR_2T test statistic for n variance changes given m coefficient changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @return suplr - corrected test statistic (see eq. 9 of PYZ 2020)
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr2 <- function(y, m, n, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(vrobust = TRUE,
              prewhit = FALSE,
              typekbv = 2,
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
  h <- round(trm*bigt)
  # ----- 1) log-likelihood function under the null
  if (p==0){
    datevec <- dating_purescSSR(y, z, m, h)$datevec 
    brc     <- as.matrix(datevec[,m])
    zbar    <- pzbar(z, m, brc)
    reg0    <- zbar
  }else{
    datevec <- dating_partscSSR(y, z, x, m, h)$datevec 
    brc     <- as.matrix(datevec[,m])
    zbar    <- pzbar(z, m, brc)
    reg0    <- cbind(zbar, x)
  }
  res0  <-  y - reg0%*%invpd(t(reg0)%*%reg0)$xinv%*%t(reg0)%*%y
  vvar0 <- t(res0)%*%res0/bigt
  tao0  <- ((res0^2)/c(vvar0))-1
  lr0   <- -(bigt/2)*(log(2*pi)+1+log(vvar0))
  
  # ----- 2) log-likelihood function under the alternative
  brcase_tmp          <- numcase2(m,n)
  suplrx              <- matrix(0,brcase_tmp$num,1)
  brcdt               <- matrix(0,brcase_tmp$num,m)
  brvdt               <- matrix(0,brcase_tmp$num,n)
  
  for (idx in 1:brcase_tmp$num){
    brcase_out <- brcase_tmp$cvbrind[[idx]]
    K <- ncol(brcase_out)
    cbrind <- as.matrix(brcase_out[1,])
    vbrind <- as.matrix(brcase_out[2,])
    estdate_out <- estdate(y,z,q,x,p,K,bigt,h,m,n,cbrind,vbrind)
    plog_out  <- ploglik(estdate_out$res,n,estdate_out$brv)
    brcdt[idx,] <- t(brc)
    brvdt[idx,] <- t(estdate_out$brv)
    
    if (con$vrobust==0){
      phi <- t(tao0)%*%tao0/(bigt-1)
    }else if (con$vrobust==1){
      phi <- correct1(tao0, plog_out$tao, con$prewhit, con$typekbv,con$kerntype)
    }
    lr1 <- plog_out$loglik
    suplrx[idx,1] <- (2/phi)*2*(lr1-lr0)
  }
  suplr   <- max(suplrx)
  maxind  <- which.max(suplrx)
  brcstar <- t(as.matrix(brcdt[maxind,]))
  brvstar <- t(as.matrix(brvdt[maxind,]))
  suplr   <- as.matrix(suplr/n)
  
  colnames(suplr) <- "supLRT_2"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,1,n,trm)$cvL)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,1,n,trm)$cvL,getcv(0.05,1,n,trm)$cvL, 
                        getcv(0.025,1,n,trm)$cvL, getcv(0.01,1,n,trm)$cvL)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr2_out <- list(suplr = suplr, brcstar = brcstar, brvstar = brvstar, cv = cv)
  class(pslr2_out) <- "test"
  return(pslr2_out)
}



#' @title Test stat for  0 vs m breaks in mean given N=n breaks in var
#' 
#' @description Computes the sup LR_3T test statistic for m coefficient changes given n variance changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr3 <- function(y, m, n, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(robust = TRUE,
              prewhit = FALSE,
              typek = 2,
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
  if (n == 0){
    # Use pslr0() i.e. Bai & Perron (1998)
    controlc <- list(robust = con$robust,
                     prewhit = con$prewhit,
                     typekc = con$typek) 
    out <- pslr0(y, m, trm, z, x, controlc)
    suplr <- as.matrix(out$suplr)
    brcstar <- out$brcstar
    brvstar <- matrix(0,n,1)
  }else if (n>0){
    # ----- Obtain h from trm
    h <- round(trm*bigt)
    # ----- 1)log-likelihood function under the null
    if (p==0){
      reg0 <- z  
    }else{
      reg0 <- cbind(z,x)
    }
    res0    <- y - reg0%*%invpd(t(reg0)%*%reg0)$xinv%*%t(reg0)%*%y
    bigvec  <- matrix(0, bigt*(n+1),1)
    for (i in 1:(n+1)){
      bigvec[((i-1)*bigt+1):(i*bigt),] <- res0^2
    }
    datevec <- dating_loglik(bigvec,h,n,bigt)$datevec
    brv     <- as.matrix(datevec[,n])
    segmake_out <- segmake(n, brv,0, n, matrix(0,n,1), matrix(1,n,1), ncol(z))
    brv   <- segmake_out$brv
    R     <- segmake_out$R
    res0  <-  estimbr(y,z,q,x,p,bigt,n,brv,R,n,brv,1)$res
    lr0   <- ploglik(res0, n, brv)$loglik
    
    # ----- 2) log-likelihood function under the alternative
    brcase_tmp          <- numcase2(m,n)
    suplrx              <- matrix(0,brcase_tmp$num,1)
    brcdt               <- matrix(0,brcase_tmp$num,m)
    brvdt               <- matrix(0,brcase_tmp$num,n)
    
    for (idx in 1:brcase_tmp$num){
      brcase_out <- brcase_tmp$cvbrind[[idx]]
      K <- ncol(brcase_out)
      cbrind <- as.matrix(brcase_out[1,])
      vbrind <- as.matrix(brcase_out[2,])
      estdate_out <- estdate(y,z,q,x,p,K,bigt,h,m,n,cbrind,vbrind)
      brc <-estdate_out$brc
      brv <- estdate_out$brv
      res1 <- estdate_out$res
      brcdt[idx,]     <- t(brc)
      brvdt[idx,]     <- t(brv)
      lr1             <- ploglik(res1,n,brv)$loglik   
      suplrx[idx,1]   <- 2*(lr1-lr0)
    }
    suplr   <- max(suplrx)
    maxind  <- which.max(suplrx)
    brcstar <- as.matrix(brcdt[maxind,])
    brvstar <- as.matrix(brvdt[maxind,])
    brcase_out_opt <- brcase_tmp$cvbrind[[maxind]]
    K_opt <- ncol(brcase_out_opt)
    cbrind_opt <- as.matrix(brcase_out_opt[1,])
    vbrind_opt <- as.matrix(brcase_out_opt[2,])
    estdate_out_opt <- estdate(y,z,q,x,p,K_opt,bigt,h,m,n,cbrind_opt,vbrind_opt)
    res1 <- estdate_out_opt$res
    if (con$robust==1){
      zbar <- pzbar(z, m, brcstar)
      if (p==0){
        reg1 <- zbar  
      }else{
        reg1 <- cbind(zbar, x)
      }
      ibigv <- diag(bigt)
      vseg  <- as.matrix(c(0,brvstar,bigt))
      for (k in 1:(n+1)){
        i <- vseg[k,] + 1
        j <- vseg[k+1,]
        vvar <- t(as.matrix(res1[i:j,]))%*%as.matrix(res1[i:j,])/(j-i+1) 
        ibigv[i:j,i:j]  <- diag(j-i+1)*c(invpd(vvar)$xinv)
      }
      
      ys    <- ibigv^(1/2)%*%y
      zbars <- ibigv^(1/2)%*%zbar
      reg0s <- ibigv^(1/2)%*%reg0
      reg1s <- ibigv^(1/2)%*%reg1
      
      beta1s  <- invpd(t(reg1s)%*%reg1s)$xinv%*%t(reg1s)%*%ys
      res0s   <- y - reg0%*%invpd(t(reg0s)%*%reg0s)$xinv%*%t(reg0s)%*%ys
      res1s   <- y - reg1%*%invpd(t(reg1s)%*%reg1s)$xinv%*%t(reg1s)%*%ys
      if (p==0){
        vmat0 <- matrix(0,bigt,q)
        vmat1 <- matrix(0,bigt,q)
        for (i in 1:q){
          vmat0[,i] <- reg0[,i]*res0s
          vmat1[,i] <- reg0[,i]*res1s  
        }
        hac     <- correct1(vmat0,vmat1,con$prewhit,con$typek,con$kerntype)
        lambda  <- plambda(brcstar,m,bigt)
        vdel    <- bigt*invpd(t(reg1)%*%reg1)$xinv%*%kronecker(lambda,hac)%*%invpd(t(reg1)%*%reg1)$xinv
        delta1s <- beta1s
      }else{
        xs    <- ibigv^(1/2)%*%x 
        regms <- zbars - xs%*%invpd(t(xs)%*%xs)$xinv%*%t(xs)%*%zbars
        vmat0 <- matrix(0, bigt, q*(m+1))
        vmat1 <- matrix(0, bigt, q*(m+1))
        for (i in 1:(q*(m+1))){
          vmat0[,i] <- regms[,i]*res0s 
          vmat1[,i] <- regms[,i]*res1s
        }
        hac     <- correct1(vmat0,vmat1,con$prewhit,con$typek,con$kerntype)
        vdel    <- bigt*invpd(t(regms)%*%regms)$xinv%*%hac%*%invpd(t(regms)%*%regms)$xinv
        delta1s <- as.matrix(beta1s[1:((m+1)*q),1])
      }
      rsub <- matrix(0,m,m+1)
      for (j in 1:m){
        rsub[j,j]     <- -1
        rsub[j,j+1]   <- 1
      }
      rmat <- kronecker(rsub,diag(q))
      
      fstar <- t(delta1s)%*%t(rmat)%*%invpd(rmat%*%vdel%*%t(rmat))$xinv%*%rmat%*%delta1s
      suplr <- (bigt-(m+1)*q-p)*fstar/bigt 
    }
    suplr <- as.matrix(suplr/m)
  }
  
  colnames(suplr) <- "supLRT_3"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,q,m,trm)$cvL)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,q,m,trm)$cvL,getcv(0.05,q,m,trm)$cvL, 
                        getcv(0.025,q,m,trm)$cvL, getcv(0.01,q,m,trm)$cvL)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr3_out <- list(suplr = suplr, brcstar = brcstar, brvstar = brvstar, cv = cv)
  class(pslr3_out) <- "test"
  return(pslr3_out)
}



#' @title Test stat for  0 vs m breaks in mean and 0 vs n breaks in var
#' 
#' @description Computes the sup LR_4T test statistic for m coefficient changes and n variance changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr4 <- function(y, m, n, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(robust = TRUE,
              vrobust = TRUE,
              prewhit = FALSE,
              typekbv = 2,
              typekbc = 2,
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
  h <- round(trm*bigt)
  # ----- 1) log-likelihood function under the null
  if (p==0){
    reg0 <- z
  }else{
    reg0 <- cbind(z, x)
  }
  res0  <- y - reg0%*%invpd(t(reg0)%*%reg0)$xinv%*%t(reg0)%*%y
  vvar0 <- t(res0)%*%res0/bigt
  tao0  <- ((res0^2)/c(vvar0))-1
  lr0   <- -(bigt/2)*(log(2*pi)+1+log(vvar0))
  
  # ----- 2)log-likelihood function under the alternative
  brcase_tmp          <- numcase2(m,n)
  suplrx              <- matrix(0,brcase_tmp$num,1)
  suplrv              <- matrix(0,brcase_tmp$num,1)
  brcdt               <- matrix(0,brcase_tmp$num,m)
  brvdt               <- matrix(0,brcase_tmp$num,n)
  for (idx in 1:brcase_tmp$num){
    brcase_out    <- brcase_tmp$cvbrind[[idx]]
    K             <- ncol(brcase_out)
    cbrind        <- as.matrix(brcase_out[1,])
    vbrind        <- as.matrix(brcase_out[2,])
    estdate_out   <- estdate(y,z,q,x,p,K,bigt,h,m,n,cbrind,vbrind)
    brc           <- estdate_out$brc
    brv           <- estdate_out$brv
    res1          <- estdate_out$res
    plog_out      <- ploglik(res1,n,brv)
    lr1           <- plog_out$loglik
    tao1          <- plog_out$tao
    
    # variance part
    zbar          <- pzbar(z, m, brc)
    if (p==0){
      reg1        <- zbar  
    }else{
      reg1        <- cbind(zbar, x)
    }
    resv          <- y - reg1%*%invpd(t(reg1)%*%reg1)$xinv%*%t(reg1)%*%y
    vvarv         <- (t(resv)%*%resv)/bigt
    lrv0          <- -(bigt/2)*(log(2*pi)+1+log(vvarv))
    lrv1          <- ploglik(resv, n, brv)$loglik
    if (con$vrobust==0){
      phi         <- (t(tao0)%*%tao0)/(bigt-1)
    }else if (con$vrobust==1){
      phi         <- correct1(tao0,tao1,con$prewhit,con$typekbv,con$kerntype)
    }
    
    suplrx[idx,1] <- 2*(lr1-lr0)-((phi-2)/phi)*2*(lrv1-lrv0)
    suplrv[idx,1] <- (2/phi)*2*(lrv1-lrv0)
    brcdt[idx,]   <- t(brc)
    brvdt[idx,]   <- t(brv)
  }
  suplr     <- max(suplrx)
  maxind    <- which.max(suplrx)
  suplrvar  <- suplrv[maxind,1]
  brcstar   <- as.matrix(brcdt[maxind,])
  brvstar   <- as.matrix(brvdt[maxind,])
  brcase_out_opt <- brcase_tmp$cvbrind[[maxind]]
  K_opt <- ncol(brcase_out_opt)
  cbrind_opt <- as.matrix(brcase_out_opt[1,])
  vbrind_opt <- as.matrix(brcase_out_opt[2,])
  estdate_out_opt <- estdate(y,z,q,x,p,K_opt,bigt,h,m,n,cbrind_opt,vbrind_opt)
  res1 <- estdate_out_opt$res
  if (con$robust==1){
    zbar    <- pzbar(z, m, brcstar)  
    if (p==0){
      reg1  <- zbar 
    }else{
      reg1  <- cbind(zbar, x)
    }
    ibigv <- diag(bigt)
    vseg  <- as.matrix(c(0,brvstar,bigt))
    for (k in 1:(n+1)){
      i     <- vseg[k,]+1
      j     <- vseg[k+1,]
      vvar  <- t(as.matrix(res1[i:j,]))%*%as.matrix(res1[i:j,])/(j-i+1)
      ibigv[i:j,i:j] <- diag(j-i+1)*c(invpd(vvar)$xinv)
    }
    ys      <- ibigv^(1/2)%*%y
    zbars   <- ibigv^(1/2)%*%zbar
    reg0s   <- ibigv^(1/2)%*%reg0
    reg1s   <- ibigv^(1/2)%*%reg1
    
    beta1s  <- invpd(t(reg1s)%*%reg1s)$xinv%*%t(reg1s)%*%ys
    res0s   <- ys - reg0s%*%invpd(t(reg0s)%*%reg0s)$xinv%*%t(reg0s)%*%ys
    res1s   <- ys - reg1s%*%invpd(t(reg1s)%*%reg1s)$xinv%*%t(reg1s)%*%ys
    
    if (p==0){
      vmat0   <- matrix(0,bigt,q)
      vmat1   <- matrix(0,bigt,q)
      for (i in 1:q){
        vmat0[,i] <- reg0s[,i]*res0s
        vmat1[,i] <- reg0s[,i]*res1s
      }
      hac     <- correct1(vmat0,vmat1,con$prewhit,con$typekbc,con$kerntype)
      lambda  <- plambda(brcstar,m,bigt)
      vdel    <- bigt*invpd(t(reg1s)%*%reg1s)$xinv%*%kronecker(lambda,hac)%*%invpd(t(reg1s)%*%reg1s)$xinv
      delta1s <- beta1s
    }else{
      xs      <- ibigv^(1/2)%*%x 
      regms   <- zbars - xs%*%invpd(t(xs)%*%xs)$xinv%*%t(xs)%*%zbars
      vmat0   <- matrix(0,bigt,q*(m+1))
      vmat1   <- matrix(0,bigt,q*(m+1))
      for (i in 1:(q*(m+1))){
        vmat0[,i] <- regms[,i]*res0s
        vmat1[,i] <- regms[,i]*res1s
      }
      hac     <- correct1(vmat0,vmat1,con$prewhit,con$typekbc,con$kerntype)
      vdel    <- bigt*invpd(t(regms)%*%regms)$xinv%*%hac%*%invpd(t(regms)%*%regms)$xinv
      delta1s <- beta1s[1:((m+1)*q),1]
    }
    
    rsub <- matrix(0,m,m+1)
    for (j in 1:m){
      rsub[j,j]     <- -1
      rsub[j,j+1]   <- 1
    }
    rmat <- kronecker(rsub,diag(q))
    
    fstar     <- t(delta1s)%*%t(rmat)%*%invpd(rmat%*%vdel%*%t(rmat))$xinv%*%rmat%*%delta1s
    suplrcoef <- (bigt-(m+1)*q-p)*fstar/bigt
    suplr     <- suplrcoef + suplrvar
    
  }
  suplr <- as.matrix(suplr/(n+m))
  
  colnames(suplr) <- "supLRT_4"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv4(con$alpha,q,trm)[m,n])
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv4(0.1,q,trm)[m,n],getcv4(0.05,q,trm)[m,n], 
                        getcv4(0.025,q,trm)[m,n], getcv4(0.01,q,trm)[m,n])))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr4_out <- list(suplr = suplr, brcstar = brcstar, brvstar = brvstar, cv = cv)
  class(pslr4_out) <- "test"
  return(pslr4_out)
}



#' @title Test stat for  0 vs M breaks in mean given n=0 breaks in var
#' 
#' @description Computes the UDmaxLRT test statistic for M coefficient changes given no variance changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr00 <- function(y, M, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(robust = TRUE,
              prewhit = FALSE,
              typekc = 2,
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
  
  slr0 <- matrix(0,M,1)
  slr0_ls <- list()
  for (m in 1:M){
    slr0_ls[[m]] <- pslr0(y, m, trm, z, x , con)
    slr0[m,1] <- slr0_ls[[m]]$suplr
  }
  # compute UDmaxLRT_0 test stat
  UDmaxLRT_0 <- as.matrix(max(slr0))
  
  colnames(UDmaxLRT_0) <- "UDmaxLRT_0"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,q,1,trm)$cvUD)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,q,1,trm)$cvUD,getcv(0.05,q,1,trm)$cvUD, 
                        getcv(0.025,q,1,trm)$cvUD, getcv(0.01,q,1,trm)$cvUD)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr00_out <- list(UDmaxLRT = UDmaxLRT_0, testtrace = slr0_ls, cvUDmax = cv)
  class(pslr00_out) <- "test"
  return(pslr00_out)
}



#' @title Test stat for  0 vs N breaks in variance given m=0 breaks in mean
#' 
#' @description Computes the UDmaxLRT_1 test statistic for n variance changes given no coefficient changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr5 <- function(y, N, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(vrobust = TRUE,
              prewhit = FALSE,
              typekbv = 2,
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
  
  
  slr1 <- matrix(0,N,1)
  slr1_ls <- list()
  for (n in 1:N){
    slr1_ls[[n]] <- pslr1(y, n, trm, z, x, con)
    slr1[n,1] <- slr1_ls[[n]]$suplr
  }
  # compute UDmaxLRT_1 test stat
  UDmaxLRT_1 <- as.matrix(max(slr1))
  
  colnames(UDmaxLRT_1) <- "UDmaxLRT_1"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,1,1,trm)$cvUD)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,1,1,trm)$cvUD,getcv(0.05,1,1,trm)$cvUD, 
                        getcv(0.025,1,1,trm)$cvUD, getcv(0.01,1,1,trm)$cvUD)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr5_out <- list(UDmaxLRT = UDmaxLRT_1, testtrace = slr1_ls, cvUDmax = cv)
  class(pslr5_out) <- "test"
  return(pslr5_out)
}



#' @title Test stat for  0 vs N breaks in variance given M=m breaks in mean
#' 
#' @description Computes the UDmaxLR_2T test statistic for n variance changes given m coefficient changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr6 <- function(y, m, N, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(vrobust = TRUE,
              prewhit = FALSE,
              typekbv = 2,
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
  
  
  if (m>0){
    slr2 <- matrix(0,N,1)
    slr2_ls <- list()
    for (n in 1:N){
      slr2_ls[[n]] <- pslr2(y, m, n, trm, z, x, con)
      slr2[n,1] <- slr2_ls[[n]]$suplr
    }
    # compute UDmaxLRT_2 test stat
    UDmaxLRT_2 <- as.matrix(max(slr2))
    
    colnames(UDmaxLRT_2) <- "UDmaxLRT_2"
    if (is.null(con$alpha)==FALSE){
      cv <- as.matrix(getcv(con$alpha,1,1,trm)$cvUD)
      colnames(cv) <- paste0((1-con$alpha)*100,"%")
    }else{
      cv <- t(as.matrix(c(getcv(0.1,1,1,trm)$cvUD,getcv(0.05,1,1,trm)$cvUD, 
                          getcv(0.025,1,1,trm)$cvUD, getcv(0.01,1,1,trm)$cvUD)))
      colnames(cv) <- c("90%","95%","97.5%","99%")
    }
    # ----- Organize output
    pslr6_out <- list(UDmaxLRT = UDmaxLRT_2, testtrace = slr2_ls, cvUDmax = cv)
    class(pslr6_out) <- "test"
  }else if(m==0){
    stop("For m=0, use pslr5() test.")
  }
  return(pslr6_out)
}  
  
#' @title Test stat for  0 vs M breaks in mean given N=n breaks in var
#' 
#' @description Computes the UDmaxLR_3T test statistic for m coefficient changes given n variance changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr7 <- function(y, M, n, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(robust = TRUE,
              prewhit = FALSE,
              typek = 2,
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
  
  slr3 <- matrix(0,M,1)
  slr3_ls <- list()
  for (m in 1:M){
    slr3_ls[[m]] <- pslr3(y, m, n, trm, z, x , con)
    slr3[m,1] <- slr3_ls[[m]]$suplr
  }
  # compute UDmaxLRT_3 test stat
  UDmaxLRT_3 <- as.matrix(max(slr3))
  
  colnames(UDmaxLRT_3) <- "UDmaxLRT_3"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,q,1,trm)$cvUD)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,q,1,trm)$cvUD,getcv(0.05,q,1,trm)$cvUD, 
                        getcv(0.025,q,1,trm)$cvUD, getcv(0.01,q,1,trm)$cvUD)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr7_out <- list(UDmaxLRT = UDmaxLRT_3, testtrace = slr3_ls, cvUDmax = cv)
  class(pslr7_out) <- "test"
  return(pslr7_out)
}


#' @title Test stat for  0 vs M breaks in mean and 0 vs N breaks in var
#' 
#' @description Computes the UDmaxLR_4T test statistic for max of M coefficient changes and max of N variance changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr8 <- function(y, M, N, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(robust = TRUE,
              vrobust = TRUE,
              prewhit = FALSE,
              typekbv = 2,
              typekbc = 2, 
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
  
  slr4 <- matrix(0,M,N)
  slr4_ls <- list()
  count <- 1
  for (m in 1:M){
    for (n in 1:N){
      slr4_ls[[count]] <- pslr4(y, m, n, trm, z, x, con)
      slr4[m,n] <- slr4_ls[[count]]$suplr
      count <- count + 1
    }  
  }
  # compute UDmaxLRT_4 test stat
  UDmaxLRT_4 <- as.matrix(max(slr4))
  
  colnames(UDmaxLRT_4) <- "UDmaxLRT_4"
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getdmax4(con$alpha, q, trm, M, N))
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getdmax4(0.1, q, trm, M, N),getdmax4(0.05, q, trm, M, N), 
                        getdmax4(0.025, q, trm, M, N), getdmax4(0.01, q, trm, M, N))))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr8_out <- list(UDmaxLRT = UDmaxLRT_4, testtrace = slr4_ls, cvUDmax = cv)
  class(pslr8_out) <- "test"
  return(pslr8_out)
}

#' @title Test stat for  m vs m+1 breaks in mean given n breaks in var
#' 
#' @description Computes the SeqLR_9T test statistic for m coefficient changes versus m + 1 coefficient changes given n variance change
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr9 <- function(y, m, n, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(robust = TRUE,
              prewhit = FALSE,
              typek = 2,
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
  h <- round(trm*bigt)
  # ----- 1) break date estimation under the null
  if (p==0){
    reg0 <- z
  }else{
    reg0 <-cbind(z, x)
  }
  res0  <- y - reg0%*%invpd(t(reg0)%*%reg0)$xinv%*%t(reg0)%*%y
  vvar0 <- t(res0)%*%res0/bigt
  tao0  <- (res0^2/c(vvar0))-1
  lr0   <- -(bigt/2)*(log(2*pi)+1+log(vvar0))
  
  brcase_tmp          <- numcase2(m,n)
  suplrx              <- matrix(0,brcase_tmp$num,1)
  brctemp             <- matrix(0,brcase_tmp$num,m)
  brvtemp             <- matrix(0,brcase_tmp$num,n)
  
  for (idx in 1:brcase_tmp$num){
    brcase_out    <- brcase_tmp$cvbrind[[idx]]
    K             <- ncol(brcase_out)
    cbrind        <- as.matrix(brcase_out[1,])
    vbrind        <- as.matrix(brcase_out[2,])
    estdate_out   <- estdate(y,z,q,x,p,K,bigt,h,m,n,cbrind,vbrind)
    brc           <- estdate_out$brc
    brv           <- estdate_out$brv
    res1          <- estdate_out$res
    brctemp[idx,] <- t(brc)
    brvtemp[idx,] <- t(brv)
    plog_out      <- ploglik(res1,n,brv)
    lr1           <- plog_out$loglik
    tao1          <- plog_out$tao
    lrv1          <- ploglik(res0,n,brv)$loglik
    if (con$robust==0){
      phi <- t(tao0)%*%tao0/(bigt-1)
    }else if (con$robust==1){
      phi <- correct1(tao0, tao1, con$prewhit, con$typek,con$kerntype)
    }
    suplrx[idx,1] <- 2*(lr1-lr0)-((phi-2)/phi)*2*(lrv1-lr0)
  }
  maxind          <- which.max(suplrx)
  brc0            <- as.matrix(brctemp[maxind,])
  brv0            <- as.matrix(brvtemp[maxind,])
  if (nrow(brc0)==0){
    zbar  <- z;  
  }else{
    zbar  <- pzbar(z,nrow(brc0),brc0)
  }
  if (p==0){
    reg1  <- zbar
  }else{
    reg1  <- cbind(zbar, x)
  }
  res1 <- y - reg1%*%invpd(t(reg1)%*%reg1)$xinv%*%t(reg1)%*%y
  
  # ----- 2) construct segment-wise LR tests
  nseg          <- m+1
  lrtest        <- matrix(0,nseg,1)
  dv            <- matrix(0,nseg+1,1)
  if (length(brc0)!=0){
    dv[2:nseg,1]  <- brc0  
  }
  dv[nseg+1,1]  <- bigt
  ds            <- matrix(0,nseg,1)
  
  vseg          <- as.matrix(c(0,brv0,bigt))
  ys            <- matrix(0,nrow(y),ncol(y))
  zs            <- matrix(0,nrow(z),ncol(z))
  if (p>=1){
    xs          <- matrix(0,nrow(x),ncol(x))
  }
  
  for (k in 1:(n+1)){
    i         <- vseg[k,]+1
    j         <- vseg[k+1,]
    vvar      <- t(res1[i:j,1])%*%res1[i:j,1]/(j-i+1)
    ys[i:j,]  <- y[i:j,]/c(sqrt(vvar))
    zs[i:j,]  <- z[i:j,]/c(sqrt(vvar))
    if (p>=1){
      xs[i:j,]  <- x[i:j,]/c(sqrt(vvar))
    }
  }
  for (is in 1:nseg){
    lengthi <- dv[(is+1),1] - dv[is,1]
    if (lengthi>=2*h){
      starti    <- dv[is,1]+1
      endi      <- dv[is+1,1]
      segy      <- as.matrix(y[starti:endi,1])
      segz      <- as.matrix(z[starti:endi,])
      segys     <- as.matrix(ys[starti:endi,1])
      segzs     <- as.matrix(zs[starti:endi,])
      
      if (p==0){
        segx      <- matrix(0,0,0)
        segxs     <- matrix(0,0,0)
        segreg    <- segz
        segregs   <- segzs  
      }else{
        segx      <- as.matrix(x[starti:endi,])
        segxs     <- as.matrix(xs[starti:endi,])
        segreg    <- cbind(segz, segx)
        segregs   <- cbind(segzs, segxs)
      }          
      
      indbrv      <- (brv0>starti)*(brv0<endi)
      ni          <- sum(indbrv)
      brvi        <- matrix(0,0,0)
      if (n>=1){
        for (j in 1:n){
          if (indbrv[j,1]==1){
            brvi    <- as.matrix(c(brvi,brv0[j,1]-dv[is,1]))
          }
        }
      }
      segres0     <- segy - segreg%*%invpd(t(segregs)%*%segregs)$xinv%*%t(segregs)%*%segys
      seglr0      <- ploglik(segres0,ni,brvi)$loglik
      
      if (p==0){
        dateveci <- dating_purescSSR(segys,segzs, 1, h)$datevec 
      }else{
        dateveci <- dating_partscSSR(segys, segzs, segxs, 1, h)$datevec 
      }
      
      brci      <- as.matrix(dateveci[1,1])
      segres1   <- estimbr(segy,segz,q,segx,p,lengthi,1,brci,matrix(0,0,0),ni,brvi,0)$res
      
      if (con$robust==0){
        seglr1          <- ploglik(segres1,ni,brvi)$loglik
        lrtest[is,1]    <- 2*(seglr1-seglr0)
        ds[is,1]        <- dv[is,1]+brci
      }else if (con$robust==1){
        segzbar        <- pzbar(segz,1,brci)          
        segzbars       <- pzbar(segzs,1,brci)
        if (p==0){
          segreg1      <- segzbar
          segreg0      <- segz
          segreg1s     <- segzbars
          segreg0s     <- segzs
        }else{
          segreg1      <- cbind(segzbar, segx)
          segreg0      <- cbind(segz, segx)       
          segreg1s     <- cbind(segzbars, segxs)
          segreg0s     <- cbind(segzs, segxs)
        }
        beta1s         <- invpd(t(segreg1s)%*%segreg1s)$xinv%*%t(segreg1s)%*%segys
        segres0s       <- segy - segreg0%*%invpd(t(segreg0s)%*%segreg0s)$xinv%*%t(segreg0s)%*%segys
        segres1s       <- segy - segreg1%*%invpd(t(segreg1s)%*%segreg1s)$xinv%*%t(segreg1s)%*%segys
        
        if (p==0){
          vmat0        <- matrix(0,lengthi,q)
          vmat1        <- matrix(0,lengthi,q)
          for (i in 1:q){
            vmat0[,i] <- segreg0[,i]*segres0s
            vmat1[,i] <- segreg0[,i]*segres1s 
          }
          hac         <- correct1(vmat0,vmat1,con$prewhit,con$typek,con$kerntype)
          lambda      <- plambda(brci,1,lengthi)
          vdel        <- lengthi*invpd(t(segreg1)%*%segreg1)$xinv%*%kronecker(lambda,hac)%*%invpd(t(segreg1)%*%segreg1)$xinv
          delta1s     <- beta1s
        }else{
          segregms    <- segzbars - ((segxs%*%invpd(t(segxs)%*%segxs)$xinv%*%t(segxs)%*%segys)%*%matrix(1,1,ncol(segzbars)))
          vmat0       <- matrix(0,lengthi,q*2)
          vmat1       <- matrix(0,lengthi,q*2)
          for (i in 1:(q*2)){
            vmat0[,i] <- segregms[,i]*segres0s
            vmat1[,i] <- segregms[,i]*segres1s
          }
          hac         <- correct1(vmat0,vmat1,con$prewhit,con$typek,con$kerntype)
          vdel        <- lengthi*invpd(t(segregms)%*%segregms)$xinv%*%hac%*%invpd(t(segregms)%*%segregms)$xinv
          delta1s     <- as.matrix(beta1s[1:(q*2),1])
        }
        rsub           <- matrix(0,1,2)
        rsub[1,1]      <- -1
        rsub[1,2]      <- 1
        rmat           <- kronecker(rsub,diag(q))
        
        fstari         <- t(delta1s)%*%t(rmat)%*%invpd(rmat%*%vdel%*%t(rmat))$xinv%*%rmat%*%delta1s
        lrtest[is,1]   <- (lengthi-2*q-p)*fstari/lengthi
        # ds[is,1]    This should be updated... missing only when robust==1...
      }
    }else{
      lrtest[is,1]    = 0
    }
    suplr   <- as.matrix(max(lrtest))
    news    <- which.max(lrtest)
    newd    <- ds[news,1]
  }
  colnames(suplr) <- paste0("supSeq(",m+1,",",n,"|",m,",",n,")")
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,q,m,trm)$cvSeq)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,q,m,trm)$cvSeq,getcv(0.05,q,m,trm)$cvSeq, 
                        getcv(0.025,q,m,trm)$cvSeq, getcv(0.01,q,m,trm)$cvSeq)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr9_out <- list(supSeq = suplr, newd = newd, cv = cv)
  class(pslr9_out) <- "test"
  return(pslr9_out)
}



#' @title Test stat for  n vs n+1 breaks in var given m breaks in mean
#' 
#' @description Computes the SeqLR_10T test statistic for n variance changes versus n + 1 variance changes given m coefficient changes
#' 
#' @details Note: This code is a translation of the one originally written by Pierre 
#' Perron and Yohei Yamamoto for MATLAB. Original code files can be found on 
#' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
#' 
#' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
#' 
#' @export
pslr10<- function(y, m, n, trm, z, x = matrix(0,0,0), control = list()){
  # ----- Set control values
  con <- list(vrobust = TRUE,
              prewhit = FALSE,
              typek = 2,
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
  h <- round(trm*bigt)
  # ----- 1)break date estimation under the null
  if (p==0){
    reg0 <- z
  }else{
    reg0 <-cbind(z, x)
  }
  res0  <- y - reg0%*%invpd(t(reg0)%*%reg0)$xinv%*%t(reg0)%*%y
  vvar0 <- t(res0)%*%res0/bigt
  tao0  <- (res0^2/c(vvar0))-1
  lr0   <- -(bigt/2)*(log(2*pi)+1+log(vvar0))
  
  brcase_tmp          <- numcase2(m,n)
  suplrx              <- matrix(0,brcase_tmp$num,1)
  brctemp             <- matrix(0,brcase_tmp$num,m)
  brvtemp             <- matrix(0,brcase_tmp$num,n)
  
  for (idx in 1:brcase_tmp$num){
    brcase_out    <- brcase_tmp$cvbrind[[idx]]
    K             <- ncol(brcase_out)
    cbrind        <- as.matrix(brcase_out[1,])
    vbrind        <- as.matrix(brcase_out[2,])
    estdate_out   <- estdate(y,z,q,x,p,K,bigt,h,m,n,cbrind,vbrind)
    brc           <- estdate_out$brc
    brv           <- estdate_out$brv
    res1          <- estdate_out$res
    brctemp[idx,] <- t(brc)
    brvtemp[idx,] <- t(brv)
    plog_out      <- ploglik(res1,n,brv)
    lr1           <- plog_out$loglik
    tao1          <- plog_out$tao
    lrv1          <- ploglik(res0,n,brv)$loglik
    if (con$vrobust==0){
      phi <- t(tao0)%*%tao0/(bigt-1)
    }else if (con$vrobust==1){
      phi <- correct1(tao0, tao1, con$prewhit, con$typek,con$kerntype)
    }
    suplrx[idx,1] <- 2*(lr1-lr0)-((phi-2)/phi)*2*(lrv1-lr0)
  }
  maxind          <- which.max(suplrx)
  brc0            <- as.matrix(brctemp[maxind,])
  brv0            <- as.matrix(brvtemp[maxind,])
  
  # ----- 2) construct segment-wise LR tests
  nseg          <- n+1
  lrtest        <- matrix(0,nseg,1)
  dv            <- matrix(0,nseg+1,1)
  if (length(brv0)!=0){
    dv[2:nseg,1]  <- brv0  
  }
  dv[nseg+1,1]  <- bigt
  ds            <- matrix(0,nseg,1)
  
  for (is in 1:nseg){
    lengthi     <- dv[(is+1),1] - dv[is,1]
    if (lengthi>=2*h){
      starti    <- dv[is,1]+1
      endi      <- dv[is+1,1]
      segy      <- as.matrix(y[starti:endi,1])
      segz      <- as.matrix(z[starti:endi,])
      if (p==0){
        segx      <- matrix(0,0,0)
      }else{
        segx      <- as.matrix(x[starti:endi,])
      }          
      indbrc      <- (brc0>starti)*(brc0<endi)
      mi          <- sum(indbrc)
      brci        <- matrix(0,0,0)
      if (m>=1){
        for (j in 1:m){
          if (indbrc[j,1]==1){
            brci    <- as.matrix(c(brci,brc0[j,1]-dv[is,1]))
          }
        } 
      }
      if (mi==0){
        segzbar   <- segz
      }else{
        segzbar   <- pzbar(segz, mi, brci)
      }
      if (p==0){
        segreg0   <- segzbar
      }else{
        segreg0   <- cbind(segzbar, segx)
      }
      segbeta0    <- invpd(t(segreg0)%*%segreg0)$xinv%*%t(segreg0)%*%segy
      segres0     <- segy - segreg0%*%segbeta0
      segvvar0    <- t(segres0)%*%segres0/lengthi
      seglr0      <- -(lengthi/2)*(log(2*pi)+1+log(segvvar0))
      segtao0     <- ((segres0^2)/c(segvvar0))-1
      
      bigveci     <- matrix(0,lengthi*2,1)
      for (i in 1:2){
        bigveci[((i-1)*lengthi+1):(i*lengthi),1] <- segres0^2
      }
      dateveci    <- dating_loglik(bigveci,h,1,lengthi)$datevec
      brvi        <- as.matrix(dateveci[1,1])
      segres1     <- estimbr(segy,segz,q,segx,p,lengthi,mi,brci,matrix(0,0,0),1,brvi,0)$res
      plog_out    <- ploglik(segres1,1,brvi)
      seglr1      <- plog_out$loglik
      segtao1     <- plog_out$tao
      
      # segment-wise LR test
      if (con$vrobust==0){
        segphi  <- t(segtao0)%*%segtao0/(lengthi-1)
      }else if (con$vrobust==1){
        segphi  <- correct1(segtao0,segtao1,con$prewhit,con$typek,con$kerntype)
      }
      
      lrtest[is,1]  <- (2/segphi)*2*(seglr1-seglr0)
      ds[is,1]      <- dv[is,1] + brvi
    }else{
      lrtest[is,1] <- 0
    }
    suplr   <- as.matrix(max(lrtest))
    news    <- which.max(lrtest)
    newd    <- ds[news,1]
  }
  
  colnames(suplr) <- paste0("supSeq(",m,",",n+1,"|",m,",",n,")")
  # get critical values 
  if (is.null(con$alpha)==FALSE){
    cv <- as.matrix(getcv(con$alpha,q,n,trm)$cvSeq)
    colnames(cv) <- paste0((1-con$alpha)*100,"%")
  }else{
    cv <- t(as.matrix(c(getcv(0.1,q,n,trm)$cvSeq,getcv(0.05,q,n,trm)$cvSeq, 
                        getcv(0.025,q,n,trm)$cvSeq, getcv(0.01,q,n,trm)$cvSeq)))
    colnames(cv) <- c("90%","95%","97.5%","99%")
  }
  # ----- Organize output
  pslr10_out <- list(supSeq = suplr, newd = newd, cv = cv)
  class(pslr10_out) <- "test"
  return(pslr10_out)
}


