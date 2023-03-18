#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title Compute SSR recursively
//' 
//' @description This procedure computes recursive residuals from a data set that starts
//' at date "start" and ends at date "last". It returns a vector of sum of
//' squared residuals (SSR) of length last-start+1 (stored for convenience in a vector of length T).
//' 
//' @param start starting entry of the sample used.
//' @param y dependent variable
//' @param z matrix of regressors of dimension q
//' @param h minimal length of a segment
//' @param last ending date of the last segment considered.
//' 
//' @details Note: This code is an adaptation of the one originally written by Yohei 
//' Yamamoto and Pierre Perron for MATLAB. Original codes can be found on 
//' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
//' 
//' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
//' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
//' 
//' @export
// [[Rcpp::export]]
arma::vec ssr_vec(int start, arma::vec y, arma::mat z, int h, int last){
  // initialize the vectors
  arma::vec vecssr(last, arma::fill::zeros);
  // initialize the recursion with the first h data points.
  arma::mat inv1    = inv(trans(z.rows(start-1,start+h-2))*z.rows(start-1,start+h-2));
  arma::mat delta1  = inv1*(trans(z.rows(start-1,start+h-2))*y.subvec(start-1,start+h-2));
  arma::mat res     = y.subvec(start-1,start+h-2) - z.rows(start-1,start+h-2)*delta1;
  vecssr(start+h-2) = as_scalar(trans(res)*res);
  // loop to construct the recursive residuals and update the SSR
  for (int r = start+h-1; r<last; r++){
    arma::mat v       = y(r) - z.row(r)*delta1;
    arma::mat invz    = inv1*trans(z.row(r));
    double f          = 1 + as_scalar(z.row(r)*invz);
    arma::mat delta2  = delta1 + (invz*v)/f;
    arma::mat inv2    = inv1 - (invz*trans(invz))/f;
    inv1              = inv2;
    delta1            = delta2;
    vecssr(r)         = as_scalar(vecssr(r-1) + v*v/f);
  }
  return(vecssr);
}


//' @export
// [[Rcpp::export]]
arma::vec ssrbigvec(arma::mat y, arma::mat z, int h){
  int bigt = y.n_rows;
  arma::vec bigvec(bigt*(bigt+1)/2, arma::fill::zeros);
  arma::vec ssr_vec_tmp;
  for (int i =1; i<=(bigt-h+1); i ++){
    ssr_vec_tmp = ssr_vec(i,y,z,h,bigt);
    bigvec.rows(((i-1)*bigt+i-(i-1)*i/2)-1,(i*bigt-(i-1)*i/2)-1) = ssr_vec_tmp.rows(i-1,bigt-1);
  }
  return(bigvec);
}


//' @title MLE function
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' @export
// [[Rcpp::export]]
arma::vec mlef(int start, arma::vec y, arma::mat z, int q, arma::mat x, int p, int h, int last){
  Rcpp::Environment mbreaks("package:mbreaks");
  Rcpp::Function invpd = mbreaks["invpd"];
  double pi = arma::datum::pi;
  arma::vec loglik(last, arma::fill::zeros);
  if (p>0){
    y = y - x*inv(trans(x)*x)*trans(x)*y;
    z = z - x*inv(trans(x)*x)*trans(x)*z;
  }
  int i       = start-1;
  int j       = start+h-2;
  arma::mat segz    = z.rows((i-1)+1,j);
  arma::vec segy    = y.rows((i-1)+1,j);
  List segz2invpd_ls = invpd(trans(segz)*segz);
  arma::mat segz2invpd = segz2invpd_ls["xinv"];
  arma::vec b       = segz2invpd*trans(segz)*segy;
  arma::vec  res    = segy - segz*b;
  arma::mat vvar    = trans(res)*res/h;
  arma::mat vstar   = vvar + 1;
  arma::mat diagmat(segy.n_rows,segy.n_rows,arma::fill::eye);
  int itr     = 1;
  List vstar_ls;
  arma::mat vstarinpd;
  arma::mat ibigv;
  while ((max(abs(vectorise(vvar-vstar)))>1e-6) & (itr<1000)){
    vstar = vvar;
    vstar_ls = invpd(vstar);
    vstarinpd = as<arma::mat>(vstar_ls["xinv"]);
    ibigv = kron(diagmat, vstarinpd);
    segz2invpd_ls = invpd(trans(segz)*ibigv*segz);
    segz2invpd = as<arma::mat>(segz2invpd_ls["xinv"]);
    b       = segz2invpd*trans(segz)*ibigv*segy;
    res    = segy - segz*b;
    vvar    = trans(res)*res/h;
    itr   = itr + 1;
  }
  List vvar_ls = invpd(vvar);
  arma::mat vvarinvpd = vvar_ls["xinv"];
  arma::mat hdiagmat(h,h,arma::fill::eye);
  ibigv   = kron(hdiagmat,vvarinvpd);
  List ihstar_ls = invpd(trans(segz)*ibigv*segz);
  arma::mat ihstar = ihstar_ls["xinv"];
  arma::vec bstar   = b;
  arma::vec ystar   = segy;
  arma::mat zstar   = segz;
  arma::mat gstar;
  arma::mat tempz;
  arma::vec tempy;
  arma::mat icstar;
  List icstar_ls;
  vstar   = vvar;
  double idub = i;
  while (j<last){ // should we start at j+1 (currently counting value at j twice...)
    gstar   = trans(ihstar*z.rows((j-1)+1,j));
    icstar_ls = invpd(vstar + z.rows((j-1)+1,j)*ihstar*trans(z.rows((j-1)+1,j)));
    icstar = as<arma::mat>(icstar_ls["xinv"]);
    b   = bstar + gstar*icstar*(y.rows((j-1)+1,j) - z.rows((j-1)+1,j)*bstar);
    tempz   = join_cols(zstar, z.rows((j-1)+1,j));
    tempy   = join_cols(ystar, y.rows((j-1)+1,j));
    res     = tempy - tempz*b;
    vvar    = trans(res)*res/(j-i+1);
    itr     = 1;
    vstar   = vvar + 10;
    while ((max(abs(vectorise(vvar-vstar)))>1e-6) & (itr<1000)){
      vstar = vvar;
      icstar_ls = invpd(vstar + z.rows((j-1)+1,j)*ihstar*trans(z.rows((j-1)+1,j)));
      icstar = as<arma::mat>(icstar_ls["xinv"]);
      b   = bstar + gstar*icstar*(y.rows((j-1)+1,j) - z.rows((j-1)+1,j)*bstar);
      res     = tempy - tempz*b;
      vvar    = trans(res)*res/(j-i+1);
      itr   = itr + 1;
    }
    double jdub = j; 
    double lentmp = (jdub-idub+1)/2;
    loglik.row(j) = (lentmp)*(log(2*pi)+1)+lentmp*log(vvar);
    bstar     = b;
    zstar     = tempz;
    ystar     = tempy;
    vstar     = vvar;
    ihstar    = ihstar - gstar*icstar*trans(gstar);
    j = j + 1;
  }
  return(loglik);
}




//' @title MLE function
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' @export
// [[Rcpp::export]]
arma::vec mlebigvec(arma::vec y, arma::mat z, int q, arma::mat x, int p, int h, int bigt){
  arma::vec bigvec(bigt*(bigt+1)/2, arma::fill::zeros);
  arma::vec loglik;
  for (int i =1; i<=(bigt-h+1); i ++){
    loglik = mlef(i,y,z,q,x,p,h,bigt);
    bigvec.rows(((i-1)*bigt+i-(i-1)*i/2)-1,(i*bigt-(i-1)*i/2)-1) = loglik.rows(i-1,bigt-1);
  }  
  return(bigvec);
}




//' @title optimal break partitions for a given segment
//' 
//' @description procedure to obtain an optimal one break partitions for a segment that 
//' starts at date start and ends at date last. It returns the optimal break and the 
//' associated SSR.
//' 
//' @param start: beginning of the segment considered
//' @param b1: first possible break date
//' @param b2: last possible break date
//' @param last: end of segment considered
//' 
//' @details Note: This code is translated from MATLABS code written by Yohei 
//' Yamamoto and Pierre Perron. Original codes can be found on 
//' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
//' 
//' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
//' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
//' 
//' @export
// [[Rcpp::export]]
List parti(int start, int b1, int b2, int last, arma::vec bigvec, int bigt){
  arma::vec dvec(bigt, arma::fill::zeros);
  int ini = (start-1)*bigt-(start-2)*(start-1)/2+1;
  int jj;
  int k;
  for (int j = b1; j<=b2; j++){
    jj = j-start;
    k = j*bigt-(j-1)*j/2+last-j;
    dvec(j-1) = bigvec((ini+jj)-1) + bigvec(k-1);
  }
  double ssrmin = min(dvec.rows(b1-1,b2-1));
  int minindcdvec = dvec.rows(b1-1,b2-1).index_min();
  int dx = (b1-1) + minindcdvec + 1;
  List output;
  output["ssrmin"] = ssrmin;
  output["dx"] = dx;
  return(output);
}

//' @title optimal break partitions for a given segment using log-likelihood
//' 
//' @description procedure to obtain an optimal one break partitions for a segment that 
//' starts at date start and ends at date last. It returns the optimal break and the 
//' associated log-likelihood (parti2() in original MATLAB code)
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' 
//' @export
// [[Rcpp::export]]
List parti_loglik(int start, int b1, int b2, int last, arma::vec bigvec, int bigt){
  double pi = arma::datum::pi;
  arma::vec dvec(bigt, arma::fill::zeros);
  double llr1;
  double llr2;
  for (int j = b1; j<=b2; j ++){
    llr1 = -0.5*(j-start+1)*((log(2*pi)+1)+log(sum(bigvec.rows(start-1,j-1))/(j-start+1)));
    llr2 = -0.5*(last-j)*((log(2*pi)+1)+log(sum(bigvec.rows((1*bigt+j+1)-1,(1*bigt+last)-1))/(last-j)));
    dvec(j-1) = llr1 + llr2;
  }
  double lrmax = max(dvec.rows(b1-1,b2-1));
  int maxindcdvec = dvec.rows(b1-1,b2-1).index_max();
  int dx = (b1-1) + maxindcdvec + 1; 
  List output;
  output["lrmax"] = lrmax;
  output["dx"] = dx;
  return(output);
}

//' @title Compute global break dates for pure structural change model
//' 
//' @description This is the main procedure which calculates the break points that globally
//'  minimizes the SSR. It returns optimal dates and associated SSR for all numbers of breaks less than or equal to m.
//' 
//' @details Note: This code is an adaptation of the one originally written by Yohei 
//' Yamamoto and Pierre Perron for MATLAB. Original code files can be found on 
//' Pierre Perron's website: https://blogs.bu.edu/perron/codes/
//' 
//' @param y A (\code{T x 1}) vector with endogenous variable.
//' @param z A (\code{T x q}) matrix with explanatory variables subject to change.
//' @param m An integer determining the number of breaks to find.
//' @param h An integer determining the minimum length of a regime. 
//' 
//' @references Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.
//' @references Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
//' 
//' @export
// [[Rcpp::export]]
List dating_purescSSR(arma::vec y, arma::mat z, int m, int h){
  // ----- Set some values
  int bigt  = y.n_rows;
  arma::mat datevec(m, m, arma::fill::zeros);
  arma::mat optdat(bigt, m, arma::fill::zeros);
  arma::mat optssr(bigt, m, arma::fill::zeros);
  arma::vec glb(m, arma::fill::zeros);
  arma::vec dvec(bigt, arma::fill::zeros);
  arma::vec bigvec = ssrbigvec(y, z, h);
  // Determine global optimum break dates
  if (m==1){
    List ssrmin_datx = parti(1, h, bigt-h, bigt, bigvec, bigt);
    datevec(0,0) = ssrmin_datx["dx"];
    glb(0) = ssrmin_datx["ssrmin"];
  }else{
    for (int j1 = (2*h); j1<=bigt; j1++){
      List ssrmin_datx = parti(1, h, j1-h, j1, bigvec, bigt);
      optdat(j1-1,0) = ssrmin_datx["dx"];
      optssr(j1-1,0) = ssrmin_datx["ssrmin"];
    }
    datevec(0,0) = optdat(bigt-1,0);
    glb(0) = optssr(bigt-1,0);
    for (int ib = 2; ib<=m; ib++){
      if (ib==m){
        int jlast = bigt;
        for (int jb = (ib*h); jb<=(jlast-h); jb++){
          dvec(jb-1) = optssr(jb-1,ib-2) + bigvec(((jb+1)*bigt-jb*(jb+1)/2)-1);
        }
        optssr(jlast-1,ib-1) = min(dvec.rows((ib*h)-1,(jlast-h)-1));
        int minindcdvec =  dvec.rows((ib*h)-1,(jlast-h)-1).index_min();
        optdat(jlast-1,ib-1) = (ib*h-1) + minindcdvec + 1;
      }else{
        for (int jlast = ((ib+1)*h); jlast<=bigt; jlast++){
          for (int jb  = (ib*h); jb<=(jlast-h); jb++){
            dvec(jb-1) = optssr(jb-1,ib-2) + bigvec((jb*bigt-jb*(jb-1)/2+jlast-jb)-1);
          }
          optssr(jlast-1,ib-1) = min(dvec.rows((ib*h)-1,(jlast-h)-1));
          int minindcdvec = dvec.rows((ib*h)-1,(jlast-h)-1).index_min();
          optdat(jlast-1,ib-1) = (ib*h-1) + minindcdvec + 1;
        }
      }
      datevec(ib-1,ib-1) = optdat(bigt-1,ib-1);
      for (int i = 1; i<=(ib-1); i++){
        int xx = ib-i;
        datevec(xx-1,ib-1) = optdat(datevec(xx,ib-1)-1,xx-1);
      }
      glb(ib-1) = optssr(bigt-1,ib-1);
    }
  }
  // organize output 
  List output; 
  output["glb"] = glb;
  output["datevec"] = datevec;
  output["bigvec"] = bigvec;
  return(output);
}



//' @title Compute global break dates using log-likelihood 
//' 
//' @description This is the main procedure which calculates the break points that globally maximize the loglikelihood function. It returns optimal dates and associated log likelihood for all numbers of breaks less than or equal to m. 
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' 
//' @export
// [[Rcpp::export]]
List dating_loglik(arma::vec bigvec, int h, int m, int bigt){
  double pi = arma::datum::pi;
  arma::mat datevec(m, m, arma::fill::zeros);
  arma::mat optdat(bigt, m, arma::fill::zeros);
  arma::mat optlr(bigt, m, arma::fill::zeros);                    
  arma::vec dvec(bigt, arma::fill::zeros);
  arma::vec glob(m, arma::fill::zeros);
  if (m==1){
    List dating_out = parti_loglik(1,h,bigt-h,bigt,bigvec,bigt);
    datevec(0,0)    = dating_out["dx"];
    glob(0)         = dating_out["lrmax"];
  }else{
    for (int j1 = (2*h); j1<=bigt; j1++){
      List dating_out      =  parti_loglik(1,h,j1-h,j1,bigvec.rows(0,(2*bigt)-1),bigt);
      optlr(j1-1,0)   = dating_out["lrmax"];
      optdat(j1-1,0)  = dating_out["dx"];
    }                                                                      
    datevec(0,0)  = optdat(bigt-1,0);
    glob(0)       = optlr(bigt-1,0);
    for (int ib = 2; ib<=m; ib++){
      if (ib==m){
        int jlast = bigt;
        for (int jb = (ib*h); jb<=(jlast-h); jb++){
          dvec(jb-1) = optlr(jb-1,ib-2) - 0.5*(bigt-jb+1)*((log(2*pi)+1)+log(sum(bigvec.rows((m*bigt+jb+1)-1,(bigt*(m+1))-1))/(bigt-jb)));
        }
        optlr(jlast-1,ib-1)  = max(dvec.rows((ib*h)-1,(jlast-h)-1));
        int maxindcdev = dvec.rows((ib*h)-1,(jlast-h)-1).index_max();
        optdat(jlast-1,ib-1) = (ib*h-1) + maxindcdev + 1;
      }else{
        for (int jlast = ((ib+1)*h); jlast<=bigt; jlast++){
          for (int jb = (ib*h); jb<=(jlast-h); jb++){
            dvec(jb-1) = optlr(jb-1,ib-2) - 0.5*(jlast-jb+1)*((log(2*pi)+1)+log(sum(bigvec.rows((ib*bigt+jb+1)-1,(ib*bigt+jlast)-1))/(jlast-jb)));
          }
          optlr(jlast-1,ib-1) = max(dvec.rows((ib*h)-1,(jlast-h)-1));
          int maxindcdev = dvec.rows((ib*h)-1,(jlast-h)-1).index_max();
          optdat(jlast-1,ib-1) = (ib*h-1) + maxindcdev + 1;
        }
      }
      datevec(ib-1,ib-1) = optdat(bigt-1,ib-1);
      for (int i = 1; i<=(ib-1); i++){
        int xx = ib-i;
        datevec(xx-1,ib-1) = optdat(datevec(xx,ib-1)-1, xx-1);
      }
      glob(ib-1) = optlr(bigt-1,ib-1);
    }
  }
  // organize output 
  List output; 
  output["glob"] = glob;
  output["datevec"] = datevec;
  return(output);
}



//' @title Compute global break dates using log-likelihood 
//' 
//' @description This is the main procedure which calculates the break points that globally maximize the loglikelihood function. It returns optimal dates and associated log likelihood for all numbers of breaks less than or equal to m. 
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' 
//' @export
// [[Rcpp::export]]
List dating_MLE(arma::vec y, arma::mat z, int q, arma::mat x, int p, int h, int m, int bigt){
  arma::mat datevec(m, m, arma::fill::zeros); 
  arma::mat optdat(bigt, m, arma::fill::zeros);
  arma::mat optmle(bigt, m, arma::fill::zeros);     
  arma::vec dvec(bigt, arma::fill::zeros);
  arma::vec glob(m, arma::fill::zeros);
  arma::vec bigvec = mlebigvec(y, z, q, x, p, h, bigt);
  if (m == 1){
    List parti_out = parti(1,h,bigt-h,bigt,bigvec,bigt);
    datevec(0,0)  = parti_out["dx"];
    glob(0)     = parti_out["ssrmin"];
  }else{
    for (int j1 = (2*h); j1<=bigt; j1++){
      List parti_out = parti(1,h,j1-h,j1,bigvec,bigt);
      optmle(j1-1,0)  = parti_out["ssrmin"];
      optdat(j1-1,0)  = parti_out["dx"];
    }
    datevec(0,0)   = optdat(bigt-1,0);
    glob(0)      = optmle(bigt-1,0);
    for (int ib = 2; ib <=m; ib++){
      if (ib == m){
        int jlast = bigt;
        for (int jb = (ib*h); jb<=(jlast-h); jb++){
          dvec(jb-1)  = optmle(jb-1,ib-2) + bigvec(((jb+1)*bigt-jb*(jb+1)/2)-1);
        }
        optmle(jlast-1,ib-1)    = min(dvec.rows((ib*h)-1,(jlast-h)-1));
        int  minindcdvec        = dvec.rows((ib*h)-1,(jlast-h)-1).index_min();
        optdat(jlast-1,ib-1)    = (ib*h-1) + minindcdvec + 1;
      }else{
        for (int jlast = ((ib+1)*h); jlast<=bigt; jlast++){
          for (int jb = (ib*h); jb<=(jlast-h); jb++){
            dvec(jb-1)   = optmle(jb-1,ib-2) + bigvec((jb*bigt-jb*(jb-1)/2+jlast-jb)-1);
          }
          optmle(jlast-1,ib-1)  = min(dvec.rows((ib*h)-1,(jlast-h)-1));
          int minindcdvec       = dvec.rows((ib*h)-1,(jlast-h)-1).index_min();
          optdat(jlast-1,ib-1)  = (ib*h-1) + minindcdvec + 1;
        }
      }
      datevec(ib-1,ib-1) = optdat(bigt-1,ib-1);
      for (int i = 1; i<=(ib-1); i++){
        int xx = ib-i;
        datevec(xx-1,ib-1) = optdat(datevec(xx,ib-1)-1, xx-1);
      }
      glob(ib-1) = optmle(bigt-1,ib-1);
    }
  }
  glob = -glob;
  // organize output 
  List output; 
  output["glob"] = glob;
  output["datevec"] = datevec;
  output["bigvec"] = bigvec;
  return(output);
}



