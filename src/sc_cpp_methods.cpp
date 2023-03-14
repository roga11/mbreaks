#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// // ==============================================================================
// //' @title Standard normal errors using box Muller 
// //' 
// //' @param \code{T} Integer determining the length of the process to be simulated
// //' @param \code{cov_mat}  q x q covariance matrix 
// //' 
// //' @keywords internal
// //' 
// //' @export
// // [[Rcpp::export]]
// arma::mat randSN(int T, arma::mat cov_mat){
//   int q = cov_mat.n_rows;
//   double pi = arma::datum::pi;
//   arma::mat U1(T, q, arma::fill::randu);
//   arma::mat U2(T, q, arma::fill::randu);
//   arma::mat eps = trans(trans(sqrt(-2*log(U1))%cos(2*pi*U2)));
//   arma::mat C = chol(cov_mat, "lower");
//   arma::mat eps_corr = trans(C*trans(eps));
//   return(eps_corr);
// }

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
arma::mat sc_datemat(arma::mat y, arma::mat z, int h){
  int Tsize = y.n_rows;
  arma::mat date_mat(Tsize, Tsize, arma::fill::zeros);
  for (int xs = 0; xs<(Tsize-h+1); xs++){
    arma::vec ssr_vec_tmp = ssr_vec(xs+1,y,z,h,Tsize);
    date_mat.submat(xs, xs, xs, Tsize-1) = trans(ssr_vec_tmp.subvec(xs,Tsize-1));
  }
  return(date_mat);
}


//' @title MLE function
//' 
//' @references Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020), "Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model" \emph{Quantitative Economics}, vol 11, 1019-1057.
//' @export
// [[Rcpp::export]]
arma::vec mlef(int start, arma::vec y, arma::mat z, int q, arma::mat x, int p, int h, int last){
  Rcpp::Environment mbreaks("package:mbreaks");
  Rcpp::Function invpd = mbreaks["invpd"];
  Rcpp::Function vec = mbreaks["vec"];
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
  while ((max(abs(vec(vvar-vstar)))>1e-6) & (itr<1000)){
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
    while ((max(abs(vec(vvar-vstar)))>1e-6) & (itr<1000)){
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

