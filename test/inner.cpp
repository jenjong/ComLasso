
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export]]

arma::mat Res_l2_C (arma::mat y, arma::mat fx) 
{
  arma::mat z = y-fx;
  return(z);
}
// [[Rcpp::export]]
arma::mat d_getgcorr_l2  (arma::mat tX, arma::mat d_fx)
{
  arma::mat z = tX*d_fx;
  return(z);
}
