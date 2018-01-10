
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export]]
double abs_j(double a)
{
  if (a>0) 
  {
    return(a); 
  } else {
    return(-a);  
  }
}
// [[Rcpp::export]]
List init_fun_C(IntegerVector idx_gs, IntegerVector idx_ge,
                       NumericVector grad_vec)
{
  IntegerVector result(2);
  int idx_gs_length = idx_gs.size();
  int istar_s, istar_e; 
  double v0;
  double v = R_NegInf;
  int j1, j2;
  for (int istar = 0 ; istar< idx_gs_length; istar ++)
  {
    istar_s = idx_gs[istar];
    istar_e = idx_ge[istar];
    for (int i=istar_s; i < istar_e ; i++)
    {
      for (int j=istar_s; j < istar_e ; j++)
      {
        if (i==j) continue;
        v0 = -grad_vec[i-1] + grad_vec[j-1];
        if (v0 > v) 
        {
          v = v0;
          j1 = i;
          j2 = j;
        }
      }
    }
  }
  result[0] = j1;
  result[1] = j2;
  return(List::create(Named("v")= v, Named("jvec") = result));
}

// [[Rcpp::export]]
List d2_fun_C( IntegerVector act_vec,
               IntegerVector idx_gs, IntegerVector idx_ge,
               NumericVector grad_vec,
               NumericVector corr_vec, NumericVector rderiv,
               double lambda, int a, double tol)
{
  IntegerVector result(2);
  int act_vec_length = act_vec.size();
  double v = R_PosInf;
  double v0;
  double cjj,djj;
  int istar;
  int istar_s, istar_e;
  int j1,j2;
  for (int i_a = 0 ; i_a < act_vec_length; i_a ++)
  {
    istar = act_vec[i_a];
    istar_s = idx_gs[istar-1];
    istar_e = idx_ge[istar-1];
    for (int i=istar_s; i < istar_e ; i++)
    {
      for (int j=istar_s; j < istar_e ; j++)
      {
        if (i == j ) continue;
        cjj = -corr_vec[i-1]+ corr_vec[j-1] - 2*rderiv[1+a+1-1];
        djj = 2*lambda + grad_vec[i-1] - grad_vec[j-1];
        v0 = djj/cjj;
        if ( (cjj>=0)&(djj<0)) break;
        if ((cjj<0)&(djj<0)) break;
        if ((cjj<0)&(djj>0)) v0 = R_PosInf;
        if (v0 < v) 
        {
          v = v0;
          j1 = i;
          j2 = j;
        }
      }
    }
  }
  
  double g1, g2;
  int  jstar1, jstar2;
  double j1_range_max, j1_range_min, j2_range_max, j2_range_min;
  double range_diff;
  int js1,js2;
  double mu_tmp = NA_REAL;
  double tmp;
  if (v<R_PosInf)
  {
    g1 = grad_vec[j1-1] + corr_vec[j1-1]*v;
    g2 = grad_vec[j2-1] + corr_vec[j2-1]*v;
// jstar1: positive activation ; jstar2: negative activation 
    if (g1 < g2)    
    {
      jstar1 = j1;
      jstar2 = j2;
    } else {
      jstar2 = j1;
      jstar1 = j2;
    }
    
    j1_range_max = lambda + rderiv[1+a+1-1]*v - g1;
    j1_range_min = -(lambda + rderiv[1+a+1-1]*v) - g1;
    j2_range_max = lambda + rderiv[1+a+1-1]*v - g2;
    j2_range_min = -(lambda + rderiv[1+a+1-1]*v) - g2;
    range_diff = abs(j1_range_min- j2_range_max);
    
    if ( range_diff  < tol )  
    {
        js1 = j1;
        js2 = j2;
    } else {
        js1 = j2;
        js2 = j1;
    }
    result[0] = jstar1;
    result[1] = jstar2;
// mu_tmp is computed
    mu_tmp = lambda + rderiv[1+a+1-1]*v -grad_vec[js2-1]-corr_vec[js2-1]*v;
// mu_tmp should be equal to the following value
// mu_tmp <- -lambda -rderiv[1+a+1]*v -grad_vec[js1] -corr_vec[js1]*v 
  }
  return(List::create(Named("v")= v, Named("jvec") = result,
                      Named("mu")= mu_tmp));
}

    
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
