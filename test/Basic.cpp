#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List lucefun1(IntegerVector n1, IntegerVector nvec,List a, List ilist, List jlist, NumericVector muvec,
                IntegerVector avec, NumericVector theta, NumericVector gamma, double nu)
                
{
  /*List isActiveList, NumericVector numActive)*/
  /*maxiteration*/
  int max_iter = 25;
  double eps = 1e-5;
  double diff_eps = 10000;
  int cal_iter = 0;
  /* initializing  */
  int p = muvec.size();
  NumericVector mu = rep(0.0, p);
  NumericVector old_mu = rep(0.0, p);
  
  for (int i=0; i<p; i++)
  {
  mu[i] = muvec[i];
  old_mu[i] = muvec[i];
  }
  
  List i_list = List(p);
  List j_list = List(p);
  
  for ( int i=0 ; i<p; i++)
  {
  i_list[i] = rep(0, avec[i]);
  j_list[i] = rep(0, avec[i]);
  SEXP tmp1 = ilist[i];
  i_list[i] = tmp1;  
  SEXP tmp2 = jlist[i];
  j_list[i] = tmp2;
  }
  
  int n = n1.size();
  List a_list = List(n);
  for ( int i=0 ; i<n; i++)
  {
  a_list[i] = rep(0, n1[i]);
  SEXP tmp3 = a[i];
  a_list[i] = tmp3;  
  }
  
  double cm = 0;
  double bm1 = 0;
  double bm10 = 0;
  double bm3 = 0;
  double sbm = 0;
  double inter_v1 = 0;
  int ip = 0;
  int u = 0;
  double sumvar = 0;
  
/* -------link information is required --------*/  
    /* ---- calculate am    ----*/  
    
  /*NumericVector amvec = rep(0.0,p);  
  for ( int i = 0 ; i<p ; i++)
  {
    amvec[i] = numActive[i]*nu/2;
  }
  */
  double am = (p-1)*nu/2;
  
  
  /* ---- calculate bm2    ----*/  
  /* note that bm2 does not depend on mu*/
  NumericVector bm2 = rep(0.0, p);
  int idx = 0;
  double tmp1 = 0;
  double tmp2 = 0;
  double tmp3 = 0;
  double tmp4 = 0;  
  
  
  for (int i = 0 ; i<(p-1) ; i++)
  {
    for (int j = (i+1) ; j<p ; j++)
    {
     /* add N set */  
     
     
     /* add N set */  
     tmp1 =  bm2[i];
     /*std::printf("bm2 i:  %f \n ", tmp);*/
     bm2[i] = tmp1 - nu*theta[idx] - gamma[idx];
     
     tmp2 =  bm2[j];
     /*std::printf("bm2 j:  %f \n ", tmp);*/
     bm2[j] = tmp2 + nu*theta[idx] + gamma[idx];
     
     idx = idx + 1;
    }
  }
  
  
 for (int iter = 0 ; iter < max_iter; iter++) 
 {
    for (int istar = 0 ; istar<(p) ; istar++)
    {
      /* ---- calculate bm1    ----*/  
      bm1 = 0;
      IntegerVector ivec = i_list[istar];
       /* -- ip declare --*/  
      ip = ivec.size();
    
      for (int i=0 ; i<ip ; i++)
      {
        /* -- u declare -- */  
        u = ivec[i]-1; /* ivec[i] is (location) index of R*/
        /*std::printf("u:  %d \n ", u);*/
        
        IntegerVector ir = a_list[u]; /* ir is (location) index of R*/
        for (int j=0; j<n1[u] ; j++ ) /* n1[i] is (length) index of R*/
        {
            bm10 = 0;
            for (int k=j; k<n1[u]; k++)
            {
            bm10 += mu[ (ir[k]-1) ];  
            }
          bm1 += 1/bm10; 
          if (istar == (ir[j]-1)) break;
          /*std::printf("ir_j:  %d \n ", ir[j]);
          std::printf("bm10:  %f \n ", bm10);*/
        }
      }
      
      
      /* -- cm declare -- */
      cm = -avec[istar];
      
     bm3 = 0;
    
      /* -------link information is required --------*/
      for (int i=0 ; i<p ; i++)
      {
        if (i!=istar) 
        {
          bm3 = bm3 -nu*mu[i];
        }
      }
    
      sbm = bm1 + bm2[istar] + bm3;
      inter_v1 = pow(sbm,2)-8*am*cm;
      mu[istar] =  (-sbm + sqrt(inter_v1))/(4*am);
      /*std::printf("istar:  %d \n ", istar);*/
    } 
    /* normalization */
    /*
    sumvar = mu[(p-1)];
    for (int i=0; i<p; i++)
    {
      mu[i]= mu[i]/sumvar;    
    }
    */
    diff_eps = max(abs(mu-old_mu));  
    
    /* convergence 
    for (int i=0; i<(p-1); i++)
    {
      if ( abs(mu[i]-old_mu[i]) > diff_eps ) diff_eps = abs(mu[i]-old_mu[i]);
    }
    */
    if (diff_eps<eps) break;

    for (int i=0; i<(p-1); i++)
    {
       old_mu[i] = mu[i];
    }
    
    cal_iter = iter;
    
 }
 
  /*return(bm1);*/
  return(List::create(Named("mu")= mu, Named("iter") = cal_iter));
  }



