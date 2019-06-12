//#define ARMA_USE_SUPERLU
//#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <stdio.h>
#include <time.h>
#include "aux_functions.h"

#define INF arma::datum::inf
#define FMIN(a, b) ((a) < (b) ? (a) : (b))
#define FMAX(a, b) ((a) < (b) ? (a) : (b))

using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



// [[Rcpp::export]]
List roprobit_internal(arma::sp_mat X,
                       arma::sp_mat XXinv,
                       int niterR,
                       int thinR,
                       Rcpp::NumericVector ChoiceSetLengthR,
                       Rcpp::NumericVector ROLLengthR) {
  // initialize variables
  arma::uword k=0, i=0, nIDs=ChoiceSetLengthR.length(),
    Nchoices_i=1, Nranked_i=1, r=1, iter=0, niter=niterR, thin=thinR;
  double upper_bound=INF, lower_bound=-INF, 
    MaxUnranked_i=-INF, MinRanked_i=INF,
    u=0, Xb_k=0;
  
  // generate helper variables
  int Nsamples = floor(niter / thin);
  arma::colvec beta = arma::zeros(X.n_cols,1);
  arma::mat betavalues = arma::zeros(Nsamples, X.n_cols);
  arma::colvec Y = X*beta;
  arma::colvec Xb = Y;
  //arma::mat XXinv = arma::spsolve(trans(X)*X, arma::eye(X.n_cols, X.n_cols)); // requires superLU solver ...
  arma::sp_mat Proj = XXinv * trans(X);
  arma::vec MaxUnranked(nIDs, fill::zeros); // cannot fill directly with INF
  arma::vec MinRanked(nIDs, fill::zeros);
  for (i=0; i<nIDs; i++) {
    MaxUnranked[i] = -INF;
    MinRanked[i] = INF;
  }
  arma::uvec ChoiceSetLength = Rcpp::as<arma::uvec>( ChoiceSetLengthR );
  arma::uvec ROLLength = Rcpp::as<arma::uvec>( ROLLengthR );

  for (iter=0; iter<niter; iter++) {
    k = 0;
    for (i=0; i<nIDs; i++) {
      // store variables to reduce memory access
      Nchoices_i = ChoiceSetLength[i];
      Nranked_i = ROLLength[i];
      MaxUnranked_i = ( Nranked_i<Nchoices_i ? Y.subvec((k+Nranked_i-1), (k+Nchoices_i-1)).min() : -INF);
      MinRanked_i = MinRanked[i];
      upper_bound = INF;
      lower_bound = -INF;
      for (r=1; r<=Nchoices_i; r++) {
        Xb_k = Xb[k];
        // generate bounds
        if (r==1) {
          upper_bound = INF;
          lower_bound = Y[k+1] - Xb_k;
        } else if (r<Nranked_i) {
          upper_bound = Y[k-1] - Xb_k;
          lower_bound = Y[k+1] - Xb_k;
        } else if (r==Nranked_i) {
          upper_bound = Y[k-1] - Xb_k;
          lower_bound = MaxUnranked_i - Xb_k;
        } else {
          upper_bound = MinRanked_i - Xb_k;
          lower_bound = -INF;
        }
        // draw truncated error terms
        u = ( lower_bound<upper_bound ? truncn2(0, 1, lower_bound, upper_bound) : 0 );
        // Rcpp::Rcout << "i=" << i << ", r=" << r << ", upper_bound=" << upper_bound << ", lower_bound=" << lower_bound << ", u=" << u << "\n";
        Y[k] = Xb_k + u;
        // update maximum unranked and minimum ranked utility
        if (r==Nranked_i) {
          MinRanked_i = Y[k];
          MinRanked[i] = MinRanked_i;
        }
        // increment vector position
        k += 1;
      }
    }
    
    // estimate linear model
    beta = Proj*Y;
    // Rcpp::Rcout << beta;
    if ( (iter % thin == 0) ) betavalues.row(iter/thin) = trans(beta);
    Xb = X*beta;
  }
  
  Rcpp::Rcout << "Done.\n";
  
  return List::create(  
    // parameter draws
    Named("betadraws") = betavalues,
    Named("foo") = 1
  );
}

