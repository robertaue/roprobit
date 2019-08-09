//#define ARMA_USE_SUPERLU
#define ARMA_NO_DEBUG
//#define DEBUG
#include <RcppArmadillo.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "aux_functions.h"
#include "demean.h"
#include "truncnorm_approx.h"

#define INF arma::datum::inf
#define FMIN(a, b) ((a) < (b) ? (a) : (b))
#define FMAX(a, b) ((a) < (b) ? (a) : (b))
#define TRUNCNORM(mu, sigma, lower, upper) truncn2(mu, sigma, lower, upper)
//#define TRUNCNORM(mu, sigma, lower, upper) truncnorm_approx(runiform(gen), mu, sigma, lower, upper)

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
                       Rcpp::NumericMatrix initparm, // initial parameters
                       Rcpp::NumericVector ChoiceSetLengthR,
                       Rcpp::NumericVector ROLLengthR,
                       Rcpp::NumericVector GroupIDsR,
                       int nCores,
                       bool demeanY) {
  // initialize variables
  arma::uword nIDs=ChoiceSetLengthR.length(), niter=niterR, thin=thinR;
  arma::uword chunk5percent = floor( (double) niter / (double) 20);
  
  omp_set_num_threads(nCores);
#ifdef DEBUG
  std::cout.setf(std::ios::unitbuf); // for debugging
#endif
  
  // generate helper variables
  int Nsamples = floor(niter / thin);
  //arma::colvec beta = arma::zeros(X.n_cols,1);
  arma::colvec beta = Rcpp::as<arma::colvec>( initparm );
  arma::mat betavalues = arma::zeros(Nsamples, X.n_cols);
  arma::colvec Y = X*beta;
  arma::colvec Xb = Y;
  arma::Col<int> GroupIDs = Rcpp::as<arma::Col<int>>( GroupIDsR);
  //arma::mat XXinv = arma::spsolve(trans(X)*X, arma::eye(X.n_cols, X.n_cols)); // requires superLU solver ...
  arma::sp_mat Proj = XXinv * trans(X);
  arma::mat XXinv_dense(XXinv); // quick & dirty, need dense matrix for cholesky decomp.
  arma::vec MaxUnranked(nIDs, fill::zeros); // cannot fill directly with INF
  //arma::vec MinRanked(nIDs, fill::zeros);
  for (arma::uword i=0; i<nIDs; i++) {
    MaxUnranked[i] = -INF;
    //MinRanked[i] = INF;
  }
  arma::uvec ChoiceSetLength = Rcpp::as<arma::uvec>( ChoiceSetLengthR );
  arma::uvec ROLLength = Rcpp::as<arma::uvec>( ROLLengthR );
  arma::uvec StartPosition = cumsum(ChoiceSetLength) - ChoiceSetLength;
  /*Rcpp::Rcout << ChoiceSetLength.head(5) << std::endl;
  Rcpp::Rcout << ROLLength.head(5) << std::endl;
  Rcpp::Rcout << StartPosition.head(5) << std::endl;*/
  
  Rcpp::Rcout << "Drawing " << niter << " MC samples (each . = 5%%) ";
  
#pragma omp parallel
{
  // initialize thread-private variables
  arma::uword k=0, i=0, Nchoices_i=1, Nranked_i=1, r=1, iter=0;
  double upper_bound=INF, lower_bound=-INF, 
    MaxUnranked_i=-INF, MinRanked_i=INF,
    u=0, Xb_k=0;
  
  // initialize the random number generator in this thread
  std::random_device rd;  // a random number that be used to seed the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  // Uniform random number generator, see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
  std::uniform_real_distribution<double> runiform(0, 1); 
  
  // initialize latent errors U
  // initially, b=0 so that U can be drawn from N(0,1) distribution
  // I use an approximation to the quantiles of the normal distribution, 
  // see https://math.stackexchange.com/questions/28751/proof-of-upper-tail-inequality-for-standard-normal-distribution/28754#28754
#pragma omp single
  for (i=0; i<nIDs; i++) {
    Nchoices_i = ChoiceSetLength[i];
    for (r=1; r<=Nchoices_i; r++) {
      Y[k] = qnorm_approx(1 - (double) (2*r-1)/ (double) (2*Nchoices_i));
      //Y[k] = U[k];
      k++;
    }
  }
  
  for (iter=0; iter<niter; iter++) {
    
#ifdef DEBUG 
    printf("tid=%d: iter=%d\n", omp_get_thread_num(), iter);
#endif
    
#pragma omp for private(i)
    for (i=0; i<nIDs; i++) {
#ifdef DEBUG
      printf("tid=%d: i=%d\n", omp_get_thread_num(), i);
#endif
      k = StartPosition[i];
      // store variables to reduce memory access
      Nchoices_i = ChoiceSetLength[i];
      Nranked_i = ROLLength[i];
      MaxUnranked_i = ( Nranked_i<Nchoices_i ? Y.subvec((k+Nranked_i), (k+Nchoices_i-1)).max() : -INF);
      //MinRanked_i = MinRanked[i];
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
        u = ( lower_bound<upper_bound ? TRUNCNORM(0.0, 1.0, lower_bound, upper_bound) : upper_bound );
        //if (u<lower_bound || u>upper_bound) printf("iter=%d, i=%d, r=%d, lower_bound=%f, upper_bound=%f, u=%f\n", iter, i, r, lower_bound, upper_bound, u);
        // Rcpp::Rcout << "i=" << i << ", r=" << r << ", upper_bound=" << upper_bound << ", lower_bound=" << lower_bound << ", u=" << u << "\n";
        Y[k] = Xb_k + u;
        // update maximum unranked and minimum ranked utility
        if (r==Nranked_i) {
          MinRanked_i = Y[k];
          //MinRanked[i] = MinRanked_i;
        }
        // increment vector position
        k += 1;
      }
    }
    
#ifdef DEBUG
    printf("tid=%d: end iter=%d\n", omp_get_thread_num(), iter);
#endif
 
 // estimate linear model   
#pragma omp single
{
    if (demeanY) demean(Y.memptr(), GroupIDs.memptr(), Y.n_elem);
    beta = mvrnormArma(Proj*Y, XXinv_dense, 1); //beta = Proj*Y;
    if ( (iter % thin == 0) ) betavalues.row(iter/thin) = trans(beta);
    Xb = X*beta;
    
#ifndef DEBUG
    if ( !fmod(iter, chunk5percent) ) Rcpp::Rcout << '.';
#endif
} // end single
  } // end iter
} // end parallel
  
  
  Rcpp::Rcout << " Done.\n";
  
  return List::create(  
    // parameter draws
    Named("betadraws") = betavalues,
    Named("Y") = Y,
    Named("Xb") = Xb
  );
}

