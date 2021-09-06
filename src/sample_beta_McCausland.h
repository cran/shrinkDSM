#ifndef SAMPLE_BETA_MCCAUSLAND_H
#define SAMPLE_BETA_MCCAUSLAND_H

#include <RcppArmadillo.h>

void sample_beta_McCausland(arma::mat& beta_nc, 
                            const arma::vec& X, 
                            const arma::mat& z, 
                            const arma::vec& theta_sr, 
                            const arma::vec& V, 
                            const arma::vec& beta_mean, 
                            const arma::field<arma::uvec>& finder, 
                            int J, 
                            int K) ;

#endif
