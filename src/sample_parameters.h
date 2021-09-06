#ifndef SAMPLE_PARAMETERS_H
#define SAMPLE_PARAMETERS_H

#include <RcppArmadillo.h>


void sample_alpha(arma::vec& beta_mean_samp,
                  arma::vec& theta_sr_samp,
                  arma::vec& y,
                  arma::colvec& tau2,
                  arma::colvec& xi2,
                  arma::vec& sig2,
                  arma::mat& Z_tau,
                  arma::mat& beta_nc_samp,
                  arma::field<arma::uvec>& finder,
                  int J,
                  int K,
                  int n_tau);

void sample_phi(arma::vec& phi_samp,
                arma::vec& X_demean_samp,
                arma::vec& var_samp,
                arma::vec& f_samp,
                arma::vec& sigma2_phi,
                arma::vec& int_tau,
                arma::field<arma::uvec>& G_finder,
                int G);

void resample_phi_g(arma::vec& phi_samp,
                    arma::vec& f_samp,
                    arma::vec& sigma2_phi,
                    arma::vec& h_samp,
                    int T,
                    int G);

void sample_f(arma::vec& f_samp,
              arma::vec& X_demean_samp,
              arma::vec& phi_tau_samp,
              arma::vec& var_samp,
              arma::vec& h_samp,
              arma::field<arma::uvec>& finder,
              int J
              );
                    
void sample_comps(arma::vec& mu_samp, 
                  arma::vec& var_samp, 
                  const arma::vec& eps_samp, 
                  const arma::vec& weights, 
                  const arma::vec& sigma2, 
                  const arma::vec& means, 
                  bool initial);
#endif
