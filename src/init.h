#ifndef INIT_H
#define INIT_H

#include <RcppArmadillo.h>

void init(arma::vec& it,
          arma::vec& tau_obs,
          arma::vec& int_tau,
          arma::vec& cen_tau,
          arma::mat& Z_tau,
          arma::vec& G_tau,
          arma::vec& lam_st,
          arma::vec& nev,
          arma::vec& y,
          const arma::vec& delta,
          const arma::vec& G,
          const arma::mat& Z,
          const bool grouped,
          const bool Z_given = false
          );

#endif
