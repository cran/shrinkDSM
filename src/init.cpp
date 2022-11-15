// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <stdexcept>
#include <string>
using namespace Rcpp;

void init(arma::vec& it,
          arma::vec& tau_obs,
          arma::vec& int_tau,
          arma::vec& cen_tau,
          arma::mat& Z_tau,
          arma::vec& G_tau,
          arma::vec& lam_st,
          arma::vec& nev, // Up to here are "outputs"
          arma::vec& y,
          const arma::vec& delta,
          const arma::vec& G,
          const arma::mat& Z,
          const bool grouped,
          const bool tv_inputs){

  int n = y.n_elem; // number of individuals
  int n_int = it.n_elem; // number of intervals

  if(it(n_int-1)<arma::max(y)){

    int orig_size = it.size();
    it.resize(orig_size + 1);
    it(orig_size) = arma::max(y);
    n_int++;

  }


  int dim_Z = Z.n_cols; // number of covariates

  // vector of lengths of intervals
  arma::vec int_length(it.n_elem, arma::fill::zeros);
  int_length(0) = it(0);
  for (int j = 1; j < n_int; j++){
    int_length(j) = it(j) - it(j-1);
  }

  arma::vec nrisk(n_int, arma::fill::zeros), ncen(n_int, arma::fill::zeros);
  arma::umat ind_obs(n_int, n, arma::fill::zeros), ind_surv(n_int, n, arma::fill::zeros);
  nev = arma::vec(n_int, arma::fill::zeros);
  nrisk(0) = n;
  for(int j = 0; j < n_int; j++){


    /*
     indicator if individual exits in period j
     rows are time periods; columns are individuals
     */
    if(j == 0){
      ind_obs.row(j) = (y <= it(0)).t();
    }else if(j==(n_int-1)){
      ind_obs.row(j) = (y>it(j-1)).t();
    }else{
      ind_obs.row(j) = ((y <= it(j)) % (y > it(j-1))).t();
    }

    /*(ind_st, ind_end
     ind_surv... indicator if individual survives period j
     rows are time periods; columns are individuals
     nrisk... number of individuals at risk at the beginning of period j
     */
    if (j < (n_int - 1)){
      ind_surv.row(j) = (y > it(j)).t();
      nrisk(j+1) = nrisk(j) - arma::sum(ind_obs.row(j));
    }

    /*
     ncen...number of individuals who exit in period j and are censored
     nev ...number of individuals who exit in period j and are uncensored
     */

    ncen(j) = arma::as_scalar(ind_obs.row(j) * (delta == 0));
    nev(j) = arma::as_scalar(ind_obs.row(j) * (delta == 1));

  }


  int n_tau = arma::sum(nrisk);
  tau_obs = arma::vec(n_tau, arma::fill::zeros);
  int_tau = arma::vec(n_tau, arma::fill::zeros);
  cen_tau = arma::vec(n_tau, arma::fill::zeros);
  lam_st = arma::vec(n_tau, arma::fill::zeros);
  G_tau = arma::vec(n_tau, arma::fill::zeros);

  if (tv_inputs == false) Z_tau = arma::mat(n_tau, dim_Z, arma::fill::zeros);



  int ind_st = 0;
  arma::vec obs_time(n_int, arma::fill::zeros) ;

  for (int j = 0; j < n_int; j++){
    // beginning of period j
    int beg;
    if (j == 0){
      beg = 0;
    } else {
      beg = it(j-1);
    }

    // number of individuals who exit in period j
    int nj = nev(j) + ncen(j);

    int nriskj;
    if (j < (n_int - 1)){
      nriskj = nrisk(j+1);
    } else {
      nriskj = 0;
    }
    int ind_end;

    if (nj > 0){
      // indicator for individuals who exit in period j
      arma::uvec h = (arma::find(ind_obs.row(j) == 1));

      /*
       end indicator is starting indicator +
       number of individuals in period j -
       1 (because starting indicator is 1 individual)
       */
      ind_end = ind_st + nj - 1;

      // Residual observed survival time in period j for individuals who exit in period j
      tau_obs(arma::span(ind_st, ind_end)) = y.elem(h) - beg;

      // Indicator for period in which current row is in
      int_tau(arma::span(ind_st, ind_end)) = j * arma::ones(nj);
      // Indicator if row is censored in current period (survivors of period j )
      cen_tau(arma::span(ind_st, ind_end)) = 1 - delta(h);

      if (!tv_inputs) Z_tau.rows(ind_st, ind_end) = Z.rows(h);
      if (grouped && !tv_inputs) G_tau.rows(ind_st, ind_end) = G.rows(h);

      obs_time(j) = arma::sum(tau_obs(arma::span(ind_st, ind_end)));
      ind_st = ind_end + 1;
    }

    if (nriskj > 0){
      arma::uvec h = (arma::find(ind_surv.row(j) == 1));
      ind_end = ind_st + nriskj - 1;

      tau_obs(arma::span(ind_st, ind_end)) = int_length(j) * arma::ones(nriskj);

      int_tau(arma::span(ind_st, ind_end)) = j * arma::ones(nriskj);
      cen_tau(arma::span(ind_st, ind_end)) = arma::ones(nriskj);

      if (!tv_inputs) Z_tau.rows(ind_st, ind_end) = Z.rows(h);
      if(grouped && !tv_inputs) G_tau.rows(ind_st, ind_end) = G.rows(h);

      ind_st = ind_end + 1;
      obs_time(j) = obs_time(j) + int_length(j) * nriskj;


    }
  }

  int i1 = 0;
  arma::vec haz(n_int, arma::fill::zeros);

  for (int j = 0; j < n_int; j++){
    int i2 = i1 + nrisk(j) - 1;
    if (obs_time(j) == 0){
      haz(j) = arma::as_scalar(haz(haz > 0));
    } else {
      if (nev(j) > 0 ){
        haz(j) = nev(j)/obs_time(j);
      } else {
        haz(j) = 0.1/obs_time(j);
      }
      lam_st(arma::span(i1, i2)) = haz(j) * arma::ones(nrisk(j));
      i1 = i2 + 1;
    }
  }

  if(tv_inputs) {
    Z_tau = Z;
    G_tau = G;
    int rows_Z = Z_tau.n_rows;
    if ( rows_Z != n_tau ) {
      throw std::invalid_argument(
          "Time-varying input data has incorrect dimensions. Expected: " +
            std::to_string(n_tau) + " got: " + std::to_string(rows_Z) + ".");
    }
    if ( grouped ) {
      int rows_G = G_tau.n_rows;
      if ( rows_G != n_tau ) {
        throw std::invalid_argument("Grouping input has incorrect dimensions. Expected: " +
            std::to_string(n_tau) + " got: " + std::to_string(rows_G) + ".");
      }
    }
  }
}

// [[Rcpp::export]]
Rcpp::List init_export(arma::vec& it,
          arma::vec& y,
          const arma::vec& delta,
          const arma::vec& G,
          const arma::mat& Z,
          arma::mat& Z_tau,
          const bool grouped,
          bool tv_inputs) {

  arma::vec tau_obs;
  arma::vec int_tau;
  arma::vec cen_tau;
  arma::vec G_tau;
  // Log hazard
  arma::vec log_lam_samp;
  arma::vec nev;

  init(it, tau_obs, int_tau, cen_tau, Z_tau, G_tau, log_lam_samp, nev, y, delta, G, Z, grouped, tv_inputs);

  return Rcpp::List::create(_["tau_obs"] = tau_obs,
                            _["int_tau"] = int_tau,
                            _["cen_tau"] = cen_tau,
                            _["G_tau"] = G_tau,
                            _["Z_tau"]= Z_tau,
                            _["log_lam_samp"] = log_lam_samp,
                            _["nev"] = nev);
}
