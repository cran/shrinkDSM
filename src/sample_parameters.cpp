// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <shrinkTVP.h>
using namespace Rcpp;

/*
 Sample static parameters in non-centered parameterization
 Regression parameters without the time variation. See `draw_state_pe_noncen` for the FFBS part
 (\tilde \beta in FS)
 First part is beta second part is sqrt theta
 */
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
                  int n_tau) {

  // Transform data into classic linear regression conditionals
  // Create Z \tilde\beta
  arma::mat W2(n_tau, K);
  for (int j = 0; j < J; j++){
    arma::uvec h = finder(j);
    arma::mat Z_curr = Z_tau.rows(h);
    W2.rows(h) = Z_curr.each_row() % beta_nc_samp.col(j+1).t();
  }

  /* Create design matrix [Z Z\tilde\beta]
   for static components of NC paramet.
   */
  arma::mat W = arma::join_rows(Z_tau, W2);

  arma::vec param_vec(2*K, arma::fill::zeros);

  arma::vec prior_var = arma::join_cols(tau2, xi2);
  shrinkTVP::sample_lin_reg_stab(param_vec,
                                 y,
                                 W,
                                 sig2,
                                 prior_var);

  std::for_each(param_vec.begin(), param_vec.end(), shrinkTVP::res_protector);

  beta_mean_samp = param_vec.rows(0, K-1);
  theta_sr_samp = param_vec.rows(K, 2*K-1);
}


void sample_f(arma::vec& f_samp,
              arma::vec& X_demean_samp,
              arma::vec& phi_tau_samp,
              arma::vec& var_samp,
              arma::vec& h_samp,
              arma::field<arma::uvec>& finder,
              int J){

  // Create storage objects as to prevent garbage collection
  arma::vec X_demean_curr;
  arma::vec phi_tau_curr;
  arma::vec var_samp_curr;
  arma::vec prior_var_curr(1, arma::fill::zeros);

  // This hack keeps me up at night
  arma::vec res(1, arma::fill::zeros);

  for (int j = 0; j < J; j++){
    arma::uvec h = finder(j);
    X_demean_curr = X_demean_samp.rows(h);
    phi_tau_curr = phi_tau_samp.rows(h);
    var_samp_curr = var_samp.rows(h);
    prior_var_curr(0) = std::exp(h_samp(j));


    // TO-DO look closer into what stochvol returns

    shrinkTVP::sample_lin_reg_stab(res,
                                   X_demean_curr,
                                   phi_tau_curr,
                                   var_samp_curr,
                                   prior_var_curr);
    f_samp(j) = res(0);
  }
}


void sample_phi(arma::vec& phi_samp,
                arma::vec& X_demean_samp,
                arma::vec& var_samp,
                arma::vec& f_samp,
                arma::vec& sigma2_phi,
                arma::vec& int_tau,
                arma::field<arma::uvec>& G_finder,
                int G){

  arma::uvec h;
  arma::vec X_demean_curr;
  arma::vec var_samp_curr;
  arma::uvec t;
  arma::vec f_curr;

  for (int g = 0; g < G; g++){
    h = G_finder(g);
    X_demean_curr = X_demean_samp.rows(h);
    var_samp_curr = var_samp.rows(h);
    t =  arma::conv_to<arma::uvec>::from(int_tau.rows(h));
    f_curr = f_samp.rows(t);

    double sigma2_phi_g_post =
      1/arma::as_scalar(arma::sum(arma::pow(f_curr, 2)/var_samp_curr) + 1/sigma2_phi(g));
    double mu_phi_g_post = sigma2_phi_g_post *
      arma::sum(X_demean_curr % f_curr / var_samp_curr);

    phi_samp(g) = R::rnorm(mu_phi_g_post, std::sqrt(sigma2_phi_g_post));
  }
}


void resample_phi_g(arma::vec& phi_samp,
                    arma::vec& f_samp,
                    arma::vec& sigma2_phi,
                    arma::vec& h_samp,
                    int T,
                    int G){

  arma::uword phi_m_h = (arma::abs(phi_samp)).index_max();
  double phi_m_samp = phi_samp(phi_m_h);

  arma::uvec phi_other_h = arma::find(phi_samp != phi_m_samp);
  arma::vec phi_star_other_samp = phi_samp.rows(phi_other_h) * 1/phi_m_samp;
  arma::vec f_star_samp = f_samp * phi_m_samp;

  arma::vec sigma2_phi_other = sigma2_phi.rows(phi_other_h);

  double p = (G - T) * 0.5;
  double a = arma::sum(arma::pow(phi_star_other_samp, 2)/sigma2_phi_other) + 1/sigma2_phi(phi_m_h);
  double b = arma::sum(arma::exp(arma::log(arma::pow(f_star_samp, 2)) - h_samp));

  double phi_m_samp_new;
  phi_m_samp_new = std::sqrt(shrinkTVP::do_rgig1(p, b, a));
  phi_m_samp_new = std::copysign(phi_m_samp_new, phi_m_samp);
  phi_samp = phi_samp * phi_m_samp_new/phi_m_samp;
  f_samp = f_samp * phi_m_samp/phi_m_samp_new;

}

void sample_comps(arma::vec& mu_samp,
                  arma::vec& var_samp,
                  const arma::vec& eps_samp,
                  const arma::vec& weights,
                  const arma::vec& sigma2,
                  const arma::vec& means,
                  bool initial){
  int n = eps_samp.n_elem;
  int R = weights.n_elem;
  arma::mat sigh = arma::repmat(arma::sqrt(sigma2).t(), n, 1);
  arma::mat muh  = arma::repmat(means.t(), n, 1);
  arma::mat ph   = arma::repmat(weights.t(), n, 1);
  arma::mat xh   = arma::repmat(eps_samp, 1, R);

  arma::mat probtau(arma::size(ph));
  if (initial == true){
    probtau = ph;
  } else {
    arma::mat lpt = - arma::log(sigh) - 0.5 * arma::pow((xh-muh)/sigh, 2) + arma::log(ph) + 10;
    arma::vec maxes = arma::max(lpt, 1);
    arma::mat lh = lpt - arma::repmat(maxes, 1, R);
    arma::vec sums = arma::exp(lh) * arma::ones(R);
    probtau = arma::exp(lh)/arma::repmat(sums, 1, R);
  }
  arma::mat cumsums = arma::cumsum(probtau, 1);
  arma::mat u = arma::repmat(Rcpp::as<arma::vec>(Rcpp::runif(n, 0, 1)), 1, R);

  arma::uvec r = arma::sum(u > cumsums, 1);

  mu_samp = means(r);
  var_samp = sigma2(r);
}


