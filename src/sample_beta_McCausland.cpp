// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <shrinkTVP.h>
using namespace Rcpp;

/*
McCausland Notation:
yₜ = Xₜβ + Zₜαₜ + Gₜuₜ, t = 1, ..., n
αₜ₊₁ = Wₜβ + Tₜαₜ + Hₜuₜ
α₁ ~ 𝒩(a₁, P₁)
*/

void sample_beta_McCausland(arma::mat& beta_nc,
                            const arma::vec& X,
                            const arma::mat& z,
                            const arma::vec& theta_sr,
                            const arma::vec& V,
                            const arma::vec& beta_mean,
                            const arma::field<arma::uvec>& finder,
                            int J,
                            int K) {
  // helpers
  arma::mat I_K = arma::eye(K, K);
  arma::mat theta_sr_diag = arma::diagmat(theta_sr);

  // Storage objects (for calculating alpha)
  arma::cube L_upper_store(K, K, J+1);
  arma::cube steptwo_store(K, K, J+1);
  arma::cube m_store(K, 1, J+1);

  // Augmented data
  arma::vec X_star = X - z * beta_mean;
  arma::mat z_star = z * theta_sr_diag;

  // Omega Mat
  arma::mat Omega = 2*I_K;
  arma::mat Omega_offdiag = -1 * I_K;

  // Vector c
  arma::vec c(K, arma::fill::zeros);

  // Sigma_inverse
  arma::mat Sigma_inv(K, K, arma::fill::zeros);

  // Calculate necessary objects for state 0
  arma::mat L_lower = shrinkTVP::robust_chol(Omega);
  arma::mat L_upper = L_lower.t();

  arma::mat steptwo = solve(arma::trimatl(L_lower), Omega_offdiag);
  arma::mat stepthree = steptwo.t() * steptwo;

  arma::mat a0 = solve(arma::trimatl(L_lower), c);
  arma::mat m = solve(arma::trimatu(L_upper), a0);

  // Store objects
  L_upper_store.slice(0) = L_upper;
  steptwo_store.slice(0) = steptwo;
  m_store.slice(0) = m;
  // Start forward loop
  for (int t = 1; t < J + 1; t++) {
    arma::uvec h = finder(t - 1);
    Omega = (z_star.rows(h)).t() * arma::diagmat(1/V(h)) * z_star.rows(h) + (1 + (t != J))*I_K;
    c = (z_star.rows(h)).t() *  arma::diagmat(1/V(h)) * X_star(h);
    Sigma_inv = Omega - stepthree;

    L_lower = shrinkTVP::robust_chol(Sigma_inv);
    L_upper = L_lower.t();

    steptwo = solve(arma::trimatl(L_lower), Omega_offdiag);

    stepthree = steptwo.t() * steptwo;

    m = arma::solve(arma::trimatu(L_upper), arma::solve(arma::trimatl(L_lower), (c - Omega_offdiag * m)));
    L_upper_store.slice(t) = L_upper;
    steptwo_store.slice(t) = steptwo;
    m_store.slice(t) = m;
  }

  arma::vec eps = Rcpp::rnorm(K, 0, 1);
  arma::mat l = arma::solve(arma::trimatu(L_upper_store.slice(J)), eps);
  beta_nc.col(J) = l + m_store.slice(J);

  for (int t = J-1; t >= 0; t--){
    eps = Rcpp::rnorm(K, 0, 1);
    arma::vec q = eps - steptwo_store.slice(t) * beta_nc.col(t+1);
    l = arma::solve(arma::trimatu(L_upper_store.slice(t)), q);
    beta_nc.col(t) = l + m_store.slice(t);
  }

  std::for_each(beta_nc.begin(), beta_nc.end(), shrinkTVP::res_protector);
}
