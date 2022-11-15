// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include "sample_beta_McCausland.h"
#include "init.h"
#include "sample_parameters.h"
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <math.h>
#include <stochvol.h>
#include <shrinkTVP.h>

using namespace Rcpp;

// [[Rcpp::export]]
List do_shrinkDSM(arma::vec y,
                  arma::mat z,
                  std::string mod_type,
                  arma::vec delta,
                  arma::vec it,
                  arma::vec group,
                  int niter,
                  int nburn,
                  int nthin,
                  double d1,
                  double d2,
                  double e1,
                  double e2,
                  arma::vec sigma2_phi,
                  bool learn_lambda2_B,
                  bool learn_kappa2_B,
                  double lambda2_B,
                  double kappa2_B,
                  bool learn_a_xi,
                  bool learn_a_tau,
                  double a_xi,
                  double a_tau,
                  bool learn_c_xi,
                  bool learn_c_tau,
                  double c_xi,
                  double c_tau,
                  bool a_eq_c_xi,
                  bool a_eq_c_tau,
                  double a_tuning_par_xi,
                  double a_tuning_par_tau,
                  double c_tuning_par_xi,
                  double c_tuning_par_tau,
                  double beta_a_xi,
                  double beta_a_tau,
                  double alpha_a_xi,
                  double alpha_a_tau,
                  double beta_c_xi,
                  double beta_c_tau,
                  double alpha_c_xi,
                  double alpha_c_tau,
                  double Bsigma_sv,
                  double a0_sv,
                  double b0_sv,
                  bool display_progress,
                  arma::vec adaptive,
                  arma::vec target_rates,
                  arma::vec max_adapts,
                  arma::ivec batch_sizes,
                  bool tv_inputs,
                  List phi_prior) {

  // progress bar setup
  arma::vec prog_rep_points = arma::round(arma::linspace(0, niter, 50));
  Progress p(50, display_progress);

  // K = number of covariates
  int K = z.n_cols;

  // G = number of groups
  // grouped = true if G > 1
  arma::vec groups = arma::unique(group);
  int G = groups.n_elem;
  bool grouped;
  if(G>1){
    grouped = true;
  }else{
    grouped = false;
  }

  int nsave = std::floor((niter - nburn)/nthin);

  // Parameters of Mixture Components from FS&F 2006
  arma::vec weights = {1.6776755e-001, 1.4703656e-001, 1.2530633e-001, 1.1597262e-001, 1.0706560e-001, 1.0375856e-001,  1.0148516e-001,  8.8013418e-002,  3.9624418e-002,  3.9697961e-003 };
  arma::vec means = {1.8196721e+000, 1.2390286e+000, 7.6400964e-001, -3.0628675e-001, -6.7317508e-001, 4.3109958e-002, 3.9101728e-001, -1.0557713e+000, 3.2861388e+000, 5.0884679e+000};
  arma::vec sigma2 = {1.0989187e+000, 4.2227759e-001, 1.9753006e-001, 7.6600902e-002, 9.4682094e-002, 7.7834110e-002, 1.0689729e-001, 1.4602798e-001, 2.0192335e+000, 4.4966689e+000};

  // Create stacked observations and covariates and calculate necessary dimensions
  arma::vec tau_obs;
  arma::vec int_tau;
  arma::vec cen_tau;
  arma::vec G_tau;
  arma::mat Z_tau;
  // Log hazard
  arma::vec log_lam_samp;
  arma::vec nev;
  init(it, tau_obs, int_tau, cen_tau, Z_tau, G_tau, log_lam_samp, nev, y, delta, group, z, grouped, tv_inputs);
  const int J = it.n_elem;
  int n_tau = tau_obs.n_elem;

  arma::vec phi_tau_samp;
  arma::vec f_c_tau_samp;
  if (grouped) {
    phi_tau_samp = arma::vec(n_tau, arma::fill::zeros);
    f_c_tau_samp = arma::vec(n_tau, arma::fill::zeros);
  }

  // Create finder field as to not run arma::find in every iteration
  arma::field<arma::uvec> finder(J);
  for (int j = 0; j < J; j++){
    finder(j) = arma::find(int_tau == j);
  }

  // Create second finder for G_tau
  arma::field<arma::uvec> G_finder(G);
  if(grouped){
    for (int g = 0; g < G; g++){
      G_finder(g) = arma::find(G_tau == g);
    }
  }

  // Initial values for MCMC

  // Static beta in NC paramet.
  arma::vec beta_mean_samp(K);
  beta_mean_samp.fill(0.1);

  // Sqrt of variances in NC paramet.
  arma::vec theta_sr_samp(K);
  theta_sr_samp.fill(0.2);

  // Shrinkage adapatation parameter for xi2
  double a_xi_samp;
  if (learn_a_xi && (mod_type != "ridge")) {
    a_xi_samp = 0.1;
  } else {
    a_xi_samp = a_xi;
  }


  // Shrinkage adaptation parameter for tau2
  double a_tau_samp;
  if (learn_a_tau && (mod_type != "ridge")) {
    a_tau_samp = 0.1;
  } else {
    a_tau_samp = a_tau;
  }

  // Tail parameter for xi2
  double c_xi_samp;
  if (learn_c_xi && (mod_type == "triple")) {
    c_xi_samp = 0.1;
  } else {
    c_xi_samp = c_xi;
  }

  // Tail parameter for tau2
  double c_tau_samp;
  if (learn_c_tau && (mod_type == "triple")) {
    c_tau_samp = 0.1;
  } else {
    c_tau_samp = c_tau;
  }

  // Needed for data augmentation of TG sampler
  double d2_samp;
  double e2_samp;

  // lambda2/kappa2 data augmented for TG sampler
  arma::vec lambda2_til_samp(K);
  arma::vec kappa2_til_samp(K);
  if (mod_type == "triple") {
    lambda2_til_samp = arma::vec(K);
    lambda2_til_samp.fill(20);
    kappa2_til_samp = arma::vec(K);
    kappa2_til_samp.fill(20);
  }


  // Global shrinkage parameter of xi2
  double kappa2_B_samp;
  if (learn_kappa2_B && (mod_type != "ridge")) {
    kappa2_B_samp = 20;
  } else {
    kappa2_B_samp = kappa2_B;
  }

  // Global shrinkage parameter of tau2
  double lambda2_B_samp;
  if (learn_lambda2_B && (mod_type != "ridge")) {
    lambda2_B_samp = 20;
  } else {
    lambda2_B_samp = lambda2_B;
  }

  // Prior variances and alternate parameterization counterparts
  arma::vec tau2_samp(K);
  arma::vec xi2_samp(K);
  arma::vec tau2_til_samp(K);
  arma::vec xi2_til_samp(K);
  if (mod_type == "double") {
    tau2_samp = arma::vec(K, arma::fill::ones);
    xi2_samp = arma::vec(K, arma::fill::ones);
  } else if (mod_type == "triple") {
    tau2_til_samp = arma::vec(K, arma::fill::ones);
    xi2_til_samp = arma::vec(K, arma::fill::ones);

    shrinkTVP::calc_xi2_tau2(xi2_samp,
                             xi2_til_samp,
                             kappa2_til_samp,
                             kappa2_B_samp,
                             c_xi_samp,
                             a_xi_samp);

    shrinkTVP::calc_xi2_tau2(tau2_samp,
                             tau2_til_samp,
                             lambda2_til_samp,
                             lambda2_B_samp,
                             c_tau_samp,
                             a_tau_samp);
  } else {
    tau2_samp.fill(2.0/lambda2_B_samp);
    xi2_samp.fill(2.0/kappa2_B_samp);
  }

  // beta centered paramet.
  arma::mat beta_c_samp(K, J + 1, arma::fill::zeros);
  // time varying beta in nc paramet.
  arma::mat beta_nc_samp(K, J+1, arma::fill::zeros);

  // Component parameters, including eps_samp, which are inital draws from Gumbel
  arma::vec mu_samp(n_tau, arma::fill::zeros);
  arma::vec var_samp(n_tau, arma::fill::ones);
  arma::vec eps_samp =  arma::zeros(n_tau);
  sample_comps(mu_samp, var_samp, eps_samp, weights, sigma2, means, true);

  // Residual survival times
  arma::vec tau_res_samp(n_tau);
  for (int i = 0; i < n_tau; i++){
    tau_res_samp(i) = R::rexp(1/log_lam_samp(i));
  }
  tau_res_samp %= cen_tau;
  // Total survial times
  arma::vec tau_samp = tau_res_samp + tau_obs;

  // All necessary preliminaries for stochvol
  double h0_samp = 1;
  arma::vec h_samp(J, arma::fill::ones);
  arma::vec sv_para = {0,
                       .7,
  1};
  arma::uvec r(J); r.fill(5);
  using stochvol::PriorSpec;
  const PriorSpec prior_spec = {  // prior specification object for the update_*_sv functions
    PriorSpec::Latent0(),  // stationary prior distribution on priorlatent0
    PriorSpec::Mu(PriorSpec::Constant(0)),  // mu is constant
    PriorSpec::Phi(PriorSpec::Beta(a0_sv, b0_sv)),  // stretched beta prior on phi
    PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma_sv))  // normal(0, Bsigma) prior on sigma
  };  // heavy-tailed, leverage, regression turned off
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert {  // very expert settings for the Kastner, Fruehwirth-Schnatter (2014) sampler
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    1e-8,  // B011inv,
    1e-12,  //B022inv,
    3,  // MHsteps,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // unused for independence prior for sigma
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };

  // Values used for adaptive MH
  arma::mat batches;
  arma::vec curr_sds;
  arma::ivec batch_nrs;
  arma::ivec batch_pos;

  if (arma::sum(adaptive) > 0) {
    batch_pos = arma::ivec(4, arma::fill::zeros);
    batches = arma::mat(arma::max(batch_sizes), 4, arma::fill::zeros);
    curr_sds = {a_tuning_par_xi,
                a_tuning_par_tau,
                c_tuning_par_xi,
                c_tuning_par_tau};
    batch_nrs = arma::ivec(4, arma::fill::ones);
  }

  // Storage objects
  arma::cube beta_c_save(J+1, K, nsave);
  arma::mat theta_sr_save(K, nsave, arma::fill::none);
  arma::mat beta_mean_save(K, nsave, arma::fill::none);


  // Conditional storage objects
  arma::mat xi2_save;
  arma::mat tau2_save;
  if (mod_type != "ridge") {
    xi2_save = arma::mat(K, nsave, arma::fill::zeros);
    tau2_save = arma::mat(K, nsave, arma::fill::zeros);
  }

  arma::mat lambda2_save;
  arma::mat kappa2_save;
  if (mod_type == "triple") {
    lambda2_save = arma::mat(K, nsave, arma::fill::zeros);
    kappa2_save = arma::mat(K, nsave, arma::fill::zeros);
  }

  arma::vec kappa2_B_save;
  arma::vec lambda2_B_save;
  if (learn_kappa2_B && (mod_type != "ridge")) {
    kappa2_B_save = arma::vec(nsave, arma::fill::zeros);
  }
  if (learn_lambda2_B && (mod_type != "ridge")) {
    lambda2_B_save = arma::vec(nsave, arma::fill::zeros);
  }

  arma::vec a_xi_save;
  arma::vec a_tau_save;
  arma::vec c_xi_save;
  arma::vec c_tau_save;

  if (learn_a_xi && (mod_type != "ridge")) {
    a_xi_save = arma::vec(nsave, arma::fill::zeros);
  }
  if (learn_a_tau && (mod_type != "ridge")) {
    a_tau_save = arma::vec(nsave, arma::fill::zeros);
  }

  if (learn_c_xi && (mod_type == "triple")) {
    c_xi_save = arma::vec(nsave, arma::fill::zeros);
  }
  if (learn_c_tau && (mod_type == "triple")) {
    c_tau_save = arma::vec(nsave, arma::fill::zeros);
  }

  arma::vec a_xi_sd_save;
  arma::vec a_tau_sd_save;
  arma::vec c_xi_sd_save;
  arma::vec c_tau_sd_save;
  arma::vec a_xi_acc_rate_save;
  arma::vec a_tau_acc_rate_save;
  arma::vec c_xi_acc_rate_save;
  arma::vec c_tau_acc_rate_save;

  if (bool(adaptive(0)) && learn_a_xi && (mod_type != "ridge")) {
    a_xi_sd_save = arma::vec(std::floor(niter/batch_sizes(0)), arma::fill::zeros);
    a_xi_acc_rate_save = arma::vec(std::floor(niter/batch_sizes(0)), arma::fill::zeros);
  }

  if (bool(adaptive(1)) && learn_a_tau && (mod_type != "ridge")) {
    a_tau_sd_save = arma::vec(std::floor(niter/batch_sizes(1)), arma::fill::zeros);
    a_tau_acc_rate_save = arma::vec(std::floor(niter/batch_sizes(1)), arma::fill::zeros);
  }

  if (bool(adaptive(2)) && learn_c_xi && (mod_type == "triple")) {
    c_xi_sd_save = arma::vec(std::floor(niter/batch_sizes(2)), arma::fill::zeros);
    c_xi_acc_rate_save = arma::vec(std::floor(niter/batch_sizes(2)), arma::fill::zeros);
  }

  if (bool(adaptive(3)) && learn_c_tau && (mod_type == "triple")) {
    c_tau_sd_save = arma::vec(std::floor(niter/batch_sizes(3)), arma::fill::zeros);
    c_tau_acc_rate_save = arma::vec(std::floor(niter/batch_sizes(3)), arma::fill::zeros);
  }


  // Get all objects for factors/loadings and the associated shrinkage
  std::string mod_type_phi = as<std::string>(phi_prior["mod_type_phi"]);
  bool learn_a_phi = as<bool>(phi_prior["learn_a_phi"]);
  double a_phi = as<double>(phi_prior["a_phi"]);
  bool learn_c_phi = as<bool>(phi_prior["learn_c_phi"]);
  double c_phi = as<double>(phi_prior["c_phi"]);
  bool a_phi_eq_c_phi = as<bool>(phi_prior["a_phi_eq_c_phi"]);
  bool learn_lambda2_B_phi = as<bool>(phi_prior["learn_lambda2_B_phi"]);
  double lambda2_B_phi = as<double>(phi_prior["lambda2_B_phi"]);
  double e1_phi = as<double>(phi_prior["e1_phi"]);
  double e2_phi = as<double>(phi_prior["e2_phi"]);
  double beta_a_phi = as<double>(phi_prior["beta_a_phi"]);
  double alpha_a_phi = as<double>(phi_prior["alpha_a_phi"]);
  double beta_c_phi = as<double>(phi_prior["beta_c_phi"]);
  double alpha_c_phi = as<double>(phi_prior["alpha_c_phi"]);
  bool a_phi_adaptive = as<bool>(phi_prior["a_phi_adaptive"]);
  bool c_phi_adaptive = as<bool>(phi_prior["c_phi_adaptive"]);
  double a_phi_tuning_par = as<double>(phi_prior["a_phi_tuning_par"]);
  double c_phi_tuning_par = as<double>(phi_prior["c_phi_tuning_par"]);
  arma::vec target_rates_phi = as<arma::vec>(phi_prior["target_rates_phi"]);
  arma::vec max_adapts_phi = as<arma::vec>(phi_prior["max_adapts_phi"]);
  arma::ivec batch_sizes_phi = as<arma::ivec>(phi_prior["batch_sizes_phi"]);
  arma::vec adaptive_phi = as<arma::vec>(phi_prior["adaptive_phi"]);

  // Initial values for sampling
  // Factor loadings
  arma::vec phi_samp(G, arma::fill::ones);

  // Factor
  arma::vec f_samp(J, arma::fill::zeros);


  // Shrinkage adaptation parameter for phi
  double a_phi_samp;
  if (learn_a_phi && (mod_type_phi != "ridge")) {
    a_phi_samp = 0.1;
  } else {
    a_phi_samp = a_phi;
  }

  // Tail parameter for phi
  double c_phi_samp;
  if (learn_c_phi && (mod_type_phi == "triple")) {
    c_phi_samp = 0.1;
  } else {
    c_phi_samp = c_phi;
  }

  // Needed for data augmentation of TG sampler
  double e2_phi_samp;

  // lambda2 data augmented for TG sampler for phi
  arma::vec lambda2_til_phi_samp(G);
  if (mod_type_phi == "triple") {
    lambda2_til_phi_samp = arma::vec(G);
    lambda2_til_phi_samp.fill(20);
  }

  // Global shrinkage parameter of phi
  double lambda2_B_phi_samp;
  if (learn_lambda2_B_phi && (mod_type_phi != "ridge")) {
    lambda2_B_phi_samp = 20;
  } else {
    lambda2_B_phi_samp = lambda2_B_phi;
  }

  // Prior variances and alternate parameterization counterparts
  arma::vec tau2_phi_samp(G);
  arma::vec tau2_til_phi_samp(G);
  if (mod_type_phi == "double") {
    tau2_phi_samp = arma::vec(G, arma::fill::ones);
  } else if (mod_type_phi == "triple") {
    tau2_til_phi_samp = arma::vec(G, arma::fill::ones);

    shrinkTVP::calc_xi2_tau2(tau2_phi_samp,
                             tau2_til_phi_samp,
                             lambda2_til_phi_samp,
                             lambda2_B_phi_samp,
                             c_phi_samp,
                             a_phi_samp);
  } else {
    tau2_phi_samp.fill(2.0/lambda2_B_phi_samp);
  }

  // Values used for adaptive MH
  bool is_adaptive_phi = a_phi_adaptive || c_phi_adaptive;
  arma::mat batches_phi;
  arma::vec curr_sds_phi;
  arma::ivec batch_nrs_phi;
  arma::ivec batch_pos_phi;

  if (is_adaptive_phi) {
    batch_pos_phi = arma::ivec(2, arma::fill::zeros);
    batches_phi = arma::mat(arma::max(batch_sizes_phi), 2, arma::fill::zeros);
    curr_sds_phi = {a_phi_tuning_par,
                    c_phi_tuning_par};
    batch_nrs_phi = arma::ivec(2, arma::fill::ones);
  }




  arma::mat phi_save;
  arma::cube f_save;
  arma::cube h_save;


  if (grouped){
    f_save = arma::cube(J, 1, nsave, arma::fill::none);
    phi_save = arma::mat(G, nsave, arma::fill::none);
    h_save = arma::cube(J, 1, nsave, arma::fill::none);
  }
  // Begin Gibbs loop
  for (int iter = 0; iter < niter; iter++){

    // Augmented survival times
    arma::vec X_samp = -arma::log(tau_samp) - mu_samp;

    if (grouped) {
      X_samp -= f_c_tau_samp;
    }


    // Draw non-centered states
    sample_beta_McCausland(beta_nc_samp,
                           X_samp,
                           Z_tau,
                           theta_sr_samp,
                           var_samp,
                           beta_mean_samp,
                           finder,
                           J,
                           K);

    /*
     Sample static components of NC paramet.
     */

    sample_alpha(beta_mean_samp,
                 theta_sr_samp,
                 X_samp,
                 tau2_samp,
                 xi2_samp,
                 var_samp,
                 Z_tau,
                 beta_nc_samp,
                 finder,
                 J,
                 K,
                 n_tau);

    // Weave into centered parameterization
    shrinkTVP::to_CP(beta_c_samp,
                     beta_nc_samp,
                     theta_sr_samp,
                     beta_mean_samp);

    // Resample static components of NC paramet
    shrinkTVP::resample_alpha(beta_mean_samp,
                              theta_sr_samp,
                              beta_c_samp,
                              beta_nc_samp,
                              xi2_samp,
                              tau2_samp);

    shrinkTVP::to_NCP(beta_nc_samp,
                      beta_c_samp,
                      theta_sr_samp,
                      beta_mean_samp);

    bool True = true;
    std::string Fail = "fail";
    int one = 1;

    if (mod_type == "double") {
      shrinkTVP::sample_DG_TVP(beta_mean_samp,
                               theta_sr_samp,
                               tau2_samp,
                               xi2_samp,
                               lambda2_B_samp,
                               kappa2_B_samp,
                               a_xi_samp,
                               beta_a_xi,
                               alpha_a_xi,
                               a_tau_samp,
                               beta_a_tau,
                               alpha_a_tau,
                               d1,
                               d2,
                               e1,
                               e2,
                               learn_kappa2_B,
                               learn_lambda2_B,
                               learn_a_xi,
                               learn_a_tau,
                               a_tuning_par_xi,
                               a_tuning_par_tau,
                               adaptive,
                               batches,
                               curr_sds,
                               target_rates,
                               max_adapts,
                               batch_nrs,
                               batch_sizes,
                               batch_pos,
                               iter,
                               True, // These are values for shrinkTVPs error handling
                               Fail, // Could include later, if wanted
                               one);
    } else if (mod_type == "triple") {
      shrinkTVP::sample_TG_TVP(beta_mean_samp,
                               theta_sr_samp,
                               tau2_samp,
                               xi2_samp,
                               tau2_til_samp,
                               xi2_til_samp,
                               lambda2_til_samp,
                               kappa2_til_samp,
                               lambda2_B_samp,
                               kappa2_B_samp,
                               a_xi_samp,
                               beta_a_xi,
                               alpha_a_xi,
                               a_tau_samp,
                               beta_a_tau,
                               alpha_a_tau,
                               d2_samp,
                               e2_samp,
                               c_xi_samp,
                               c_tau_samp,
                               beta_c_xi,
                               alpha_c_xi,
                               beta_c_tau,
                               alpha_c_tau,
                               learn_kappa2_B,
                               learn_lambda2_B,
                               learn_a_xi,
                               learn_a_tau,
                               learn_c_xi,
                               learn_c_tau,
                               a_tuning_par_xi,
                               a_tuning_par_tau,
                               c_tuning_par_xi,
                               c_tuning_par_tau,
                               a_eq_c_xi,
                               a_eq_c_tau,
                               adaptive,
                               batches,
                               curr_sds,
                               target_rates,
                               max_adapts,
                               batch_nrs,
                               batch_sizes,
                               batch_pos,
                               iter,
                               True,
                               Fail,
                               one);
    }
    // Sample survival times
    /*
     To this end we first weave back into the centered parameterization and
     use this to create the log hazard rates (log_lam_samp)
     */
    shrinkTVP::to_CP(beta_c_samp,
                     beta_nc_samp,
                     theta_sr_samp,
                     beta_mean_samp);

    arma::uvec h_tau = arma::conv_to<arma::uvec>::from(int_tau + 1);
    arma::mat beta_c_tau = (beta_c_samp.cols(h_tau)).t();
    log_lam_samp = arma::sum(Z_tau % beta_c_tau, 1);

    /*
     In the case that groups are supplied, the latent factor and the group
     specific factor loadings have to be sampled and added to log_lam_samp.
     */
    if (grouped){
      // Create demeaned X to isolate factors and loadings on RHS
      arma::vec X_demean_samp = X_samp - log_lam_samp + f_c_tau_samp;

      // Sample the factor loadings
      sample_phi(phi_samp,
                 X_demean_samp,
                 var_samp,
                 f_samp,
                 tau2_phi_samp,
                 int_tau,
                 G_finder,
                 G);

      // fill up phi_tau with the values just sampled
      // TO-DO put into sample_f?

      for(int g = 0; g < G; g++){
        arma::uvec h = G_finder(g);
        phi_tau_samp.rows(h).fill(phi_samp(g));
      }

      // Sample factor
      sample_f(f_samp,
               X_demean_samp,
               phi_tau_samp,
               var_samp,
               h_samp,
               finder,
               J);

      // Boosting via ASIS to improve mixing of factor loadings
      resample_phi_g(phi_samp,
                     f_samp,
                     tau2_phi_samp,
                     h_samp,
                     J,
                     G);

      // Sample volatilites of the factor
      double mu = sv_para(0);
      double phi = sv_para(1);
      double sigma = std::sqrt(sv_para(2));

      arma::vec datastand = 2. * arma::log(arma::abs(f_samp));
      stochvol::update_fast_sv(datastand, mu, phi, sigma, h0_samp, h_samp, r, prior_spec, expert);

      sv_para = {mu, phi, std::pow(sigma, 2)};

      arma::vec curr_batch;

      if (mod_type_phi == "double") {
        // Sample a with MH
        if (learn_a_phi){
          if (is_adaptive_phi) {
            curr_batch = batches_phi.col(0);
          }
          a_phi_samp = shrinkTVP::DG_MH_step(a_phi_samp,
                                             a_phi_tuning_par,
                                             lambda2_B_phi_samp,
                                             phi_samp,
                                             beta_a_phi,
                                             alpha_a_phi,
                                             adaptive_phi(0),
                                             curr_batch,
                                             curr_sds_phi(0),
                                             target_rates_phi(0),
                                             max_adapts_phi(0),
                                             batch_nrs_phi(0),
                                             batch_sizes_phi(0),
                                             batch_pos_phi(0));
          if (is_adaptive_phi) {
            batches_phi.col(0) = curr_batch;
          }
        }



        shrinkTVP::DG_sample_local_shrink(tau2_phi_samp,
                                          phi_samp,
                                          lambda2_B_phi_samp,
                                          a_phi_samp);

        if (learn_lambda2_B_phi == true) {
          lambda2_B_phi_samp = shrinkTVP::DG_sample_global_shrink(tau2_phi_samp,
                                                                  a_phi_samp,
                                                                  e1_phi,
                                                                  e2_phi);
        }

      } else if (mod_type_phi == "triple") {

        // Sample a_xi/a_tau using MH
        if (learn_a_phi){

          if (is_adaptive_phi) {
            curr_batch = batches_phi.col(0);
          }
          a_phi_samp = shrinkTVP::TG_MH_step(a_phi_samp,
                                             a_phi_tuning_par,
                                             lambda2_B_phi_samp,
                                             lambda2_til_phi_samp,
                                             phi_samp,
                                             beta_a_phi,
                                             alpha_a_phi,
                                             false,
                                             e2_phi_samp,
                                             c_phi_samp,
                                             a_phi_eq_c_phi,
                                             adaptive_phi(0),
                                             curr_batch,
                                             curr_sds_phi(0),
                                             target_rates_phi(0),
                                             max_adapts_phi(0),
                                             batch_nrs_phi(0),
                                             batch_sizes_phi(0),
                                             batch_pos_phi(0));

          if (is_adaptive_phi) {
            batches_phi.col(0) = curr_batch;
          }
        }

        shrinkTVP::TG_sample_prior_var_til(tau2_til_phi_samp,
                                           phi_samp,
                                           lambda2_til_phi_samp,
                                           lambda2_B_phi_samp,
                                           a_phi_samp,
                                           c_phi_samp);
        shrinkTVP::TG_sample_local_shrink(lambda2_til_phi_samp,
                                          phi_samp,
                                          tau2_til_phi_samp,
                                          lambda2_B_phi_samp,
                                          c_phi_samp,
                                          a_phi_samp);

        if (learn_c_phi && (a_phi_eq_c_phi == false)){
          if (is_adaptive_phi) {
            curr_batch = batches_phi.col(1);
          }
          c_phi_samp = shrinkTVP::TG_MH_step(c_phi_samp,
                                             c_phi_tuning_par,
                                             lambda2_B_phi_samp,
                                             lambda2_til_phi_samp,
                                             phi_samp,
                                             beta_c_phi,
                                             alpha_c_phi,
                                             true,
                                             e2_phi_samp,
                                             a_phi_samp,
                                             a_phi_eq_c_phi,
                                             adaptive_phi(1),
                                             curr_batch,
                                             curr_sds_phi(1),
                                             target_rates_phi(1),
                                             max_adapts_phi(1),
                                             batch_nrs_phi(1),
                                             batch_sizes_phi(1),
                                             batch_pos_phi(1));


          if (is_adaptive_phi) {
            batches_phi.col(1) = curr_batch;
          }

        } else if (a_phi_eq_c_phi == true) {
          c_phi_samp = a_phi_samp;
        }


        if (learn_lambda2_B_phi == true) {
          e2_phi_samp = shrinkTVP::TG_sample_d2(lambda2_B_phi_samp,
                                                a_phi_samp,
                                                c_phi_samp);

          lambda2_B_phi_samp = shrinkTVP::TG_sample_global_shrink(tau2_til_phi_samp,
                                                                  lambda2_til_phi_samp,
                                                                  phi_samp,
                                                                  a_phi_samp,
                                                                  c_phi_samp,
                                                                  e2_phi_samp,
                                                                  false);
        }

        shrinkTVP::calc_xi2_tau2(tau2_phi_samp,
                                 tau2_til_phi_samp,
                                 lambda2_til_phi_samp,
                                 lambda2_B_phi_samp,
                                 c_phi_samp,
                                 a_phi_samp);

        std::for_each(tau2_phi_samp.begin(), tau2_phi_samp.end(), shrinkTVP::res_protector);
      }


      /*
       Once everything to do with the factors and the factor loadings has
       been sampled we update log_lam_samp
       */
      for(int g = 0; g < G; g++){
        arma::uvec h = G_finder(g);
        arma::uvec t =  arma::conv_to<arma::uvec>::from(int_tau.rows(h));
        arma::vec f_cur = f_samp.rows(t);
        f_c_tau_samp.rows(h) = arma::as_scalar(phi_samp.row(g)) * f_cur;
      }

      log_lam_samp += f_c_tau_samp;
    }

    // Sample auxilliary survival times
    // R C code uses p(x) = 1/lambda * exp(-1/lambda  * x) parameterization for exponential dist.
    // todo: do this only for cen_tau == 1;
    double exp_param;
    for (int i = 0; i < n_tau; i++){
      // This stops the parameter of the exponential distribution from becoming Inf
      if ((-log_lam_samp(i)) > 50) {
        exp_param = 5.184706e+21;
      } else {
        exp_param = 1/std::exp(log_lam_samp(i));
      }
      tau_res_samp(i) = R::rexp(exp_param);
    }

    tau_res_samp %= cen_tau;
    tau_samp = tau_res_samp + tau_obs;
    eps_samp = -arma::log(tau_samp) - log_lam_samp;

    // Sample components of auxilliary mixture approximation of Gumbel distribution
    sample_comps(mu_samp,
                 var_samp,
                 eps_samp,
                 weights,
                 sigma2,
                 means,
                 0);




    // Store everything
    if (((1 + iter - nburn) % nthin == 0) && (iter >= nburn)){

      int store_pos = ((1 + iter - nburn) / nthin) - 1;

      theta_sr_save.col(store_pos) = theta_sr_samp;
      beta_mean_save.col(store_pos) = beta_mean_samp;
      beta_c_save.slice(store_pos) = beta_c_samp.t();
      if(grouped){

        // Identification of factor loadings
        double modif = 1;
        // if (phi_samp[0] < 0) {
        //   modif = -1;
        // }

        f_save.slice(store_pos) = modif * f_samp;
        h_save.slice(store_pos) = h_samp;
        phi_save.col(store_pos) = modif * phi_samp;
      }
      if (mod_type != "ridge"){
        xi2_save.col(store_pos) = xi2_samp;
        tau2_save.col(store_pos) = tau2_samp;
      }
      // conditional storing
      if (mod_type == "double") {
        xi2_save.col(store_pos) = xi2_samp;
        tau2_save.col(store_pos) = tau2_samp;
      } else if (mod_type == "triple") {
        xi2_save.col(store_pos) = xi2_til_samp;
        tau2_save.col(store_pos) = tau2_til_samp;
      }

      if (mod_type == "triple") {
        kappa2_save.col(store_pos) = kappa2_til_samp;
        lambda2_save.col(store_pos) = lambda2_til_samp;
      }

      if (learn_kappa2_B && (mod_type != "ridge")) {
        kappa2_B_save(store_pos) = kappa2_B_samp;
      }
      if (learn_lambda2_B && (mod_type != "ridge")) {
        lambda2_B_save(store_pos) = lambda2_B_samp;
      }

      if (learn_a_xi && (mod_type != "ridge")) {
        a_xi_save(store_pos) = a_xi_samp;
      }

      if (learn_a_tau && (mod_type != "ridge")) {
        a_tau_save(store_pos) = a_tau_samp;
      }

      if (learn_c_xi && (mod_type == "triple")) {
        c_xi_save(store_pos) = c_xi_samp;
      }

      if (learn_c_tau && (mod_type == "triple")) {
        c_tau_save(store_pos) = c_tau_samp;
      }
    }

    // Conditionally store MH statistics
    if (learn_a_xi && bool(adaptive(0)) && (batch_pos(0) == (batch_sizes(0) - 2))){
      a_xi_sd_save(batch_nrs(0) - 1) = curr_sds(0);
      a_xi_acc_rate_save(batch_nrs(0) - 1) = arma::accu(batches.col(0))/batch_sizes(0);
    }
    if (learn_a_tau && bool(adaptive(1)) && (batch_pos(1) == (batch_sizes(1) - 2))){
      a_tau_sd_save(batch_nrs(1) - 1) = curr_sds(1);
      a_tau_acc_rate_save(batch_nrs(1) - 1) = arma::accu(batches.col(1))/batch_sizes(1);
    }
    if (learn_c_xi && bool(adaptive(2)) && (batch_pos(2) == (batch_sizes(2) - 2))){
      c_xi_sd_save(batch_nrs(2) - 1) = curr_sds(2);
      c_xi_acc_rate_save(batch_nrs(2) - 1) = arma::accu(batches.col(2))/batch_sizes(2);
    }
    if (learn_c_tau && bool(adaptive(3)) && (batch_pos(3) == (batch_sizes(3) - 2))){
      c_tau_sd_save(batch_nrs(3) - 1) = curr_sds(3);
      c_tau_acc_rate_save(batch_nrs(3) - 1) = arma::accu(batches.col(3))/batch_sizes(3);
    }
    // Random sign switch
    for (int i = 0; i < K; i++){
      if(R::runif(0,1) > 0.5){
        theta_sr_samp(i) = -theta_sr_samp(i);
      }
    }

    // Increment progress bar
    if (arma::any(prog_rep_points == iter)){
      p.increment();
    }

    // Check for user interrupts
    if (iter % 50 == 0){
      Rcpp::checkUserInterrupt();
    }
  }

  return List::create(_["beta"] = beta_c_save,
                      _["theta_sr"] = theta_sr_save.t(),
                      _["beta_mean"] = beta_mean_save.t(),
                      _["f"] = f_save,
                      _["h"] = h_save,
                      _["phi"] = phi_save.t(),
                      _["xi2"] = xi2_save.t(),
                      _["tau2"] = tau2_save.t(),
                      _["kappa2"] = kappa2_save.t(),
                      _["lambda2"] = lambda2_save.t(),
                      _["kappa2_B"] = kappa2_B_save,
                      _["lambda2_B"] = lambda2_B_save,
                      _["a_xi"] = a_xi_save,
                      _["a_tau"] = a_tau_save,
                      _["c_xi"] = c_xi_save,
                      _["c_tau"] = c_tau_save,
                      _["MH_diag"] = List::create(
                        _["a_xi_sds"] = a_xi_sd_save,
                        _["a_xi_acc_rate"] = a_xi_acc_rate_save,
                        _["a_tau_sds"] = a_tau_sd_save,
                        _["a_tau_acc_rate"] = a_tau_acc_rate_save,
                        _["c_xi_sds"] = c_xi_sd_save,
                        _["c_xi_acc_rate"] = c_xi_acc_rate_save,
                        _["c_tau_sds"] = c_tau_sd_save,
                        _["c_tau_acc_rate"] = c_tau_acc_rate_save
                      ));

}
