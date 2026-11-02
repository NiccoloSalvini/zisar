// ==========================================================================
//  ZISAR: Zero-Inflated Spatial Autoregressive Model
//  TMB C++ template
//
//  Latent field:   (I - rho*W)*eta = X*beta + u,   u ~ N(0, sigma^2 I)
//  Observation:    y_i ~ pi_i * I(0) + (1-pi_i) * NB(exp(eta_i + o_i), kappa)
//  Inflation:      logit(pi_i) = z_i' * gamma   [restricted model]
//
//  o_i is an optional offset (log exposure) on the count mean; it scales
//  unit i only and does not pass through the spatial multiplier.
// ==========================================================================

#define TMB_LIB_INIT R_init_zisar
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // ---- DATA ----
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_MATRIX(Z);
  DATA_VECTOR(offset);
  DATA_SPARSE_MATRIX(W);
  DATA_VECTOR(eigenvalues_W_re);
  DATA_VECTOR(eigenvalues_W_im);
  // Large n: log-determinant from a sparse Cholesky factor instead of the
  // eigenvalues (logdet_sparse = 1). I_n, WpWt = W + W', WtW = W'W.
  DATA_INTEGER(logdet_sparse);
  DATA_SPARSE_MATRIX(I_n);
  DATA_SPARSE_MATRIX(WpWt);
  DATA_SPARSE_MATRIX(WtW);

  // ---- PARAMETERS (fixed) ----
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(gamma);
  PARAMETER(rho_raw);
  PARAMETER(log_kappa);
  PARAMETER(log_sigma);

  // ---- RANDOM EFFECTS ----
  PARAMETER_VECTOR(eta);

  int n = y.size();

  Type rho   = tanh(rho_raw);
  Type kappa = exp(log_kappa);
  Type sigma = exp(log_sigma);
  Type sigma2 = sigma * sigma;

  ADREPORT(rho);
  ADREPORT(kappa);
  ADREPORT(sigma);

  Type nll = Type(0.0);

  // ---- (1) SAR prior on eta ----
  // log p(eta) = -n/2*log(2*pi*sigma2) + log|det(I-rho*W)|
  //              - 1/(2*sigma2) * ||(I-rho*W)*eta - X*beta||^2

  // Log-determinant
  Type log_det = Type(0.0);
  if (logdet_sparse == 1) {
    // log|det(I - rho*W)| = 0.5 * log det[(I - rho*W)'(I - rho*W)]
    Eigen::SparseMatrix<Type> Q = I_n - rho * WpWt + rho * rho * WtW;
    log_det = Type(0.5) * newton::log_determinant(Q);
  } else {
    for (int i = 0; i < n; i++) {
      Type re = eigenvalues_W_re(i);
      Type im = eigenvalues_W_im(i);
      Type mod_sq = (Type(1.0) - rho * re) * (Type(1.0) - rho * re) +
                    (rho * im) * (rho * im);
      log_det += Type(0.5) * log(mod_sq + Type(1e-30));
    }
  }

  // Residual: r = (I - rho*W)*eta - X*beta
  vector<Type> Xbeta = X * beta;
  vector<Type> W_eta = W * eta;
  vector<Type> r_vec = eta - rho * W_eta - Xbeta;

  // Negative log-prior, including n/2 * log(2*pi) for a correct logLik/AIC
  nll += Type(0.5) * Type(n) * log(Type(2.0) * Type(M_PI));
  nll += Type(n) * log_sigma;
  nll -= log_det;
  nll += Type(0.5) / sigma2 * (r_vec * r_vec).sum();

  // ---- (2) ZINB observation likelihood ----
  vector<Type> Zgamma = Z * gamma;

  for (int i = 0; i < n; i++) {
    Type mu_i = exp(eta(i) + offset(i));
    Type pi_i = invlogit(Zgamma(i));

    if (y(i) < Type(0.5)) {
      // y_i == 0: P(Y=0) = pi + (1-pi) * (kappa/(kappa+mu))^kappa
      Type log_p_nb_zero = kappa * log(kappa / (kappa + mu_i));
      Type p_nb_zero = exp(log_p_nb_zero);
      Type p_zero = pi_i + (Type(1.0) - pi_i) * p_nb_zero;
      nll -= log(p_zero + Type(1e-20));
    } else {
      // y_i > 0: P(Y=y) = (1-pi) * NB(y; mu, kappa)
      Type log_nb = lgamma(y(i) + kappa) - lgamma(kappa) - lgamma(y(i) + Type(1.0))
                  + kappa * log(kappa / (kappa + mu_i))
                  + y(i) * log(mu_i / (kappa + mu_i));
      nll -= log(Type(1.0) - pi_i + Type(1e-20));
      nll -= log_nb;
    }
  }

  return nll;
}
