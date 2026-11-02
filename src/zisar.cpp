// ==========================================================================
//  ZISAR: Zero-Inflated Spatial Autoregressive Model
//  TMB C++ template
//
//  Latent field:   (I - rho*W)*eta = X*beta + u,   u ~ N(0, sigma^2 I)
//  Observation:    y_i ~ pi_i * I(0) + (1-pi_i) * NB(exp(eta_i + o_i), kappa)
//  Inflation:      logit(pi_i) = z_i' * gamma + g_i
//                  restricted model (zero_spatial = 0): g = 0
//                  general model    (zero_spatial = 1):
//                    (I - rho0*W)*g = e0,  e0 ~ N(0, sigma0^2 I)
//
//  o_i is an optional offset (log exposure) on the count mean; it scales
//  unit i only and does not pass through the spatial multiplier.
// ==========================================================================

#define TMB_LIB_INIT R_init_zisar
#include <TMB.hpp>

// log|det(I - rho*W)|: from a sparse Cholesky factor of
// (I - rho*W)'(I - rho*W) (sparse = 1) or from the eigenvalues of W
template<class Type>
Type logdet_IW(Type rho, int sparse,
               Eigen::SparseMatrix<Type> &I_n, Eigen::SparseMatrix<Type> &WpWt,
               Eigen::SparseMatrix<Type> &WtW,
               vector<Type> &eig_re, vector<Type> &eig_im)
{
  Type log_det = Type(0.0);
  if (sparse == 1) {
    Eigen::SparseMatrix<Type> Q = I_n - rho * WpWt + rho * rho * WtW;
    log_det = Type(0.5) * newton::log_determinant(Q);
  } else {
    for (int i = 0; i < eig_re.size(); i++) {
      Type re = eig_re(i);
      Type im = eig_im(i);
      Type mod_sq = (Type(1.0) - rho * re) * (Type(1.0) - rho * re) +
                    (rho * im) * (rho * im);
      log_det += Type(0.5) * log(mod_sq + Type(1e-30));
    }
  }
  return log_det;
}

// Negative log-density of a SAR field: (I - rho*W)*f - m ~ N(0, sigma^2 I)
template<class Type>
Type sar_nll(vector<Type> f, vector<Type> m, Type rho, Type log_sigma,
             Eigen::SparseMatrix<Type> &W, Type log_det)
{
  int n = f.size();
  Type sigma2 = exp(Type(2.0) * log_sigma);
  vector<Type> Wf = W * f;
  vector<Type> r_vec = f - rho * Wf - m;
  Type nll = Type(0.5) * Type(n) * log(Type(2.0) * Type(M_PI));  // needed for a correct logLik/AIC
  nll += Type(n) * log_sigma;
  nll -= log_det;
  nll += Type(0.5) / sigma2 * (r_vec * r_vec).sum();
  return nll;
}

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
  // 1 = general model: SAR random field g in the inflation predictor
  DATA_INTEGER(zero_spatial);

  // ---- PARAMETERS (fixed) ----
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(gamma);
  PARAMETER(rho_raw);
  PARAMETER(log_kappa);
  PARAMETER(log_sigma);
  PARAMETER(rho0_raw);     // general model only (mapped out otherwise)
  PARAMETER(log_sigma0);   // general model only (mapped out otherwise)

  // ---- RANDOM EFFECTS ----
  PARAMETER_VECTOR(eta);
  PARAMETER_VECTOR(g);     // length n in the general model, 0 otherwise

  int n = y.size();

  Type rho   = tanh(rho_raw);
  Type kappa = exp(log_kappa);
  Type sigma = exp(log_sigma);

  ADREPORT(rho);
  ADREPORT(kappa);
  ADREPORT(sigma);

  Type nll = Type(0.0);

  // ---- (1) SAR prior on eta ----
  // log p(eta) = -n/2*log(2*pi*sigma2) + log|det(I-rho*W)|
  //              - 1/(2*sigma2) * ||(I-rho*W)*eta - X*beta||^2
  vector<Type> Xbeta = X * beta;
  Type log_det = logdet_IW(rho, logdet_sparse, I_n, WpWt, WtW,
                           eigenvalues_W_re, eigenvalues_W_im);
  nll += sar_nll(eta, Xbeta, rho, log_sigma, W, log_det);

  // ---- (1b) SAR prior on g (general model) ----
  vector<Type> zlin = Z * gamma;
  if (zero_spatial == 1) {
    Type rho0   = tanh(rho0_raw);
    Type sigma0 = exp(log_sigma0);
    ADREPORT(rho0);
    ADREPORT(sigma0);
    Type log_det0 = logdet_IW(rho0, logdet_sparse, I_n, WpWt, WtW,
                              eigenvalues_W_re, eigenvalues_W_im);
    vector<Type> zeros(n);
    zeros.setZero();
    nll += sar_nll(g, zeros, rho0, log_sigma0, W, log_det0);
    zlin += g;
  }

  // ---- (2) ZINB observation likelihood ----
  // On the log scale: log(pi) = -log(1 + exp(-zlin)),
  // log(1 - pi) = -log(1 + exp(zlin))
  for (int i = 0; i < n; i++) {
    Type mu_i = exp(eta(i) + offset(i));
    Type log_pi   = -logspace_add(Type(0.0), -zlin(i));
    Type log_1mpi = -logspace_add(Type(0.0),  zlin(i));
    Type log_p_nb_zero = kappa * log(kappa / (kappa + mu_i));

    if (y(i) < Type(0.5)) {
      // y_i == 0: P(Y=0) = pi + (1-pi) * (kappa/(kappa+mu))^kappa
      nll -= logspace_add(log_pi, log_1mpi + log_p_nb_zero);
    } else {
      // y_i > 0: P(Y=y) = (1-pi) * NB(y; mu, kappa)
      Type log_nb = lgamma(y(i) + kappa) - lgamma(kappa) - lgamma(y(i) + Type(1.0))
                  + log_p_nb_zero
                  + y(i) * log(mu_i / (kappa + mu_i));
      nll -= log_1mpi + log_nb;
    }
  }

  return nll;
}
