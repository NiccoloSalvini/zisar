// ==========================================================================
//  ZISAR: Zero-Inflated Spatial Autoregressive Model
//  TMB C++ template
//
//  Latent field:   (I - rho*W)*eta = X*beta + u,   u ~ N(0, sigma^2 I)
//  Observation:    y_i ~ pi_i * I(0) + (1-pi_i) * NB(exp(eta_i + o_i), kappa)
//  Inflation:      zero_spatial = 0 (restricted): logit(pi) = Z*gamma
//                  zero_spatial = 1 (SAR field):  logit(pi) = Z*gamma + g,
//                    (I - rho0*W)*g = e0,  e0 ~ N(0, sigma0^2 I)
//                  zero_spatial = 2 (SAR on the inflation predictor):
//                    logit(pi) = zeta,  (I - rho0*W)*zeta = Z*gamma
//                  zero_spatial = 3 (reduced-rank SAR field):
//                    logit(pi) = Z*gamma + B*alpha, B an n x K sparse basis
//                    of local functions on a knot lattice,
//                    (I - rho0*Wk)*alpha = e0,  e0 ~ N(0, sigma0^2 I_K)
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

// Least-squares objective 0.5 * ||(I - rho0*W) z - m||^2, whose minimiser is
// the SAR filter z = (I - rho0*W)^{-1} m; solved by newton::Newton, which puts
// the solution on the tape with derivatives from the implicit function theorem
template<class Type>
struct sar_filter_ls {
  Eigen::SparseMatrix<Type> W;
  vector<Type> m;
  Type rho0;
  sar_filter_ls(const Eigen::SparseMatrix<Type> &W, const vector<Type> &m, Type rho0) :
    W(W), m(m), rho0(rho0) {}
  template<class Vec>
  typename Vec::value_type operator()(const Vec &z) {
    typedef typename Vec::value_type T;
    Eigen::SparseMatrix<T> Wt = W.template cast<T>();
    Eigen::Matrix<T, Eigen::Dynamic, 1> zm = z.matrix();
    Eigen::Matrix<T, Eigen::Dynamic, 1> Wz = Wt * zm;
    T r2 = T(0.0);
    for (int i = 0; i < zm.size(); i++) {
      T r = zm(i) - T(rho0) * Wz(i) - T(m(i));
      r2 += r * r;
    }
    return T(0.5) * r2;
  }
};

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
  // Inflation predictor: 0 = restricted, 1 = SAR random field g,
  // 2 = SAR filter on Z*gamma, 3 = reduced-rank SAR field B*alpha
  DATA_INTEGER(zero_spatial);
  // Reduced-rank field (zero_spatial = 3): basis, knot weights and their
  // eigenvalues (dummies otherwise)
  DATA_SPARSE_MATRIX(B);
  DATA_SPARSE_MATRIX(Wk);
  DATA_VECTOR(eig_k_re);
  DATA_VECTOR(eig_k_im);

  // ---- PARAMETERS (fixed) ----
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(gamma);
  PARAMETER(rho_raw);
  PARAMETER(log_kappa);
  PARAMETER(log_sigma);
  PARAMETER(rho0_raw);     // zero_spatial 1, 2 and 3 (mapped out otherwise)
  PARAMETER(log_sigma0);   // zero_spatial 1 and 3 (mapped out otherwise)

  // ---- RANDOM EFFECTS ----
  PARAMETER_VECTOR(eta);
  PARAMETER_VECTOR(g);     // length n in the general model, 0 otherwise
  PARAMETER_VECTOR(alpha); // length K in the reduced-rank model, 0 otherwise

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

  // ---- (1b') reduced-rank SAR field on the knot lattice ----
  if (zero_spatial == 3) {
    Type rho0   = tanh(rho0_raw);
    Type sigma0 = exp(log_sigma0);
    ADREPORT(rho0);
    ADREPORT(sigma0);
    Eigen::SparseMatrix<Type> dummy;
    Type log_det0 = logdet_IW(rho0, 0, dummy, dummy, dummy, eig_k_re, eig_k_im);
    vector<Type> zeros(alpha.size());
    zeros.setZero();
    nll += sar_nll(alpha, zeros, rho0, log_sigma0, Wk, log_det0);
    vector<Type> Ba = B * alpha;
    zlin += Ba;
  }

  // ---- (1c) SAR filter on the inflation predictor ----
  if (zero_spatial == 2) {
    Type rho0 = tanh(rho0_raw);
    ADREPORT(rho0);
    sar_filter_ls<Type> F(W, zlin, rho0);
    newton::newton_config cfg;
    cfg.sparse = true;
    vector<Type> start = zlin;
    vector<Type> zeta = newton::Newton(F, start, cfg);
    zlin = zeta;
    REPORT(zlin);
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
