#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <string>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// ── Tolerances and macros ────────────────────────────────────────────────────
#define tolr        1e-5
#define ZTOL        (DBL_EPSILON * 10.0)

#ifndef GSL_MAX_DBL
#define GSL_MAX_DBL(a, b) (std::max((double)(a), (double)(b)))
#endif
#ifndef GSL_MIN_INT
#define GSL_MIN_INT(a, b) (std::min((int)(a), (int)(b)))
#endif

static inline double gsl_pow_2(double x) { return x * x; }

static inline double gsl_sf_lnpoch(double a, double x) {
   return std::lgamma(a + x) - std::lgamma(a);
}

static inline int gsl_fcmp(double a, double b, double eps) {
   double diff = a - b;
   double m = std::max(std::fabs(a), std::fabs(b));
   if (std::fabs(diff) <= eps * m) return 0;
   return (diff > 0) ? 1 : -1;
}

// ============================================================================
//  GIG SAMPLER
// ============================================================================

double _gig_mode(double lambda, double omega) {
   if (lambda >= 1.)
      return (sqrt((lambda-1.)*(lambda-1.) + omega*omega) + (lambda-1.)) / omega;
   else
      return omega / (sqrt((1.-lambda)*(1.-lambda) + omega*omega) + (1.-lambda));
}

void _rgig_ROU_shift_alt(double *res, int n, double lambda, double lambda_old,
                         double omega, double alpha) {
   double xm, nc, s, t, U, V, X;
   double a, b, c, p, q, fi, fak, y1, y2, uplus, uminus;
   int i;
   t  = 0.5  * (lambda - 1.);
   s  = 0.25 * omega;
   xm = _gig_mode(lambda, omega);
   nc = t * log(xm) - s * (xm + 1./xm);
   a  = -(2.*(lambda+1.)/omega + xm);
   b  = (2.*(lambda-1.)*xm/omega - 1.);
   c  = xm;
   p  = b - a*a/3.;
   q  = (2.*a*a*a)/27. - (a*b)/3. + c;
   fi = acos(-q / (2.*sqrt(-(p*p*p)/27.)));
   fak = 2.*sqrt(-p/3.);
   y1  = fak * cos(fi/3.)           - a/3.;
   y2  = fak * cos(fi/3. + 4./3.*M_PI) - a/3.;
   uplus  = (y1-xm) * exp(t*log(y1) - s*(y1+1./y1) - nc);
   uminus = (y2-xm) * exp(t*log(y2) - s*(y2+1./y2) - nc);
   for (i = 0; i < n; i++) {
      do {
         U = uminus + unif_rand() * (uplus - uminus);
         V = unif_rand();
         X = U/V + xm;
      } while ((X <= 0.) || ((log(V)) > (t*log(X) - s*(X+1./X) - nc)));
      res[i] = (lambda_old < 0.) ? (alpha/X) : (alpha*X);
   }
}

void _rgig_ROU_noshift(double *res, int n, double lambda, double lambda_old,
                       double omega, double alpha) {
   double xm, nc, ym, um, s, t, U, V, X;
   int i;
   t  = 0.5  * (lambda - 1.);
   s  = 0.25 * omega;
   xm = _gig_mode(lambda, omega);
   nc = t * log(xm) - s * (xm + 1./xm);
   ym = ((lambda+1.) + sqrt((lambda+1.)*(lambda+1.) + omega*omega)) / omega;
   um = exp(0.5*(lambda+1.)*log(ym) - s*(ym+1./ym) - nc);
   for (i = 0; i < n; i++) {
      do {
         U = um * unif_rand();
         V = unif_rand();
         X = U/V;
      } while ((log(V)) > (t*log(X) - s*(X+1./X) - nc));
      res[i] = (lambda_old < 0.) ? (alpha/X) : (alpha*X);
   }
}

void _rgig_newapproach1(double *res, int n, double lambda, double lambda_old,
                        double omega, double alpha) {
   double A[3], Atot, k0, k1, k2, xm, x0, a, U, V, X, hx;
   int i;
   if (lambda >= 1. || omega > 1.) Rcpp::stop("invalid parameters");
   xm = _gig_mode(lambda, omega);
   x0 = omega / (1. - lambda);
   k0 = exp((lambda-1.)*log(xm) - 0.5*omega*(xm+1./xm));
   A[0] = k0 * x0;
   if (x0 >= 2./omega) {
      k1 = 0.; A[1] = 0.;
      k2 = pow(x0, lambda-1.);
      A[2] = k2 * 2. * exp(-omega*x0/2.) / omega;
   } else {
      k1 = exp(-omega);
      A[1] = (lambda == 0.)
         ? k1 * log(2./(omega*omega))
            : k1 / lambda * (pow(2./omega, lambda) - pow(x0, lambda));
      k2 = pow(2/omega, lambda-1.);
      A[2] = k2 * 2 * exp(-1.) / omega;
   }
   Atot = A[0] + A[1] + A[2];
   for (i = 0; i < n; i++) {
      do {
         V = Atot * unif_rand();
         do {
            if (V <= A[0]) { X = x0*V/A[0]; hx = k0; break; }
            V -= A[0];
            if (V <= A[1]) {
               if (lambda == 0.) { X = omega*exp(exp(omega)*V); hx = k1/X; }
               else { X = pow(pow(x0,lambda)+(lambda/k1*V), 1./lambda); hx = k1*pow(X,lambda-1.); }
               break;
            }
            V -= A[1];
            a = (x0 > 2./omega) ? x0 : 2./omega;
            X = -2./omega * log(exp(-omega/2.*a) - omega/(2.*k2)*V);
            hx = k2 * exp(-omega/2.*X);
            break;
         } while(0);
         U = unif_rand() * hx;
      } while (log(U) > (lambda-1.)*log(X) - omega/2.*(X+1./X));
      res[i] = (lambda_old < 0.) ? (alpha/X) : (alpha*X);
   }
}

double rgig(double lambda, double chi, double psi) {
   double omega, alpha;
   double *res;
   int i, n = 1;
   if (!(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
       (chi < 0. || psi < 0) || (chi == 0. && lambda <= 0.) ||
       (psi == 0. && lambda >= 0.)) {
      printf("lambda=%lf, chi=%lf, psi=%lf\n", lambda, chi, psi);
      Rcpp::stop("invalid parameters for GIG distribution!");
   }
   res = (double *)malloc(n * sizeof(double));
   if (chi < ZTOL) {
      for (i = 0; i < n; i++)
         res[i] = (lambda > 0.) ? R::rgamma(lambda, 2./psi) : 1./R::rgamma(-lambda, 2./psi);
   } else if (psi < ZTOL) {
      for (i = 0; i < n; i++)
         res[i] = (lambda > 0.) ? 1./R::rgamma(lambda, 2./chi) : R::rgamma(-lambda, 2./chi);
   } else {
      double lambda_old = lambda;
      if (lambda < 0.) lambda = -lambda;
      alpha = sqrt(chi/psi);
      omega = sqrt(psi*chi);
      do {
         if (lambda > 2. || omega > 3.)
         { _rgig_ROU_shift_alt(res, n, lambda, lambda_old, omega, alpha); break; }
         if (lambda >= 1.-2.25*omega*omega || omega > 0.2)
         { _rgig_ROU_noshift(res, n, lambda, lambda_old, omega, alpha); break; }
         if (lambda >= 0. && omega > 0.)
         { _rgig_newapproach1(res, n, lambda, lambda_old, omega, alpha); break; }
         Rcpp::stop("parameters must satisfy lambda>=0 and omega>0.");
      } while(0);
   }
   double ret = res[0];
   free(res);
   return ret;
}

double rig(double mu) {
   double y   = randn<double>(); y *= y;
   double mu2 = gsl_pow_2(mu);
   double quad = 4.*mu*y + mu2*gsl_pow_2(y);
   double x_div_mu = (2. + mu*y - sqrt(quad)) / 2.;
   double x = mu * x_div_mu;
   if (x <= 0.) return 1e-5;
   double u = log(randu<double>());
   return (u <= -log1p(x_div_mu)) ? x : mu/x_div_mu;
}

inline double log_marginal_pdf(unsigned n, double as, double bs, double vnorm) {
   return gsl_sf_lnpoch(as, ((double)n)/2.)
   + as*log(bs)
   - (as + ((double)n)/2.) * log(bs + vnorm/2.);
}

double log_sum_exp(const arma::vec &x) {
   double mx = x.max();
   return log(sum(exp(x - mx))) + mx;
}

// ============================================================================
//  OMEGA BLOCK HELPERS (shared by both samplers, unchanged)
// ============================================================================

void update_loading_params(arma::mat &Lambda, arma::mat &T_phi,
                           arma::mat &Plam, const double a, const double b) {
   arma::mat abs_theta = abs(Lambda);
   arma::mat psi_inv   = T_phi / abs_theta;
   psi_inv.transform([](double v){ return rig(v); });
   
   T_phi = abs_theta;
   T_phi.transform([&a,&b](double v){
      double tmp = rgig(a-1., 2.*v, 2.*b);
      return GSL_MAX_DBL(tmp, 9e-6);
   });
   
   if (T_phi.has_nan()) Rcpp::stop("T_phi has nan in update_loading_params");
   
   Plam = psi_inv / square(T_phi);
   if (Plam.has_nan()) Rcpp::stop("Plam has nan in update_loading_params");
}

void update_phi(arma::mat &phimat, const arma::mat &Plam,
                const arma::mat &u,  const arma::mat &v,
                const arma::vec &ssqv) {
   const unsigned p = phimat.n_rows;
   arma::mat v_phi = v * phimat;
   for (unsigned j = 0; j < p; ++j) {
      arma::mat tmp  = v.col(j) * phimat.row(j);
      arma::mat util = u - v_phi + tmp;
      arma::rowvec w = (util.t() * v.col(j)).t();
      arma::rowvec d = 1. / sqrt(Plam.row(j) + ssqv(j));
      phimat.row(j)  = d % (d%w + arma::rowvec(phimat.n_cols, fill::randn));
      v_phi -= (tmp - v.col(j) * phimat.row(j));
   }
}

// ============================================================================
//  B-BLOCK HELPERS -- shared, generic over matrix shape
//
//  update_B_params operates purely on B.n_rows/B.n_cols, so it is reused
//  UNCHANGED for the reduced-rank V matrix (Kx x r) by just calling it with
//  V instead of B. Same for full-rank B (p x Kx). No new shrinkage math.
// ============================================================================

void sample_B(arma::mat &B, const arma::mat &Y, const arma::mat &X,
              const arma::mat &Plam_B, const arma::vec &ps,
              const bool huge = false) {
   const int p = Y.n_cols;
   const int K = X.n_cols;
   const int T = Y.n_rows;
   
   for (int j = 0; j < p; ++j) {
      double sigma_j = 1.0 / sqrt(ps(j));
      arma::vec y_norm = Y.col(j) * sigma_j;
      arma::mat X_norm = X * sigma_j;
      
      arma::vec prior_var = Plam_B.row(j).t();
      
      if (huge) {
         arma::vec u(K); u.imbue(R::norm_rand);
         u %= sqrt(prior_var);
         arma::vec delta(T); delta.imbue(R::norm_rand);
         arma::vec v_vec = X_norm * u + delta;
         arma::mat XPV = X_norm.each_row() % prior_var.t();
         arma::vec w   = arma::solve(XPV*X_norm.t() + arma::eye(T,T), y_norm - v_vec);
         B.row(j)      = (u + XPV.t()*w).t();
      } else {
         arma::mat Omega_post = X_norm.t() * X_norm;
         Omega_post.diag() += 1. / prior_var;
         arma::mat chol_op   = chol(Omega_post, "upper");
         arma::mat chol_inv  = arma::inv(trimatu(chol_op));
         arma::vec mu_post   = chol_inv * chol_inv.t() * (X_norm.t() * y_norm);
         arma::vec z(K); z.imbue(R::norm_rand);
         B.row(j) = (mu_post + chol_inv*z).t();
      }
   }
}

// Generic DL-hyperparameter update. Works on ANY matrix shape: called with
// full B (p x Kx) for the full-rank sampler, and with V (Kx x r) for the
// reduced-rank sampler's DL option. No modification needed either way.
void update_B_params(arma::mat &B, arma::mat &T_phi_B,
                     arma::mat &Plam_B, arma::vec &tau_B,
                     const double a_B, const double b_B) {
   const unsigned p  = B.n_rows;
   const unsigned Kx = B.n_cols;
   
   for (unsigned j = 0; j < p; ++j) {
      arma::rowvec coefs    = B.row(j);
      arma::rowvec coefs_sq = square(coefs);
      arma::rowvec lambda_j = T_phi_B.row(j);
      double xi_j           = tau_B(j);
      arma::rowvec psi_j(Kx);
      
      for (unsigned h = 0; h < Kx; ++h) {
         double lam = rgig(a_B - 0.5,
                           coefs_sq(h) / GSL_MAX_DBL(lambda_j(h), 1e-10),
                           2.0 * xi_j);
         lambda_j(h) = GSL_MAX_DBL(lam, 9e-6);
         
         double ps = rgig(0.5,
                          coefs_sq(h) / GSL_MAX_DBL(lambda_j(h)*lambda_j(h), 1e-20),
                          1.0);
         psi_j(h) = GSL_MAX_DBL(ps, 9e-6);
      }
      
      // xi fixed (DL_plus=FALSE), no update needed
      T_phi_B.row(j) = lambda_j;
      Plam_B.row(j) = psi_j % square(lambda_j);
   }
}

// ============================================================================
//  REDUCED-RANK B = U V^T HELPERS
//
//  U : p  x r   (ridge prior only, s_U2)
//  V : Kx x r   (ridge OR DL prior, selected by prior_B string)
// ============================================================================

// Sample U | V, Y, ps. Row-by-row r-dim Gaussian regression on Z = X V.
void sample_U(arma::mat &U, const arma::mat &Ymat, const arma::mat &Z,
              const arma::vec &ps, const double s_U2) {
   const unsigned p = Ymat.n_cols;
   const unsigned r = Z.n_cols;
   arma::mat ZtZ = symmatu(Z.t() * Z);
   arma::mat ZtY = Z.t() * Ymat;
   for (unsigned j = 0; j < p; ++j) {
      arma::mat M = ps(j) * ZtZ;
      M.diag()   += 1.0 / s_U2;
      arma::mat Rchol = chol(symmatu(M));
      arma::vec rhs   = ps(j) * ZtY.col(j);
      arma::vec mu    = solve(trimatu(Rchol),
                              solve(trimatl(Rchol.t()), rhs));
      arma::vec z(r); z.imbue(R::norm_rand);
      U.row(j) = (mu + solve(trimatu(Rchol), z)).t();
   }
}

// Sample V | U, Y, ps. Joint (Kx*r)-dim Gaussian via the Kronecker form of
// the posterior precision. PriorVarV (Kx x r) supplies the per-entry prior
// variance: constant matrix for ridge, or Plam_V (DL local*global scales)
// for the DL option. Same function serves both priors.
void sample_V(arma::mat &V, const arma::mat &U, const arma::mat &Ymat,
              const arma::mat &X, const arma::mat &XtX,
              const arma::vec &ps, const arma::mat &PriorVarV) {
   const unsigned Kx = X.n_cols;
   const unsigned r  = U.n_cols;
   
   arma::mat Ud   = U.each_col() % ps;
   arma::mat UtDU = symmatu(U.t() * Ud);
   
   arma::mat P = kron(UtDU, XtX);
   P.diag()   += 1.0 / vectorise(PriorVarV);
   
   arma::mat XtY   = X.t() * Ymat;
   arma::mat XtYD  = XtY.each_row() % ps.t();
   arma::mat b_mat = XtYD * U;
   arma::vec b_vec = vectorise(b_mat);
   
   arma::mat Rchol = chol(symmatu(P));
   arma::vec mu    = solve(trimatu(Rchol),
                           solve(trimatl(Rchol.t()), b_vec));
   arma::vec z(Kx*r); z.imbue(R::norm_rand);
   arma::vec vsol  = mu + solve(trimatu(Rchol), z);
   
   V = reshape(vsol, Kx, r);
}

// ============================================================================
//  DP RESIDUAL VARIANCE UPDATE (shared, unchanged)
// ============================================================================

// ── Stirling-gamma conjugate DP-concentration update (Zito et al. 2024) ──────
// Posterior log-density of alpha under a Stirling-gamma(a,b,m) prior, given
// k clusters in n observations (up to a constant).
static inline double log_pdf_Sg_post(double x, double a, double b, double m,
                                     double k, double n) {
   return (a + k - 1.0) * std::log(x)
   - b * (std::lgamma(x + m) - std::lgamma(x))
   -      (std::lgamma(x + n) - std::lgamma(x));
}

// Newton mode-finder for the ratio-of-uniforms bounds. `extra2logx` adds the
// +2log(x) term (for Mv). Returns the maximised log-value; guards keep x>0 and
// cap iterations. dgrad uses R::digamma for the psi terms.
static inline double sg_bound(double a, double b, double m, double k, double n,
                              bool extra2logx) {
   double x = 1.0;
   for (int it = 0; it < 50; ++it) {
      double g = -(a + k - 1.0) / x
      + b * (R::digamma(x + m) - R::digamma(x))
      +     (R::digamma(x + n) - R::digamma(x));
      double h =  (a + k - 1.0) / (x * x)
         + b * (R::trigamma(x + m) - R::trigamma(x))
         +     (R::trigamma(x + n) - R::trigamma(x));
         if (extra2logx) { g += 2.0 / x; h += -2.0 / (x * x); }
         // maximise: Newton on the gradient of (neg log-dens). h is d/dx of -g here.
         if (std::fabs(h) < 1e-12) break;
         double step = g / h;
         double xn = x - step;
         if (xn <= 1e-10 || !std::isfinite(xn)) xn = x / 2.0;
         if (std::fabs(xn - x) < 1e-8) { x = xn; break; }
         x = xn;
   }
   double val = log_pdf_Sg_post(x, a, b, m, k, n);
   if (extra2logx) val += 2.0 * std::log(x);
   return val;
}

// One draw from the Stirling-gamma posterior via ratio-of-uniforms, with a
// max-tries cap. Returns -1.0 on failure so the caller can fall back.
static inline double rSg_posterior_draw(double a, double b, double m,
                                        double k, double n) {
   double Mu = sg_bound(a, b, m, k, n, false);
   double Mv = sg_bound(a, b, m, k, n, true);
   double s  = std::exp(0.5 * (Mv - Mu));
   for (int tries = 0; tries < 100000; ++tries) {
      double u = unif_rand();
      double v = unif_rand();
      if (u <= 0.0) continue;
      double x = s * v / u;
      if (std::log(u) < 0.5 * (log_pdf_Sg_post(x, a, b, m, k, n) - Mu))
         return x;
   }
   return -1.0;  // failed; caller falls back to Escobar-West
}

double update_dir_prec(double alpha, double a, double b, unsigned n, unsigned k) {
   double phi  = R::rbeta(alpha+1., (double)n);
   double gam1 = a + k, gam2 = b - log(phi);
   double pi   = (gam1-1.) / (n*gam2);
   return (log(randu()) < log(pi)-log1p(pi))
      ? randg(distr_param(gam1,   1./gam2))
         : randg(distr_param(gam1-1.,1./gam2));
}

// ============================================================================
//  FULL-RANK SAMPLER (original, unchanged)
// ============================================================================

// Per-observation Gaussian log-likelihood for the current draw.
// Omega = phimat phimat^T + diag(ps); E is the Teff x p residual matrix
// (rows e_t = y_t - B y_{t-1}). Returns a length-Teff vec of log p(y_t | .).
arma::vec obs_loglik(const arma::mat &phimat, const arma::vec &ps,
                     const arma::mat &E) {
   const unsigned p    = phimat.n_rows;
   const unsigned Teff = E.n_rows;
   arma::mat Omega = phimat * phimat.t();
   Omega.diag() += ps;
   arma::mat R;
   if (!arma::chol(R, Omega)) {           // upper Cholesky; R^T R = Omega
      return arma::vec(Teff, fill::value(-1e10));
   }
   double logdet = 2.0 * arma::sum(arma::log(R.diag()));
   arma::mat Z = E * R.t();               // each row z_t, quad_t = z_t . z_t
   arma::vec quad = arma::sum(arma::square(Z), 1);
   double c = 0.5 * logdet - 0.5 * (double)p * std::log(2.0 * M_PI);
   return c - 0.5 * quad;
}

// [[Rcpp::export]]
Rcpp::List var_precision_sampler(
      const arma::mat &Y,
      const int lags    = 1,
      const int k       = 5,
      const int nrun    = 5000,
      const int burn    = 2500,
      const int thin    = 5,
      const double a_L  = 0.5,
      const double b_L  = 0.5,
      const double a_B  = 0.5,
      const double b_B  = 0.5,
      //  DP base measure Gamma(as, rate = bs).  Default moved from 0.1 to 1:
      //  at 0.1 the base measure puts 33% of its mass below 1e-4 and the
      //  idiosyncratic precisions collapse onto the 1e-8 floor.
      const double as   = 1.0,
      const double bs   = 1.0,
      const double a_dir= 0.1,
      const double b_dir= 0.1,
      const std::string dp_prior = "escobar_west",
      const double sg_a = 3.0,
      const double sg_b = 1.0,
      const double sg_m = 0.0,
      Rcpp::Nullable<Rcpp::NumericMatrix> Lambda_init = R_NilValue,
      Rcpp::Nullable<Rcpp::NumericVector> delta_init  = R_NilValue
) {
   const bool use_sg = (dp_prior == "stirling_gamma");
   const double sg_m_eff = (sg_m > 0.0) ? sg_m : (double)Y.n_cols;
   const unsigned T = Y.n_rows;
   const unsigned p = Y.n_cols;
   const unsigned Teff = T - lags;
   const unsigned Kx   = lags * p;
   const double a_B_eff = 1.0 / (double)(p * Kx);
   
   arma::mat Ymat(Teff, p), X(Teff, Kx);
   for (unsigned t = 0; t < Teff; ++t) {
      Ymat.row(t) = Y.row(t + lags);
      for (int lag = 0; lag < lags; ++lag)
         X.row(t).subvec(lag*p, (lag+1)*p - 1) = Y.row(t + lags - 1 - lag);
   }
   
   arma::mat B(p, Kx, fill::zeros);
   arma::mat T_phi_B(p, Kx, fill::value(1.0/Kx));
   arma::mat Plam_B(p, Kx, fill::value(1.0/Kx));
   arma::vec tau_B(p, fill::value(0.5));
   
   arma::mat phimat(p, k, fill::randn);
   if (Lambda_init.isNotNull()) {
      arma::mat L0 = Rcpp::as<arma::mat>(Rcpp::NumericMatrix(Lambda_init));
      if (L0.n_rows == p && L0.n_cols == (unsigned)k) phimat = L0;
      else Rcpp::warning("Lambda_init has wrong dimensions; using random init.");
   }
   arma::mat T_phi  = randg<arma::mat>(p, k, distr_param(a_L, 0.5));
   T_phi.transform([](double v){ return GSL_MAX_DBL(v, 9e-17); });
   arma::mat Plam   = 1. / (1. + square(T_phi));
   
   const unsigned nmix = p;
   double dir_prec = 1.0;
   double prob_empty = log(dir_prec);
   arma::uvec nj_val(nmix, fill::zeros), del;
   arma::vec  log_nj_val(nmix), probs(nmix), log_probs(nmix);
   arma::Col<int> d;
   arma::uvec *inds_eq_j = new arma::uvec[nmix];
   arma::vec ps(p, fill::ones);
   if (delta_init.isNotNull()) {
      arma::vec d0 = Rcpp::as<arma::vec>(Rcpp::NumericVector(delta_init));
      if (d0.n_elem == p) ps = arma::clamp(d0, 1e-8, arma::datum::inf);
      else Rcpp::warning("delta_init has wrong length; using default init.");
   }
   
   const unsigned n_mc = (unsigned)std::floor((double)nrun / thin);
   arma::cube  Lambda_mc(p, k, n_mc, fill::zeros);
   arma::mat   ps_mat(n_mc, p, fill::zeros);
   arma::cube  B_mc(p, Kx, n_mc, fill::zeros);
   arma::mat   loglik_mc(Teff, n_mc, fill::zeros);
   arma::vec   dir_prec_mc(n_mc, fill::zeros);
   arma::vec   nclust_mc(n_mc, fill::zeros);
   
   //  Idiosyncratic-precision floor.  Without it the DP atom for a series whose
   //  residual sum of squares is dominated by one or two quarters can collapse
   //  towards zero, Lambda' Delta^-1 Lambda + I loses positive-definiteness and
   //  chol() aborts the sampler.  1e-8 is the floor delta_init already applies.
   //  Clamping draws no random numbers, so the chain is unchanged up to the
   //  first iteration at which the floor binds.  floor_hits counts, per series,
   //  the iterations at which it bound and floor_first the first such iteration
   //  (0 = never); both are returned.
   const double PS_FLOOR = 1e-8;
   arma::uvec floor_hits(p,  fill::zeros);
   arma::uvec floor_first(p, fill::zeros);
   //  floor_iters counts the ITERATIONS in which at least one series was
   //  floored.  It is a plain counter over a flag already computed by the loop
   //  below: it draws no random numbers and alters no state the sampler reads,
   //  so every draw is identical to a build without it.
   arma::uword floor_iters = 0;

   Progress prog(nrun, true);
   
   for (int i = 0; i < nrun; ++i) {
      
      if (Progress::check_abort()) { delete[] inds_eq_j; return R_NilValue; }
      prog.increment();
      
      bool huge = (Kx > (int)Teff);
      sample_B(B, Ymat, X, Plam_B, ps, huge);
      update_B_params(B, T_phi_B, Plam_B, tau_B, a_B_eff, b_B);
      
      arma::mat E = Ymat - X * B.t();
      
      arma::mat Lmsg  = phimat.each_col() / ps;
      arma::mat Lmsg2 = symmatu(Lmsg.t() * phimat);
      Lmsg2.diag()   += 1.0;
      
      arma::mat noise(Teff, k, fill::randn);
      arma::mat chol_prec_u = chol(Lmsg2);
      arma::mat u = noise * trimatu(chol_prec_u);
      
      arma::mat tmpmt = solve(trimatl(chol_prec_u.t()), Lmsg.t());
      arma::mat v     = E + noise * tmpmt;
      arma::vec ssqv  = sum(square(v), 0).t();
      
      update_loading_params(phimat, T_phi, Plam, a_L, b_L);
      update_phi(phimat, Plam, u, v, ssqv);
      
      if (!i) {
         unsigned n_init = GSL_MIN_INT((unsigned)std::ceil(log((double)p)), 2u);
         arma::mat dat(1, p); dat.row(0) = ssqv.t();
         del.set_size(p);
         for (unsigned j = 0; j < p; ++j) del(j) = j % n_init;
         for (unsigned j = 0; j < nmix; ++j) {
            inds_eq_j[j] = find(del == j);
            nj_val(j)    = inds_eq_j[j].n_elem;
            log_nj_val(j)= (nj_val(j) > 0) ? log((double)nj_val(j)) : -1e10;
         }
      }
      
      arma::vec v_atoms(nmix, fill::zeros);
      for (unsigned j = 0; j < nmix; ++j)
         if (nj_val(j)) v_atoms(j) = sum(ssqv(inds_eq_j[j]));
         
         double log_probs_max;
         for (unsigned jj = 0; jj < p; ++jj) {
            double dens, cluster_prob;
            double pdf_empty = log_marginal_pdf(Teff, as, bs, ssqv(jj));
            
            arma::uvec non_empty = find(nj_val > 0);
            arma::uword cur      = del(jj);
            
            if (nj_val(cur) > 1) {
               arma::uvec empty_inds = find(nj_val == 0, 1, "first");
               non_empty.resize(non_empty.n_elem + 1);
               non_empty.tail(1) = empty_inds(0);
            }
            
            unsigned nclust = non_empty.n_elem;
            log_probs.set_size(nclust);
            
            for (unsigned j1 = 0; j1 < nclust; ++j1) {
               unsigned j = non_empty(j1);
               if (j != cur) {
                  if (nj_val(j) < 1) {
                     dens = pdf_empty; cluster_prob = prob_empty;
                  } else {
                     dens = log_marginal_pdf(Teff,
                                             as + ((double)(Teff*nj_val(j)))/2.,
                                             bs + v_atoms(j)/2., ssqv(jj));
                     cluster_prob = log_nj_val(j);
                  }
                  log_probs(j1) = dens + cluster_prob;
               }
            }
            unsigned j = cur;
            if (nj_val(j) == 1) {
               dens = pdf_empty; cluster_prob = prob_empty;
            } else {
               dens = log_marginal_pdf(Teff,
                                       as + ((double)(Teff*(nj_val(j)-1)))/2.,
                                       bs + (v_atoms(j)-ssqv(jj))/2., ssqv(jj));
               cluster_prob = log((double)(nj_val(j)-1));
            }
            arma::uvec indddd = find(non_empty == j, 1, "first");
            log_probs(indddd(0)) = dens + cluster_prob;
            
            double log_DEN = log_sum_exp(log_probs);
            probs = normalise(exp(log_probs - log_DEN), 1);
            
            if (gsl_fcmp(sum(probs), 1.0, 1e-5)) {
               log_probs_max = max(log_probs);
               log_probs    -= log_probs_max;
               probs         = normalise(exp(log_probs), 1);
            }
            if (sum(probs) == 0) Rcpp::stop("sum(probs)==0 in DP update");
            
            d.set_size(probs.n_elem);
            R::rmultinom(1, probs.begin(), probs.n_elem, d.begin());
            arma::uvec dd = find(d == 1, 1, "first");
            del(jj) = non_empty(dd(0));
            
            if (del(jj) != cur) {
               j = cur; --nj_val(j);
               if (nj_val(j)) { log_nj_val(j) = log((double)nj_val(j)); v_atoms(j) -= ssqv(jj); }
               j = del(jj); log_nj_val(j) = log((double)(++nj_val(j)));
               if (nj_val(j) == 1) v_atoms(j) = ssqv(jj); else v_atoms(j) += ssqv(jj);
            }
         }
         
         for (unsigned j = 0; j < nmix; ++j) {
            if (nj_val(j)) {
               inds_eq_j[j] = find(del == j);
               double atom = randg(distr_param(as + 0.5*(Teff*nj_val(j)),
                                               1./(bs + 0.5*v_atoms(j))));
               ps(inds_eq_j[j]).fill(atom);
            } else inds_eq_j[j].reset();
         }

         bool any_floored = false;
         for (unsigned j = 0; j < p; ++j)
            if (ps(j) < PS_FLOOR) {
               ++floor_hits(j);
               if (floor_first(j) == 0) floor_first(j) = i + 1;
               ps(j) = PS_FLOOR;
               any_floored = true;
            }
         if (any_floored) ++floor_iters;
         {
            unsigned n_clust = (unsigned)sum(nj_val > 0);
            double a_sg = std::max((unsigned)1, n_clust);
            if (use_sg) {
               double draw = rSg_posterior_draw(sg_a, sg_b, sg_m_eff,
                                                (double)a_sg, (double)p);
               dir_prec = (draw > 0.0) ? draw
               : update_dir_prec(dir_prec, a_dir, b_dir, p, n_clust);
            } else {
               dir_prec = update_dir_prec(dir_prec, a_dir, b_dir, p, n_clust);
            }
         }
         prob_empty = log(dir_prec);
         
         int rem = (i+1) % thin;
         if (rem == 0) {
            int q = (int)std::floor((double)(i+1)/thin) - 1;
            Lambda_mc.slice(q) = phimat;
            ps_mat.row(q)      = ps.t();
            B_mc.slice(q)      = B;
            loglik_mc.col(q)   = obs_loglik(phimat, ps, E);
            dir_prec_mc(q)     = dir_prec;
            nclust_mc(q)       = (double)sum(nj_val > 0);
         }
   }
   
   delete[] inds_eq_j;
   
   return Rcpp::List::create(
      Rcpp::Named("Lambda")    = Lambda_mc,
      Rcpp::Named("residuals") = ps_mat,
      Rcpp::Named("B")         = B_mc,
      Rcpp::Named("loglik")    = loglik_mc,
      Rcpp::Named("dir_prec")  = dir_prec_mc,
      Rcpp::Named("nclust")    = nclust_mc,
      Rcpp::Named("floor_hits")  = floor_hits,
      Rcpp::Named("floor_first") = floor_first,
      Rcpp::Named("floor_iters") = floor_iters
   );
}

// ============================================================================
//  REDUCED-RANK SAMPLER (B = U V^T), with selectable prior on V
//     prior_B == "ridge": vec(V) ~ N(0, s_V2 * I)   [default, fast, simple]
//     prior_B == "dl":    Dirichlet-Laplace shrinkage on V, same machinery
//                         as the full-rank sampler's DL prior on B, just
//                         applied to the Kx x r matrix V instead of p x Kx.
//  U always uses a ridge prior (s_U2); r is small so this rarely matters.
// ============================================================================

// [[Rcpp::export]]
Rcpp::List var_precision_sampler_rr(
      const arma::mat &Y,
      const int lags       = 1,
      const int k          = 5,
      const int rank       = 3,
      const int nrun       = 5000,
      const int burn       = 2500,
      const int thin       = 5,
      const double a_L     = 0.5,
      const double b_L     = 0.5,
      const std::string prior_B = "ridge",   // "ridge" or "dl"
      const double s_U2    = 1.0,            // ridge prior variance for U
      const double s_V2    = 1.0,            // ridge prior variance for V (if prior_B=="ridge")
      const double a_V     = 0.5,            // DL shape for V (if prior_B=="dl")
      const double b_V     = 0.5,            // DL rate  for V (if prior_B=="dl")
      //  DP base measure Gamma(as, rate = bs); default 1, as in the full-rank
      //  sampler above.
      const double as      = 1.0,
      const double bs      = 1.0,
      const double a_dir   = 0.1,
      const double b_dir   = 0.1,
      const std::string dp_prior = "escobar_west",
      const double sg_a = 3.0,
      const double sg_b = 1.0,
      const double sg_m = 0.0,
      Rcpp::Nullable<Rcpp::NumericMatrix> Lambda_init = R_NilValue,
      Rcpp::Nullable<Rcpp::NumericVector> delta_init  = R_NilValue,
      Rcpp::Nullable<Rcpp::NumericMatrix> V_init      = R_NilValue
) {
   const bool use_sg = (dp_prior == "stirling_gamma");
   const double sg_m_eff = (sg_m > 0.0) ? sg_m : (double)Y.n_cols;
   const bool use_dl  = (prior_B == "dl");
   if (!use_dl && prior_B != "ridge")
      Rcpp::stop("prior_B must be 'ridge' or 'dl'");
   
   const unsigned T = Y.n_rows;
   const unsigned p = Y.n_cols;
   const unsigned Teff = T - lags;
   const unsigned Kx   = lags * p;
   const unsigned r    = (unsigned)rank;
   
   if (r > p || r > Kx)
      Rcpp::stop("rank must be <= min(p, lags*p)");
   
   arma::mat Ymat(Teff, p), X(Teff, Kx);
   for (unsigned t = 0; t < Teff; ++t) {
      Ymat.row(t) = Y.row(t + lags);
      for (int lag = 0; lag < lags; ++lag)
         X.row(t).subvec(lag*p, (lag+1)*p - 1) = Y.row(t + lags - 1 - lag);
   }
   arma::mat XtX = symmatu(X.t() * X);
   
   arma::mat U(p,  r, fill::randn);  U *= 0.1;
   arma::mat V(Kx, r, fill::randn);  V *= 0.1;
   if (V_init.isNotNull()) {
      arma::mat V0 = Rcpp::as<arma::mat>(Rcpp::NumericMatrix(V_init));
      if (V0.n_rows == Kx && V0.n_cols == r) V = V0;
      else Rcpp::warning("V_init has wrong dimensions; using random init.");
   }
   
   // -- DL state for V, only used if prior_B == "dl" --------------------------
   const double a_V_eff = 1.0 / (double)(Kx * r);   // analogous to a_B_eff
   arma::mat T_phi_V(Kx, r, fill::value(1.0/r));
   arma::mat Plam_V(Kx, r, fill::value(1.0/r));
   arma::vec tau_V(Kx, fill::value(0.5));
   
   
   // constant ridge prior-variance matrix, reused every iteration if ridge
   arma::mat ridge_mat_V(Kx, r, fill::value(s_V2));
   
   arma::mat phimat(p, k, fill::randn);
   if (Lambda_init.isNotNull()) {
      arma::mat L0 = Rcpp::as<arma::mat>(Rcpp::NumericMatrix(Lambda_init));
      if (L0.n_rows == p && L0.n_cols == (unsigned)k) phimat = L0;
      else Rcpp::warning("Lambda_init has wrong dimensions; using random init.");
   }
   arma::mat T_phi  = randg<arma::mat>(p, k, distr_param(a_L, 0.5));
   T_phi.transform([](double v){ return GSL_MAX_DBL(v, 9e-17); });
   arma::mat Plam   = 1. / (1. + square(T_phi));
   
   const unsigned nmix = p;
   double dir_prec = 1.0;
   double prob_empty = log(dir_prec);
   arma::uvec nj_val(nmix, fill::zeros), del;
   arma::vec  log_nj_val(nmix), probs(nmix), log_probs(nmix);
   arma::Col<int> d;
   arma::uvec *inds_eq_j = new arma::uvec[nmix];
   arma::vec ps(p, fill::ones);
   if (delta_init.isNotNull()) {
      arma::vec d0 = Rcpp::as<arma::vec>(Rcpp::NumericVector(delta_init));
      if (d0.n_elem == p) ps = arma::clamp(d0, 1e-8, arma::datum::inf);
      else Rcpp::warning("delta_init has wrong length; using default init.");
   }
   
   const unsigned n_mc = (unsigned)std::floor((double)nrun / thin);
   arma::cube  Lambda_mc(p, k, n_mc, fill::zeros);
   arma::mat   ps_mat(n_mc, p, fill::zeros);
   arma::cube  B_mc(p, Kx, n_mc, fill::zeros);
   arma::cube  U_mc(p,  r, n_mc, fill::zeros);
   arma::cube  V_mc(Kx, r, n_mc, fill::zeros);
   arma::mat   loglik_mc(Teff, n_mc, fill::zeros);
   arma::vec   dir_prec_mc(n_mc, fill::zeros);
   arma::vec   nclust_mc(n_mc, fill::zeros);
   
   //  Idiosyncratic-precision floor.  Without it the DP atom for a series whose
   //  residual sum of squares is dominated by one or two quarters can collapse
   //  towards zero, Lambda' Delta^-1 Lambda + I loses positive-definiteness and
   //  chol() aborts the sampler.  1e-8 is the floor delta_init already applies.
   //  Clamping draws no random numbers, so the chain is unchanged up to the
   //  first iteration at which the floor binds.  floor_hits counts, per series,
   //  the iterations at which it bound and floor_first the first such iteration
   //  (0 = never); both are returned.
   const double PS_FLOOR = 1e-8;
   arma::uvec floor_hits(p,  fill::zeros);
   arma::uvec floor_first(p, fill::zeros);
   //  floor_iters counts the ITERATIONS in which at least one series was
   //  floored.  It is a plain counter over a flag already computed by the loop
   //  below: it draws no random numbers and alters no state the sampler reads,
   //  so every draw is identical to a build without it.
   arma::uword floor_iters = 0;

   Progress prog(nrun, true);
   
   for (int i = 0; i < nrun; ++i) {
      
      if (Progress::check_abort()) { delete[] inds_eq_j; return R_NilValue; }
      prog.increment();
      
      // ── BLOCK 1: reduced-rank B via U, V ─────────────────────────────────
      arma::mat Z = X * V;
      sample_U(U, Ymat, Z, ps, s_U2);
      
      const arma::mat &PriorVarV = use_dl ? Plam_V : ridge_mat_V;
      sample_V(V, U, Ymat, X, XtX, ps, PriorVarV);
      if (use_dl)  update_B_params(V, T_phi_V, Plam_V, tau_V, a_V_eff, b_V);
      
      Z = X * V;
      arma::mat B = U * V.t();
      
      arma::mat E = Ymat - Z * U.t();
      
      // ── BLOCK 2: Sample Omega | residuals E ───────────────────────────────
      arma::mat Lmsg  = phimat.each_col() / ps;
      arma::mat Lmsg2 = symmatu(Lmsg.t() * phimat);
      Lmsg2.diag()   += 1.0;
      
      arma::mat noise(Teff, k, fill::randn);
      arma::mat chol_prec_u = chol(Lmsg2);
      arma::mat u = noise * trimatu(chol_prec_u);
      
      arma::mat tmpmt = solve(trimatl(chol_prec_u.t()), Lmsg.t());
      arma::mat v     = E + noise * tmpmt;
      arma::vec ssqv  = sum(square(v), 0).t();
      
      update_loading_params(phimat, T_phi, Plam, a_L, b_L);
      update_phi(phimat, Plam, u, v, ssqv);
      
      // ── BLOCK 3: Sample residual variances (DP) ───────────────────────────
      if (!i) {
         unsigned n_init = GSL_MIN_INT((unsigned)std::ceil(log((double)p)), 2u);
         arma::mat dat(1, p); dat.row(0) = ssqv.t();
         del.set_size(p);
         for (unsigned j = 0; j < p; ++j) del(j) = j % n_init;
         for (unsigned j = 0; j < nmix; ++j) {
            inds_eq_j[j] = find(del == j);
            nj_val(j)    = inds_eq_j[j].n_elem;
            log_nj_val(j)= (nj_val(j) > 0) ? log((double)nj_val(j)) : -1e10;
         }
      }
      
      arma::vec v_atoms(nmix, fill::zeros);
      for (unsigned j = 0; j < nmix; ++j)
         if (nj_val(j)) v_atoms(j) = sum(ssqv(inds_eq_j[j]));
         
         double log_probs_max;
         for (unsigned jj = 0; jj < p; ++jj) {
            double dens, cluster_prob;
            double pdf_empty = log_marginal_pdf(Teff, as, bs, ssqv(jj));
            
            arma::uvec non_empty = find(nj_val > 0);
            arma::uword cur      = del(jj);
            
            if (nj_val(cur) > 1) {
               arma::uvec empty_inds = find(nj_val == 0, 1, "first");
               non_empty.resize(non_empty.n_elem + 1);
               non_empty.tail(1) = empty_inds(0);
            }
            
            unsigned nclust = non_empty.n_elem;
            log_probs.set_size(nclust);
            
            for (unsigned j1 = 0; j1 < nclust; ++j1) {
               unsigned j = non_empty(j1);
               if (j != cur) {
                  if (nj_val(j) < 1) {
                     dens = pdf_empty; cluster_prob = prob_empty;
                  } else {
                     dens = log_marginal_pdf(Teff,
                                             as + ((double)(Teff*nj_val(j)))/2.,
                                             bs + v_atoms(j)/2., ssqv(jj));
                     cluster_prob = log_nj_val(j);
                  }
                  log_probs(j1) = dens + cluster_prob;
               }
            }
            unsigned j = cur;
            if (nj_val(j) == 1) {
               dens = pdf_empty; cluster_prob = prob_empty;
            } else {
               dens = log_marginal_pdf(Teff,
                                       as + ((double)(Teff*(nj_val(j)-1)))/2.,
                                       bs + (v_atoms(j)-ssqv(jj))/2., ssqv(jj));
               cluster_prob = log((double)(nj_val(j)-1));
            }
            arma::uvec indddd = find(non_empty == j, 1, "first");
            log_probs(indddd(0)) = dens + cluster_prob;
            
            double log_DEN = log_sum_exp(log_probs);
            probs = normalise(exp(log_probs - log_DEN), 1);
            
            if (gsl_fcmp(sum(probs), 1.0, 1e-5)) {
               log_probs_max = max(log_probs);
               log_probs    -= log_probs_max;
               probs         = normalise(exp(log_probs), 1);
            }
            if (sum(probs) == 0) Rcpp::stop("sum(probs)==0 in DP update");
            
            d.set_size(probs.n_elem);
            R::rmultinom(1, probs.begin(), probs.n_elem, d.begin());
            arma::uvec dd = find(d == 1, 1, "first");
            del(jj) = non_empty(dd(0));
            
            if (del(jj) != cur) {
               j = cur; --nj_val(j);
               if (nj_val(j)) { log_nj_val(j) = log((double)nj_val(j)); v_atoms(j) -= ssqv(jj); }
               j = del(jj); log_nj_val(j) = log((double)(++nj_val(j)));
               if (nj_val(j) == 1) v_atoms(j) = ssqv(jj); else v_atoms(j) += ssqv(jj);
            }
         }
         
         for (unsigned j = 0; j < nmix; ++j) {
            if (nj_val(j)) {
               inds_eq_j[j] = find(del == j);
               double atom = randg(distr_param(as + 0.5*(Teff*nj_val(j)),
                                               1./(bs + 0.5*v_atoms(j))));
               ps(inds_eq_j[j]).fill(atom);
            } else inds_eq_j[j].reset();
         }

         bool any_floored = false;
         for (unsigned j = 0; j < p; ++j)
            if (ps(j) < PS_FLOOR) {
               ++floor_hits(j);
               if (floor_first(j) == 0) floor_first(j) = i + 1;
               ps(j) = PS_FLOOR;
               any_floored = true;
            }
         if (any_floored) ++floor_iters;
         {
            unsigned n_clust = (unsigned)sum(nj_val > 0);
            double a_sg = std::max((unsigned)1, n_clust);
            if (use_sg) {
               double draw = rSg_posterior_draw(sg_a, sg_b, sg_m_eff,
                                                (double)a_sg, (double)p);
               dir_prec = (draw > 0.0) ? draw
               : update_dir_prec(dir_prec, a_dir, b_dir, p, n_clust);
            } else {
               dir_prec = update_dir_prec(dir_prec, a_dir, b_dir, p, n_clust);
            }
         }
         prob_empty = log(dir_prec);
         
         int rem = (i+1) % thin;
         if (rem == 0) {
            int q = (int)std::floor((double)(i+1)/thin) - 1;
            Lambda_mc.slice(q) = phimat;
            ps_mat.row(q)      = ps.t();
            B_mc.slice(q)      = B;
            U_mc.slice(q)      = U;
            V_mc.slice(q)      = V;
            loglik_mc.col(q)   = obs_loglik(phimat, ps, E);
            dir_prec_mc(q)     = dir_prec;
            nclust_mc(q)       = (double)sum(nj_val > 0);
         }
   }
   
   delete[] inds_eq_j;
   
   return Rcpp::List::create(
      Rcpp::Named("Lambda")    = Lambda_mc,
      Rcpp::Named("residuals") = ps_mat,
      Rcpp::Named("B")         = B_mc,
      Rcpp::Named("U")         = U_mc,
      Rcpp::Named("V")         = V_mc,
      Rcpp::Named("rank")      = (int)r,
      Rcpp::Named("prior_B")   = prior_B,
      Rcpp::Named("loglik")    = loglik_mc,
      Rcpp::Named("dir_prec")  = dir_prec_mc,
      Rcpp::Named("nclust")    = nclust_mc,
      Rcpp::Named("floor_hits")  = floor_hits,
      Rcpp::Named("floor_first") = floor_first,
      Rcpp::Named("floor_iters") = floor_iters
   );
}