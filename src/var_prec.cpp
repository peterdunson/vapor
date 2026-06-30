#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <cmath>
#include <cfloat>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// ── Tolerances and macros ────────────────────────────────────────────────────
#define tolr        1e-5
#define ZTOL        (DBL_EPSILON * 10.0)
#define crossprod(x)  symmatu((x).t() * (x))
#define tcrossprod(x) symmatu((x) * (x).t())

// ── GSL replacements (drop external GSL dependency) ──────────────────────────
#ifndef GSL_MAX_DBL
#define GSL_MAX_DBL(a, b) (std::max((double)(a), (double)(b)))
#endif
#ifndef GSL_MIN_INT
#define GSL_MIN_INT(a, b) (std::min((int)(a), (int)(b)))
#endif

static inline double gsl_pow_2(double x) { return x * x; }

// log Pochhammer: log Gamma(a + x) - log Gamma(a)
static inline double gsl_sf_lnpoch(double a, double x) {
   return std::lgamma(a + x) - std::lgamma(a);
}

// fuzzy comparison: returns 0 if |a-b| within eps*|.|, else +/-1 (sign of a-b)
static inline int gsl_fcmp(double a, double b, double eps) {
   double diff = a - b;
   double m = std::max(std::fabs(a), std::fabs(b));
   if (std::fabs(diff) <= eps * m) return 0;
   return (diff > 0) ? 1 : -1;
}

// ============================================================================
//  GIG SAMPLER (lifted from gg_model_sub.cpp)
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
   int i, count = 0;
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
         ++count;
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
   int i, count = 0;
   t  = 0.5  * (lambda - 1.);
   s  = 0.25 * omega;
   xm = _gig_mode(lambda, omega);
   nc = t * log(xm) - s * (xm + 1./xm);
   ym = ((lambda+1.) + sqrt((lambda+1.)*(lambda+1.) + omega*omega)) / omega;
   um = exp(0.5*(lambda+1.)*log(ym) - s*(ym+1./ym) - nc);
   for (i = 0; i < n; i++) {
      do {
         ++count;
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
   int i, count = 0;
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
         ++count;
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

// ── Inverse-Gaussian sampler ─────────────────────────────────────────────────
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

// ── Log-marginal PDF for DP residual variance update ─────────────────────────
inline double log_marginal_pdf(unsigned n, double as, double bs, double vnorm) {
   return gsl_sf_lnpoch(as, ((double)n)/2.)
   + as*log(bs)
   - (as + ((double)n)/2.) * log(bs + vnorm/2.);
}

// ── Log-sum-exp ──────────────────────────────────────────────────────────────
double log_sum_exp(const arma::vec &x) {
   double mx = x.max();
   return log(sum(exp(x - mx))) + mx;
}

// ============================================================================
//  BLOCK 2 HELPERS  (Omega update — from gg_model_sub.cpp)
// ============================================================================

// Update loading params: psi_inv, T_phi, Plam
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

// Update loadings Lambda (the precision factor loadings, NOT the VAR B)
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
//  BLOCK 1 HELPER  (B update — adapted from bayesianVARs sample_PHI_factor)
// ============================================================================

// Fast equation-by-equation sampler for B under DL prior
// Y_res = Y - X*B (current residuals for Omega block)
// We sample B | Y, Omega using GLS transformation
void sample_B(arma::mat &B, const arma::mat &Y, const arma::mat &X,
              const arma::mat &Plam_B, const arma::vec &ps,
              const bool huge = false) {
   const int p = Y.n_cols;
   const int K = X.n_cols;
   const int T = Y.n_rows;
   
   for (int j = 0; j < p; ++j) {
      // scale by residual precision
      double sigma_j = 1.0 / sqrt(ps(j));
      arma::vec y_norm = Y.col(j) * sigma_j;
      arma::mat X_norm = X * sigma_j;
      
      arma::vec prior_var = Plam_B.row(j).t();  // K-dim prior variance for B[j,]
      
      if (huge) {
         // Bhattacharya et al. (2016) fast sampler: O(KT) instead of O(K^3)
         arma::vec u(K); u.imbue(R::norm_rand);
         u %= sqrt(prior_var);
         arma::vec delta(T); delta.imbue(R::norm_rand);
         arma::vec v_vec = X_norm * u + delta;
         arma::mat XPV = X_norm.each_row() % prior_var.t();
         arma::vec w   = arma::solve(XPV*X_norm.t() + arma::eye(T,T), y_norm - v_vec);
         B.row(j)      = (u + XPV.t()*w).t();
      } else {
         // standard posterior
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
         // sample lambda (local scale)
         double lam = rgig(a_B - 0.5,
                           coefs_sq(h) / GSL_MAX_DBL(lambda_j(h), 1e-10),
                           2.0 * xi_j);
         lambda_j(h) = GSL_MAX_DBL(lam, 9e-6);
         
         // sample psi (auxiliary)
         double ps = rgig(0.5,
                          coefs_sq(h) / GSL_MAX_DBL(lambda_j(h)*lambda_j(h), 1e-20),
                          1.0);
         psi_j(h) = GSL_MAX_DBL(ps, 9e-6);
      }
      
      // xi fixed (DL_plus=FALSE), no update needed
      T_phi_B.row(j) = lambda_j;
      
      // prior variance = psi * lambda^2
      Plam_B.row(j) = psi_j % square(lambda_j);
   }
}
// ============================================================================
//  DP RESIDUAL VARIANCE UPDATE  (from gg_model_sub.cpp prec_est)
// ============================================================================

double update_dir_prec(double alpha, double a, double b, unsigned n, unsigned k) {
   double phi  = R::rbeta(alpha+1., (double)n);
   double gam1 = a + k, gam2 = b - log(phi);
   double pi   = (gam1-1.) / (n*gam2);
   return (log(randu()) < log(pi)-log1p(pi))
      ? randg(distr_param(gam1,   1./gam2))
         : randg(distr_param(gam1-1.,1./gam2));
}

// ============================================================================
//  MAIN SAMPLER
// ============================================================================

// [[Rcpp::export]]
Rcpp::List var_precision_sampler(
      const arma::mat &Y,       // T x p data matrix
      const int lags    = 1,    // VAR lag order
      const int k       = 5,    // number of precision factors
      const int nrun    = 5000, // total MCMC iterations
      const int burn    = 2500, // burn-in (stored but flagged)
      const int thin    = 5,    // thinning interval
      const double a_L  = 0.5,  // DL shape for Lambda (precision loadings)
      const double b_L  = 0.5,  // DL rate  for Lambda
      const double a_B  = 0.5,  // DL shape for B (VAR coefficients)
      const double b_B  = 0.5,  // DL rate  for B
      const double as   = 0.1,  // DP base measure shape
      const double bs   = 0.1,  // DP base measure rate
      const double a_dir= 0.1,  // DP concentration prior shape
      const double b_dir= 0.1   // DP concentration prior rate
) {
   // ── Dimensions ───────────────────────────────────────────────────────────
   const unsigned T = Y.n_rows;
   const unsigned p = Y.n_cols;
   const unsigned Teff = T - lags;          // effective sample size
   const unsigned Kx   = lags * p;          // number of predictors per equation
   const double a_B_eff = 1.0 / (double)(p * Kx);  // DL shape for B: 1/K following Kastner & Huber
   
   // ── Build lagged data matrices ────────────────────────────────────────────
   // X: Teff x Kx   (lagged predictors)
   // Ymat: Teff x p  (responses)
   arma::mat Ymat(Teff, p), X(Teff, Kx);
   for (unsigned t = 0; t < Teff; ++t) {
      Ymat.row(t) = Y.row(t + lags);
      for (int lag = 0; lag < lags; ++lag)
         X.row(t).subvec(lag*p, (lag+1)*p - 1) = Y.row(t + lags - 1 - lag);
   }
   
   // ── Initialise B ─────────────────────────────────────────────────────────
   arma::mat B(p, Kx, fill::zeros);
   arma::mat T_phi_B(p, Kx, fill::value(1.0/Kx));   // lambda init to 1/K
   arma::mat Plam_B(p, Kx, fill::value(1.0/Kx));    // V_i = psi*lambda^2 ~ 1/K initially
   arma::vec tau_B(p, fill::value(0.5));             // xi init to 0.5
   
   // ── Initialise Lambda (precision factor loadings) ─────────────────────────
   arma::mat phimat(p, k, fill::randn);     // p x k loadings for Omega
   arma::mat T_phi  = randg<arma::mat>(p, k, distr_param(a_L, 0.5));
   T_phi.transform([](double v){ return GSL_MAX_DBL(v, 9e-17); });
   arma::mat Plam   = 1. / (1. + square(T_phi));
   
   // ── Initialise DP residual variances ─────────────────────────────────────
   const unsigned nmix = p;
   double dir_prec = 1.0;
   double prob_empty = log(dir_prec);
   arma::uvec nj_val(nmix, fill::zeros), del;
   arma::vec  log_nj_val(nmix), probs(nmix), log_probs(nmix);
   arma::Col<int> d;
   arma::uvec *inds_eq_j = new arma::uvec[nmix];
   arma::vec ps(p, fill::ones);
   
   // ── Storage ───────────────────────────────────────────────────────────────
   const unsigned n_mc = (unsigned)std::floor((double)nrun / thin);
   arma::cube  Lambda_mc(p, k, n_mc, fill::zeros);
   arma::mat   ps_mat(n_mc, p, fill::zeros);
   arma::cube  B_mc(p, Kx, n_mc, fill::zeros);
   
   // ── Progress bar ──────────────────────────────────────────────────────────
   Progress prog(nrun, true);
   
   // ── Main Gibbs loop ───────────────────────────────────────────────────────
   for (int i = 0; i < nrun; ++i) {
      
      if (Progress::check_abort()) { delete[] inds_eq_j; return R_NilValue; }
      prog.increment();
      
      // ── BLOCK 1: Sample B | Y, Omega ─────────────────────────────────────
      // Use current ps (residual precisions) as the equation-level variances
      bool huge = (Kx > (int)Teff);
      sample_B(B, Ymat, X, Plam_B, ps, huge);
      
      // Update DL hyperparameters for B
      update_B_params(B, T_phi_B, Plam_B, tau_B, a_B_eff, b_B);
      
      // ── Compute residuals given current B ─────────────────────────────────
      arma::mat E = Ymat - X * B.t();   // Teff x p
      
      // ── BLOCK 2: Sample Omega | residuals E ───────────────────────────────
      // Step 2a: Sample u and v (PrecFact latent variable representation)
      arma::mat Lmsg  = phimat.each_col() / ps;
      arma::mat Lmsg2 = symmatu(Lmsg.t() * phimat);
      Lmsg2.diag()   += 1.0;
      
      arma::mat noise(Teff, k, fill::randn);
      arma::mat chol_prec_u = chol(Lmsg2);
      arma::mat u = noise * trimatu(chol_prec_u);
      
      arma::mat tmpmt = solve(trimatl(chol_prec_u.t()), Lmsg.t());
      arma::mat v     = E + noise * tmpmt;
      arma::vec ssqv  = sum(square(v), 0).t();
      
      // Step 2b: Update loading params (T_phi, Plam)
      update_loading_params(phimat, T_phi, Plam, a_L, b_L);
      
      // Step 2c: Update loadings phimat
      update_phi(phimat, Plam, u, v, ssqv);
      
      // ── BLOCK 3: Sample residual variances (DP) ───────────────────────────
      if (!i) {
         // initialise clustering on first iteration
         unsigned n_init = GSL_MIN_INT((unsigned)std::ceil(log((double)p)), 2u);
         arma::mat dat(1, p); dat.row(0) = ssqv.t();
         del.set_size(p);
         // simple random init
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
            // current cluster
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
         
         // regenerate unique ps values
         for (unsigned j = 0; j < nmix; ++j) {
            if (nj_val(j)) {
               inds_eq_j[j] = find(del == j);
               double atom = randg(distr_param(as + 0.5*(Teff*nj_val(j)),
                                               1./(bs + 0.5*v_atoms(j))));
               ps(inds_eq_j[j]).fill(atom);
            } else inds_eq_j[j].reset();
         }
         dir_prec   = update_dir_prec(dir_prec, a_dir, b_dir, p,
                                      (unsigned)sum(nj_val > 0));
         prob_empty = log(dir_prec);
         
         // ── Store every thin-th sample ─────────────────────────────────────
         int rem = (i+1) % thin;
         if (rem == 0) {
            int q = (int)std::floor((double)(i+1)/thin) - 1;
            Lambda_mc.slice(q) = phimat;
            ps_mat.row(q)      = ps.t();
            B_mc.slice(q)      = B;
         }
   }
   
   delete[] inds_eq_j;
   
   return Rcpp::List::create(
      Rcpp::Named("Lambda")    = Lambda_mc,   // p x k x n_mc
      Rcpp::Named("residuals") = ps_mat,      // n_mc x p  (for fact_prec_processing)
      Rcpp::Named("B")         = B_mc         // p x Kx x n_mc
   );
}
