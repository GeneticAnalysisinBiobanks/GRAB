#include <boost/math/distributions/chi_squared.hpp>
#include <RcppArmadillo.h>

inline void gethwepval(arma::vec GVec, double& hwepval, double hwepvalCutoff) {
  int n0 = 0, n1 = 0, n2 = 0, n = GVec.size();

  for (int i = 0; i < n; i++) {
    double g = GVec[i];
    if (g == 0) {
      n0++;
    } else if (g == 2) {
      n2++;
    } else if (g == 1) {
      n1++;
    } else {
      continue;
    }
  }
  
  double nsum = n0 + n1 + n2;
  double Gfreq = (0.5 * n1 + n2) / nsum;
  arma::vec obs = arma::vec(3);
  obs(0) = n0; obs(1) = n1; obs(2) = n2;

  arma::vec exp(3);
  exp(0) = (1 - Gfreq) * (1 - Gfreq) * nsum;
  exp(1) = 2 * Gfreq * (1 - Gfreq) * nsum;
  exp(2) = Gfreq * Gfreq * nsum;

  double hwechisq = arma::accu(arma::square(obs - exp) / exp);
  boost::math::chi_squared_distribution<double> chisq_dist(1.0);
  hwepval = boost::math::cdf(boost::math::complement(chisq_dist, hwechisq));
}
