#include <stdexcept>
#include <limits>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/beta.hpp>
#include <RcppArmadillo.h>

#include "UTIL.h"

double hwe_exact(int obs_hets, int obs_hom1, int obs_hom2) {

  int rare = std::min(obs_hom1, obs_hom2);
  int rare_copies = 2*rare + obs_hets;
  int genotypes = obs_hets + obs_hom1 + obs_hom2;
  if (genotypes == 0) return std::numeric_limits<double>::quiet_NaN();

  std::vector<double> probs(rare_copies + 1, 0.0);

  int mid = rare_copies * (2*genotypes - rare_copies) / (2*genotypes);
  if ((rare_copies & 1) ^ (mid & 1)) mid++;

  int curr_hets = mid;
  int curr_homr = (rare_copies - curr_hets) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  probs[mid] = 1.0;
  double sum = probs[mid];

  // downward
  while (curr_hets > 1) {
      probs[curr_hets - 2] =
        probs[curr_hets] * curr_hets * (curr_hets - 1.0) /
        (4.0 * (curr_homr + 1) * (curr_homc + 1));

      sum += probs[curr_hets - 2];
      curr_hets -= 2;
      curr_homr++;
      curr_homc++;
  }

  // upward
  curr_hets = mid;
  curr_homr = (rare_copies - curr_hets) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;

  while (curr_hets <= rare_copies - 2) {
      probs[curr_hets + 2] =
        probs[curr_hets] * 4.0 * curr_homr * curr_homc /
        ((curr_hets + 2) * (curr_hets + 1.0));

      sum += probs[curr_hets + 2];
      curr_hets += 2;
      curr_homr--;
      curr_homc--;
  }

  // normalize
  for (double &p : probs) p /= sum;

  // p-value
  double pval = 0.0;
  for (int i = 0; i <= rare_copies; i += 2) {
      if (probs[i] <= probs[obs_hets]) {
        pval += probs[i];
      }
  }
  return std::min(1.0, pval);
}


void gethwepval(arma::vec GVec, double& hwepval, double hwepvalCutoff) {
  int n0 = 0, n1 = 0, n2 = 0, n = GVec.size();

  for (int i = 0; i < n; i++) {
    double g = GVec[i];
    if (g < hwepvalCutoff) {
      n0++;
    } else if (g > 2 - hwepvalCutoff) {
      n2++;
    } else if (g > 1 - hwepvalCutoff && g < 1 + hwepvalCutoff) {
      n1++;
    } else {
      continue;
    }
  }
  
  double nsum = n0 + n1 + n2;
  int rare = std::min(n0, n2);
  int rare_copies = 2 * rare + n1;

  if (nsum == 0) {
    hwepval = std::numeric_limits<double>::quiet_NaN();
  } else if (rare_copies < 100) {
    // exact
    hwepval = hwe_exact(n1, n0, n2);
  } else {
    // chi-square
    double Gfreq = (0.5 * n1 + n2) / nsum;

    arma::vec obs = arma::vec(3);
    obs(0) = n0; obs(1) = n1; obs(2) = n2;

    arma::vec exp(3);
    exp(0) = (1 - Gfreq) * (1 - Gfreq) * nsum;
    exp(1) = 2 * Gfreq * (1 - Gfreq) * nsum;
    exp(2) = Gfreq * Gfreq * nsum;

    if (arma::any(exp == 0)) {
      hwepval = hwe_exact(n1, n0, n2);
    } else {
      double hwechisq = arma::accu(arma::square(obs - exp) / exp);
      boost::math::chi_squared_distribution<double> chisq_dist(1.0);
      hwepval = boost::math::cdf(boost::math::complement(chisq_dist, hwechisq));
    }
  }
}
