// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]]

// [[Rcpp::export]]
std::vector<double> squares(Rcpp::NumericVector data)
{
  // Rcpp::NumericVector squares(Rcpp::NumericVector data)
  // Rcpp::NumericVector result(data.size());
  unsigned int n = data.size();
  std::vector<double> result(n);
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(4); // Use 4 threads for all consecutive parallel regions
  
#pragma omp parallel
{
  Rcpp::Rcout << omp_get_num_threads() << std::endl;
  for (int i = 0; i < data.size(); i++) {
    result.at(i) = data[i] * data[i];
  }
  
}
return result;
}
