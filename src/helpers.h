
#ifndef HELPERS_H
#define HELPERS_H

#include <Rcpp.h>

// Declaration
Rcpp::IntegerVector sortedIndex(const Rcpp::NumericVector x);
Rcpp::IntegerVector compare_self(Rcpp::NumericVector x);
double tieXcount(Rcpp::NumericVector x);
double tieXYcount(Rcpp::NumericVector x, Rcpp::NumericVector y);
int kendall_discordant(Rcpp::IntegerVector x, Rcpp::IntegerVector y);
double kendall_distance_cpp(Rcpp::NumericVector x, Rcpp::NumericVector y);
double weight3(int n, int diff);
double weight3_sz(int diff);
Rcpp::NumericVector cal_weight_vec3(int n);
Rcpp::NumericVector cal_weight_sz_vec3(int n);
double weight1(int diff);
double weight1_sz(int diff);
Rcpp::NumericVector cal_weight_vec1(int n);
Rcpp::NumericVector cal_weight_sz_vec1(int n);
double weight2_sz(int diff);
Rcpp::NumericVector cal_weight_sz_vec2(int n);
int countPairsWithDiffK(Rcpp::IntegerVector arr, int n, int k);

#endif
