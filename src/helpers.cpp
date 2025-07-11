
#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


IntegerVector sortedIndex(NumericVector x){
  IntegerVector idx = seq_along(x) - 1;  // idx is a sequence from 0 to length(x)-1.

  std::stable_sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] < x[j];});

  return idx;  // return the indexes in the increasing order of x, which range from 0 to length(x)-1.
  // for equal values, whichever appears first will have the smaller index.
}


IntegerVector compare_self(NumericVector x){
  int n_entry = x.size();
  IntegerVector match_self (n_entry);
  match_self[0] = 1;

  int idx = 1;

  for (int i = 1; i < (n_entry); i++) {
    if (x[i] != x[(i - 1)]) {
      match_self[idx] = 1;  // if the i-th value of x is different from the previous one, its match_self index is 1.
    } else {
      match_self[idx] = 0;  // if the i-th value of x is the same the previous one, its match_self index is 0.
    }
    idx++;
  }
  return match_self;
}


double tieXcount(NumericVector x){

  NumericVector x2=x;
  x2 = x2[sortedIndex(x2)];
  double tieCount = 0.0;
  double m1 = 0.0;
  int k = 0;
  int n = x.size();

  for(k=1; k<n; k++){
    if(x2[k-1] == x2[k]){
      tieCount = tieCount+1;
    }else if(tieCount > 0){
      m1 = m1 + tieCount*(tieCount+1)/2;
      // Rprintf("%f\n", tieCount);
      tieCount = 0.0;
    }
  }

  if(tieCount > 0){
    // Rprintf("%f\n", tieCount);
    m1 = m1 + tieCount*(tieCount+1)/2;
  }

  return(m1);
}


double tieXYcount(NumericVector x, NumericVector y){

  NumericVector x2=x;
  NumericVector y2=y;

  // Obtain environment containing function order()
  Environment base("package:base");

  // Make function callable from C++
  Function order_r = base["order"];

  IntegerVector order_xy = order_r(x, y);

  // Rprintf("The orders: \n");
  // for(int i=0; i<order_xy.length(); ++i){
  //   Rprintf("%i \t", i, order_xy[i]);
  // }
  // Rprintf("\n");

  x2 = x2[order_xy-1];
  y2 = y2[order_xy-1];

  // Rprintf("Ordered x \t y \n");
  // for(int i=0; i<order_xy.length(); ++i){
  //   Rprintf("%f \t %f \n", x2[i], y2[i]);
  // }

  double tieCount = 0.0;
  double m2 = 0.0;
  int k = 0;
  int n = x.size();

  for(k=1; k<n; k++){
    if((x2[k-1] == x2[k]) & (y2[k-1] == y2[k])){
      tieCount = tieCount+1;
    }else if(tieCount > 0){
      // Rprintf("%f\n", tieCount);
      m2 = m2 + tieCount*(tieCount+1)/2;
      tieCount = 0.0;
    }
  }

  if(tieCount > 0){
    // Rprintf("%f\n", tieCount);
    m2 = m2 + tieCount*(tieCount+1)/2;
  }

  return(m2);
}


int kendall_discordant(IntegerVector x, IntegerVector y){
  // count the number of discordant pairs
  NumericVector sup2 = {0.0, 0.0};
  sup2[0] = 1 + max(x);
  sup2[1] = 1 + max(y);
  double sup = max(sup2);

  IntegerVector arr(sup, 0);  // same as rep(0, sup).
  double i = 0;
  double k = 0;
  int n = x.size();
  int idx = 0;
  int dis = 0;

  while (i < n){
    while ((k < n) & (x[i] == x[k])) {
      dis = dis + i;
      idx = y[k];
      while (idx != 0) {
        dis = dis - arr[idx];
        idx = idx & (idx - 1);
      }
      k++;
    }
    while (i < k) {
      idx = y[i];
      while (idx < sup) {
        arr[idx] = arr[idx] + 1;
        idx = idx + (idx & (-1*idx));
      }
      i++;
    }
  }
  return dis;
}


double kendall_distance_cpp(NumericVector x, NumericVector y){

  NumericVector x2 = clone(x);
  NumericVector y2 = clone(y);

  IntegerVector perm_y = sortedIndex(y2);
  x2 = x2[perm_y];
  y2 = y2[perm_y];
  IntegerVector y3 = compare_self(y2);
  IntegerVector y4 = cumsum(y3);
  //return y4;

  IntegerVector perm_x = sortedIndex(x2);
  x2 = x2[perm_x];
  y4 = y4[perm_x];
  IntegerVector x3 = compare_self(x2);
  IntegerVector x4 = cumsum(x3);

  double dis = kendall_discordant(x4, y4);

  dis = dis + 0.5*(tieXcount(x)+tieXcount(y)-tieXYcount(x, y));

  return dis;
}



//-----------------------------------------------------------------------------------------

// For szKendall:

// Kendall's tau Weight function that depends on |j-i|-|u-v|
double weight3(int n, int diff){
  double weight_i_j = pow(n - 1.0 - abs(diff), -0.4);
  return(weight_i_j);
}

// Structural zero discrepancy score
double weight3_sz(int diff){
  double weight_sz_i_j = pow(1.0+abs(diff), -0.5);
  return(weight_sz_i_j);
}

// Calculate Kendall's tau weight vector that depends on |j-i|
// [[Rcpp::export]]
NumericVector cal_weight_vec3(int n){
  NumericVector weight_vec(n);

  for(int iter=0; iter<n; iter++){
    weight_vec[iter] = weight3(n, iter);
  }

  return(weight_vec);
}

// Calculate SZ-weight vector that depends on |j-i|
// [[Rcpp::export]]
NumericVector cal_weight_sz_vec3(int n){
  NumericVector weight_sz_vec(n);

  for(int iter=0; iter<n; iter++){
    weight_sz_vec[iter] = weight3_sz(iter);
  }

  return(weight_sz_vec);
}

//-----------------------------------------------------------------------------------------

// For szKendall1:

// Kendall's tau Weight function that depends on |j-i|-|u-v|
double weight1(int diff){
  double weight_i_j = pow(1.0+abs(diff), -0.4);
  return(weight_i_j);
}

// Structural zero discrepancy score
double weight1_sz(int diff){
  double weight1_sz_i_j = pow(1.0+abs(diff), -0.5);
  return(weight1_sz_i_j);
}

// Calculate weight vector that depends on |j-i|
// [[Rcpp::export]]
NumericVector cal_weight_vec1(int n){
  NumericVector weight_vec(n);

  for(int iter=0; iter<n; iter++){
    weight_vec[iter] = weight1(iter);
  }

  return(weight_vec);
}

// Calculate SZ-weight vector that depends on |j-i|-|u-v|
// [[Rcpp::export]]
NumericVector cal_weight_sz_vec1(int n){
  NumericVector weight_sz_vec(n);

  for(int iter=0; iter<n; iter++){
    weight_sz_vec[iter] = weight1_sz(iter);
  }

  return(weight_sz_vec);
}

//-----------------------------------------------------------------------------------------

// For szKendall2:

// Structural zero discrepancy score
double weight2_sz(int diff){
  double weight_sz_i_j = pow(1.0+abs(diff), -0.5);
  return(weight_sz_i_j);
}

// Calculate SZ-weight vector that depends on |j-i|
// [[Rcpp::export]]
NumericVector cal_weight_sz_vec2(int n){
  NumericVector weight_sz_vec(n);

  for(int iter=0; iter<n; iter++){
    weight_sz_vec[iter] = weight2_sz(iter);
  }

  return(weight_sz_vec);
}

//-----------------------------------------------------------------------------------------



int countPairsWithDiffK(IntegerVector arr, int n, int k)
{
  int count = 0;  // Initialize count
  int MAX = n;

  // Initialize empty hashmap.
  LogicalVector hashmap (MAX, false);

  // Insert array elements to hashmap
  for (int i = 0; i < n; i++) {
    hashmap[arr[i]] = true;
  }

  for (int i = 0; i < n; i++) {
    int a = arr[i];
    if (a - k >= 0 && hashmap[a - k])
      count++;
    if (a + k < MAX && hashmap[a + k])
      count++;
    hashmap[a] = false;
  }
  return count;
}


