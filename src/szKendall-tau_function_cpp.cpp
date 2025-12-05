
#include <Rcpp.h>
using namespace Rcpp;



//  [[Rcpp::export]]
IntegerVector sortedIndex(NumericVector x){
  IntegerVector idx = seq_along(x) - 1;  // idx is a sequence from 0 to length(x)-1. 
  
  std::stable_sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] < x[j];});
  
  return idx;  // return the indexes in the increasing order of x, which range from 0 to length(x)-1. 
  // for equal values, whichever appears first will have the smaller index. 
}


// [[Rcpp::export]]
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


//  [[Rcpp::export]]
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


//  [[Rcpp::export]]
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



// [[Rcpp::export]]
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



// [[Rcpp::export]]
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


// Weight function that depends on |j-i|
// [[Rcpp::export]]
double weight(int n, int diff){
  double weight_i_j = pow(n - 1.0 - abs(diff), -0.4);
  return(weight_i_j);
}

// Structural zero discrepancy score
// [[Rcpp::export]]
double weight_sz(int diff){
  double weight_sz_i_j = pow(1.0+abs(diff), -0.5);
  return(weight_sz_i_j);
}

// Calculate weight vector that depends on |j-i|
// [[Rcpp::export]]
NumericVector cal_weight_vec(int n){
  NumericVector weight_vec(n);
  
  for(int iter=0; iter<n; iter++){
    weight_vec[iter] = weight(n, iter);
  }
  
  return(weight_vec);
}


// Calculate SZ-weight vector that depends on |j-i|-|u-v|
// [[Rcpp::export]]
NumericVector cal_weight_sz_vec(int n){
  NumericVector weight_sz_vec(n);
  
  for(int iter=0; iter<n; iter++){
    weight_sz_vec[iter] = weight_sz(iter);
  }
  
  return(weight_sz_vec);
}


// [[Rcpp::export]]
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



// szkendall's tau for loci pair (i,j) in single-cell 1 and loci pair (u,v) in single-cell 2 for all (i,j) and (u,v)
// [[Rcpp::export]]
double szkendall_cpp(NumericVector Y1, NumericVector Y2, Nullable<IntegerVector> Y1_sz_idx, Nullable<IntegerVector> Y2_sz_idx, 
                      NumericVector weight_vec, NumericVector weight_sz_vec, String type="Nodiag"){
  if(Y1.length() != Y2.length()){
    stop("Error: The two single cell contact count vectors do not have the same length!"); 
  }else{
    
    int n = 0; 
    int len = Y1.length();  // at least 3 
    
    if(type == "Nodiag"){
      // In this case, len = n*(n-1)/2
      n = ceil(sqrt(2.0*len));  // number of bins, at least 2
    }else{
      // In this case, len = n*(n+1)/2
      n = floor(sqrt(2.0*len));  // number of bins, at least 2
    }
    
    //IntegerVector row_idx(len);
    //IntegerVector col_idx(len);
    IntegerVector rowminuscol_idx(len);
    
    if(type == "Nodiag"){
      //row_idx[0] = 1;
      //col_idx[0] = 2;
      rowminuscol_idx[0] = 1;
      for(int i=1; i<(n-1); i++){
        //row_idx[seq(sum(seq(1,i)), sum(seq(1,i+1))-1)] = seq(1,i+1);
        //col_idx[seq(sum(seq(1,i)), sum(seq(1,i+1))-1)] = rep(i+2,i+1);
        rowminuscol_idx[seq(sum(seq(1,i)), sum(seq(1,i+1))-1)] = rev(seq(1, i+1));
      }
    }else{
      //row_idx[0] = 1;
      //col_idx[0] = 1;
      rowminuscol_idx[0] = 0;
      for(int i=1; i<n; i++){
        //row_idx[seq(sum(seq(1,i)), sum(seq(1,i+1))-1)] = seq(1,i+1);
        //col_idx[seq(sum(seq(1,i)), sum(seq(1,i+1))-1)] = rep(i+1,i+1);
        rowminuscol_idx[seq(sum(seq(1,i)), sum(seq(1,i+1))-1)] = rev(seq(0, i));
      }
    }
    
    //row_idx = row_idx - 1;
    //col_idx = col_idx - 1;

    IntegerVector whole = seq(1,len);
    IntegerVector region1 = whole;
    IntegerVector region4 = whole;
    IntegerVector region2 = whole;
    IntegerVector region3 = whole;
    
    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int n4 = 0;
    
    if(Y1_sz_idx.isNull() & Y2_sz_idx.isNull()){
      
      region1 = whole;
      region1 = region1 - 1;
      region4 = NULL;
      region2 = NULL;
      region3 = NULL;
      
      n1 = region1.length();
      n2 = 0;
      n3 = 0;
      n4 = 0;
      
    }else if(Y1_sz_idx.isNotNull() & Y2_sz_idx.isNull()){
      
      IntegerVector Y1_sz_idx2(Y1_sz_idx);
      region1 = setdiff(whole, Y1_sz_idx2);
      region1 = region1.sort();
      region1 = region1 - 1;
      region4 = NULL;
      
      region2 = Y1_sz_idx2;
      region2= region2.sort();
      region2 = region2 - 1;
      region3 = NULL;
      
      n1 = region1.length();
      n2 = region2.length();
      n3 = 0;
      n4 = 0;
      
    }else if(Y1_sz_idx.isNull() & Y2_sz_idx.isNotNull()){
      
      IntegerVector Y2_sz_idx2(Y2_sz_idx); 
      region1 = setdiff(whole, Y2_sz_idx2);
      region1 = region1.sort();
      region1 = region1 - 1;
      region4 = NULL;
      
      region2 = NULL;
      region3 = Y2_sz_idx2;
      region3 = region3.sort();
      region3 = region3 - 1;
      
      n1 = region1.length();
      n2 = 0;
      n3 = region3.length();
      n4 = 0;
      
    }else{
      
      IntegerVector Y1_sz_idx2(Y1_sz_idx); 
      IntegerVector Y2_sz_idx2(Y2_sz_idx); 
      
      region1 = setdiff(whole, union_(Y1_sz_idx2, Y2_sz_idx2));
      region1 = region1.sort();
      region4 = (intersect(Y1_sz_idx2, Y2_sz_idx2));
      region4 = region4.sort();
      
      region2 = (setdiff(Y1_sz_idx2, region4));
      region2= region2.sort();
      region3 = (setdiff(Y2_sz_idx2, region4));
      region3 = region3.sort();
      
      // IntegerVector region23 = union_(region2, region3);
      // region23 = region23.sort();
      
      region1 = region1 - 1;
      region2 = region2 - 1;
      region3 = region3 - 1;
      // region23 = region23 - 1;
      region4 = region4 - 1;
      
      n1 = region1.length();
      n2 = region2.length();
      n3 = region3.length();
      n4 = region4.length();
    }
    
    
    double score = 0.0;
    
    // Rprintf("n1: %d\t n2: %d\t n3: %d\t n4: %d\n", n1, n2, n3, n4);
    
    double temp = 0.0; 
    
    // Case 1: Y_ij1 = +, Y_uv1 = +, Y_ij2 = +, Y_uv2 = + 
    if(n1 > 1){
      for(int r=0; r<n1-1; r++){
        for(int s=r+1; s<n1; s++){
          temp = (Y1[region1[r]]-Y1[region1[s]])*(Y2[region1[r]]-Y2[region1[s]]); 
          score = score + (1.0*(temp<0) + 0.5*(temp==0))*weight_vec[abs(rowminuscol_idx[region1[r]] - rowminuscol_idx[region1[s]])];
        }
      }
    }
    
    
    // Case 2: Y_ij1 = 0, Y_uv1 = +, Y_ij2 = +, Y_uv2 = + 
    if((n1 > 0) & (n2 > 0)){
      for(int r=0; r<n2; r++){
        for(int s=0; s<n1; s++){
          temp = (Y1[region2[r]]-Y1[region1[s]])*(Y2[region2[r]]-Y2[region1[s]]); 
          score = score + (1.0*(temp<0) + 0.5*(temp==0))*weight_vec[abs(rowminuscol_idx[region2[r]] - rowminuscol_idx[region1[s]])];
        }
        score = score + weight_sz_vec[abs(rowminuscol_idx[region2[r]])] * n1; 
      }
    }
    
    
    // Case 3: Y_ij1 = +, Y_uv1 = 0, Y_ij2 = +, Y_uv2 = +
    // same as case 2, don't double count! 
    
    // Case 4: Y_ij1 = +, Y_uv1 = +, Y_ij2 = 0, Y_uv2 = + 
    if((n1 > 0) & (n3 > 0)){
      for(int r=0; r<n3; r++){
        for(int s=0; s<n1; s++){
          temp = (Y1[region3[r]]-Y1[region1[s]])*(Y2[region3[r]]-Y2[region1[s]]); 
          score = score + (1.0*(temp<0) + 0.5*(temp==0))*weight_vec[abs(rowminuscol_idx[region3[r]] - rowminuscol_idx[region1[s]])];
        }
        score = score + weight_sz_vec[abs(rowminuscol_idx[region3[r]])] * n1;
      }
    }
    
    // Case 5: Y_ij1 = +, Y_uv1 = +, Y_ij2 = +, Y_uv2 = 0 
    // same as case 4. 
    
    // Case 6: Y_ij1 = 0, Y_uv1 = 0, Y_ij2 = +, Y_uv2 = + 
    if(n2 > 1){
      for(int r=0; r<n2; r++){
        score = score + weight_sz_vec[abs(rowminuscol_idx[region2[r]])] * (n2-1);
      }
    } 
    
    // Case 7: Y_ij1 = 0, Y_uv1 = +, Y_ij2 = +, Y_uv2 = 0 
    if((n2 > 0) & (n3 > 0)){
      for(int r=0; r<n2; r++){
        score = score + weight_sz_vec[abs(rowminuscol_idx[region2[r]])] * n3;
      }
      for(int s=0; s<n3; s++){
        score = score + weight_sz_vec[abs(rowminuscol_idx[region3[s]])] * n2;
      }
    }
    
    // Case 8: Y_ij1 = +, Y_uv1 = 0, Y_ij2 = 0, Y_uv2 = + 
    // same as case 7. 
    
    // Case 9: Y_ij1 = +, Y_uv1 = +, Y_ij2 = 0, Y_uv2 = 0 
    if(n3 > 1){
      for(int r=0; r<n3; r++){
        score = score + weight_sz_vec[abs(rowminuscol_idx[region3[r]])] * (n3-1);
      }
    }
    
    // Case 10: Y_ij1 = 0, Y_uv1 = +, Y_ij2 = 0, Y_uv2 = + 
    // in this case the distance is 0. 
    
    // Case 11: Y_ij1 = +, Y_uv1 = 0, Y_ij2 = +, Y_uv2 = 0 
    // same as case 10. 
    
    // Case 12: Y_ij1 = 0, Y_uv1 = 0, Y_ij2 = 0, Y_uv2 = + 
    if((n4 > 0) & (n2 > 0)){
      for(int s=0; s<n2; s++){
        score = score + weight_sz_vec[abs(rowminuscol_idx[region2[s]])] * n4;
      }
    }
    
    // Case 13: Y_ij1 = 0, Y_uv1 = 0, Y_ij2 = +, Y_uv2 = 0 
    // same as case 12.
    
    // Case 14: Y_ij1 = 0, Y_uv1 = +, Y_ij2 = 0, Y_uv2 = 0 
    if((n4 > 0) & (n3 > 0)){
      for(int s=0; s<n3; s++){
        score = score + weight_sz_vec[abs(rowminuscol_idx[region3[s]])] * n4;
      }
    }    
    
    // Case 15: Y_ij1 = +, Y_uv1 = 0, Y_ij2 = 0, Y_uv2 = 0 
    // same as case 14. 
    
    // Case 16: Y_ij1 = 0, Y_uv1 = +, Y_ij2 = 0, Y_uv2 = 0 
    // in this case the distance is 0. 
    
    return(score);
  }
}


