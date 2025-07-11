
#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


// szKendall2 for loci pair (i,j) in single-cell 1 and loci pair (u,v) in single-cell 2 for all (i,j) and (u,v)
// [[Rcpp::export]]
double szkendall2(NumericVector Y1, NumericVector Y2, Nullable<IntegerVector> Y1_sz_idx, Nullable<IntegerVector> Y2_sz_idx,
                  NumericVector weight_sz_vec, String type="Nodiag"){
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


    // This part is needed in szkendall3_cpp2 function, but they are given as input in the current function szkendall2_cpp.
    // NumericVector weight_sz_vec(n);
    //
    // for(int iter=0; iter<n; iter++){
    //   weight_sz_vec[iter] = weight2_sz(iter);
    // }

    IntegerVector whole = seq(1,len);
    IntegerVector region1 = whole;
    IntegerVector region4 = whole;
    IntegerVector region2 = whole;
    IntegerVector region3 = whole;
    IntegerVector region23 = whole;

    int n1 = 0;
    int n2 = 0;
    int n3 = 0;
    int n4 = 0;
    int n23 = 0;

    if(Y1_sz_idx.isNull() & Y2_sz_idx.isNull()){

      region1 = whole;
      region1 = region1 - 1;
      region4 = NULL;
      region2 = NULL;
      region3 = NULL;
      region23 = NULL;

      n1 = region1.length();
      n2 = 0;
      n3 = 0;
      n4 = 0;
      n23 = 0;

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
      region23 = region2;

      n1 = region1.length();
      n2 = region2.length();
      n3 = 0;
      n4 = 0;
      n23 = n2;

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
      region23 = region3;

      n1 = region1.length();
      n2 = 0;
      n3 = region3.length();
      n4 = 0;
      n23 = n3;

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

      region23 = union_(region2, region3);
      region23 = region23.sort();

      region1 = region1 - 1;
      region2 = region2 - 1;
      region3 = region3 - 1;
      region4 = region4 - 1;
      region23 = region23 - 1;

      n1 = region1.length();
      n2 = region2.length();
      n3 = region3.length();
      n4 = region4.length();
      n23 = region23.length();
    }

    double score = 0.0;

    double tmp = 0.0;
    if(n23 > 0){
      for(int r=0; r<n23; r++){
        tmp = tmp+weight_sz_vec[abs(rowminuscol_idx[region23[r]])];
      }
    }
    score = score + tmp*(n1+n4);


    tmp = 0.0;
    if(n23 > 1){
      for(int r=0; r<n23; r++){
        tmp = tmp + weight_sz_vec[abs(rowminuscol_idx[region23[r]])] * (n23-1);
      }
    }
    score = score + tmp;


    for(int r=0; r<(len-1); r++){
      for(int s=r+1; s<len; s++){
        tmp = (Y1[r]-Y1[s])*(Y2[r]-Y2[s]);
        score = score + (1.0*(tmp<0) + 0.5*(tmp==0));
      }
    }

    return(score);
  }
}


