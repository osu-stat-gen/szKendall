
#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;


// szKendall1 for loci pair (i,j) in single-cell 1 and loci pair (u,v) in single-cell 2 for all (i,j) and (u,v)
// [[Rcpp::export]]
double szkendall1(NumericVector Y1, NumericVector Y2, Nullable<IntegerVector> Y1_sz_idx, Nullable<IntegerVector> Y2_sz_idx,
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


    // This part is needed in szkendall1_cpp2 function, but they are given as input in the current function szkendall1_cpp.
    // NumericVector weight_vec(n);
    // NumericVector weight_sz_vec(n);
    //
    // for(int iter=0; iter<n; iter++){
    //   weight_vec[iter] = weight1(iter);
    //   weight_sz_vec[iter] = weight1_sz(iter);
    // }

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


