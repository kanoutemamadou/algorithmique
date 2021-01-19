#include <Rcpp.h> //to use the NumericVector object
using namespace Rcpp; //to use the NumericVector object

#include<vector> //to use std::vector<double>
#include <iostream>
#include <string>

//' Similarity matrix algorithm using C++
//'
//' @param seq1 first sequence
//' @param seq1 second sequence
//' @return a similarity matrix
//' @export
// [[Rcpp::export]] //mandatory to export the function

NumericMatrix similarity_Rcpp(CharacterVector seq1, CharacterVector seq2, int match, int mismatch)
{
    int n = seq1.size();
    int m = seq2.size();
    NumericMatrix S(n, m);
    // Rprintf("the value of %i : %f \n", n);
    for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < m; j++)
      {
      	if (seq1[i] == seq2[j]) {
      	    S(i,j) = match;
      	} else {
      	   S(i,j) = mismatch;
      	}
      }
    }
    return(S);
}


//' Score matrix algorithm using C++
//'
//'@param d deletion
//' @param seq1 first sequence
//' @param seq1 second sequence
//' @return a similarity matrix
//' @export
// [[Rcpp::export]] //mandatory to export the function

NumericMatrix scoring_Rcpp(CharacterVector seq1, CharacterVector seq2, NumericMatrix S, int match, int mismatch, int d)
{
    int n = seq1.size();
    int m = seq2.size();  
    // NumericMatrix S(n, m);
    // S = similarity_Rcpp(seq1, seq2, match, mismatch);
    NumericMatrix Fij(n+1, m+1);
    
    for(int i = 0; i < (n+1); i++) {
      Fij(i,0) = i*d;
    } 
    
    for(int j = 0; j < (m+1); j++) {
      Fij(0,j) = (j)*d;
    }     
    int Match;
    int Delete;
    int Insert;
    
    for(int i = 1; i < (n+1); i++)
    {
      for(int j = 1; j < (m+1); j++)
      {
        Match = Fij((i-1),(j-1)) + S((i-1),(j-1));
        Delete = Fij((i-1),j)+d;
        Insert = Fij(i,(j-1))+d;
        Fij(i,j) = std::max({Match,Insert,Delete});
      }
    }
    return(Fij);
}


//'get Alignement algorithm using C++
//'
//'@param d deletion
//' @param seq1 first sequence
//' @param seq1 second sequence
//' @return a similarity matrix
//' @export
// [[Rcpp::export]] //mandatory to export the function

List NeedlemanWunsch_Rcpp(CharacterVector seq1, CharacterVector seq2, int match, int mismatch, int d)
{
    int n = seq1.size();
    int m = seq2.size();
    NumericMatrix S(n, m);
    NumericMatrix Fij(n+1, m+1);
    S = similarity_Rcpp(seq1, seq2, match, mismatch);
    Fij = scoring_Rcpp(seq1, seq2, S, match, mismatch, d);
    int i = seq1.size();
    int j = seq2.size();
    std::string AlignmentA("");
    std::string AlignmentB("");
    while(i > 0 || j > 0)
      {
        if((i>0)&&(j>0)&&(Fij(i,j)==Fij(i-1,j-1)+S(i-1,j-1)))
        {
          AlignmentA.insert(0,seq1[i-1]);
          AlignmentB.insert(0,seq2[j-1]);
          i = i-1;
          j = j-1;
        }
        else if((i>0)&&(Fij(i,j)==Fij(i-1,j)+d)){
          AlignmentA.insert(0,seq1[i-1]);
          AlignmentB.insert(0,"-");
          i = i-1;
        }
        else{
          AlignmentA.insert(0,"-");
            AlignmentB.insert(0,seq2[j-1]);
          j = j-1;
        }
    }
    return List::create(Named("Alignment A") = AlignmentA, Named("Alignment B") = AlignmentB, Named("Fij")=Fij, Named("score") = Fij(n, n));
}




int similarity2(CharacterVector seq1, CharacterVector seq2, int ind1, int ind2, int match, int mismatch) {
  int S = match;
  if (seq1[ind1] != seq2[ind2]) {
    S = mismatch;
  }
  return S;
}



//' Score matrix algorithm using C++
//'
//'@param d deletion
//' @param seq1 first sequence
//' @param seq1 second sequence
//' @return a similarity matrix
//' @export
// [[Rcpp::export]] //mandatory to export the function

NumericMatrix scoringV2_Rcpp(CharacterVector seq1, CharacterVector seq2, int match, int mismatch, int d)
{
  int n = seq1.size();
  int m = seq2.size();  
  
  NumericMatrix Fij(n+1, m+1);
  Fij.fill(-1000);
  int S = 1;
  
  Fij(0,0) = 0;
  Fij(0,1) = d;
  Fij(1,0) = d;
  
  int Match;
  int Delete;
  int Insert;
  
  for(int i = 1; i < (n+1); i++)
  {
    // S = similarity2(seq1[i-1], seq2[i-1], match, mismatch);
    S = similarity2(seq1, seq2, i-1, i-1, match, mismatch);
    Fij(i,i) = std::max({Fij(i-1,i) + d, Fij(i,i-1) + d, Fij(i-1,i-1) + S});
    
    if(i+1 < (n+1)) {
      // S = similarity2(seq1[i], seq2[i-1], match, mismatch);
      S = similarity2(seq1, seq2, i, i-1, match, mismatch);
      Fij(i+1,i) = std::max({Fij(i,i) + d, Fij(i,i-1) + S});
      // S = similarity2(seq1[i-1], seq2[i], match, mismatch);
      S = similarity2(seq1, seq2, i-1, i, match, mismatch);
      Fij(i,i+1) = std::max({Fij(i,i) + d, Fij(i-1,i) + S});
    }
  }
  return Fij;
}

// fillMatrix

//' Needleman V2 algorithm using C++
//'
//' @param seq1 first sequence
//' @param seq1 second sequence
//' @param match matching
//' @param mismatch mismatching
//' @param d deletion
//' @return a similarity matrix
//' @export
// [[Rcpp::export]] //mandatory to export the function

List NeedlemanWunschV2_Rcpp(CharacterVector seq1, CharacterVector seq2, int match, int mismatch, int d)
{
  
  int n = seq1.size();
  NumericMatrix Fij(n+1, n+1);
  Fij = scoringV2_Rcpp(seq1, seq2, match, mismatch, d);
  
  int i = seq1.size();
  int j = seq2.size();
  std::string AlignmentA("");
  std::string AlignmentB("");
  int S = 0;
  int t = i;
  while(i > 0 || j > 0)
  {
    if((i>0)&&(j>0)&&(Fij(i,j)==Fij(i-1,j-1)+similarity2(seq1, seq2, i-1, j-1, match, mismatch)))
    {
      AlignmentA.insert(0,seq1[i-1]);
      AlignmentB.insert(0,seq2[j-1]);
      i = i-1;
      j = j-1;
    }
    else if((i>0)&&(Fij(i,j)==Fij(i-1,j)+d)){
      AlignmentA.insert(0,seq1[i-1]);
      AlignmentB.insert(0,"-");
      i = i-1;
    } else{
      AlignmentA.insert(0,"-");
      AlignmentB.insert(0,seq2[j-1]);
      j = j-1;
    }
  }
  return List::create(Named("Alignment A") = AlignmentA, Named("Alignment B") = AlignmentB, Named("Fij")=Fij, Named("score") = Fij(n, n));
}
