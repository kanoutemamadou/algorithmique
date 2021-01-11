#include <Rcpp.h> //to use the NumericVector object
using namespace Rcpp; //to use the NumericVector object

#include<vector> //to use std::vector<double>
#include <iostream>
#include <string>

//' Similarity matrix algorithm using C++
//'
//' @param seq1 first sequence
//' @param seq1 second sequence
//' @param n length of first sequence
//' @param m length of second sequence
//' @return a similarity matrix
//' @export
// [[Rcpp::export]] //mandatory to export the function

NumericMatrix similarity_Rcpp(CharacterVector seq1, CharacterVector seq2)
{
    int n = seq1.size();
    int m = seq2.size();
    NumericMatrix S(n, m);
    Rprintf("the value of %i : %f \n", n);
    for(int i = 0; i < n; i++)
    {
      for(int j = 0; j < m; j++)
      {
      	if (seq1[i] == seq2[j]) {
      	    S(i,j) = 1;
      	} else {
      	   S(i,j) = -1;
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
//' @param n length of first sequence
//' @param m length of second sequence
//' @return a similarity matrix
//' @export
// [[Rcpp::export]] //mandatory to export the function

NumericMatrix scoring_Rcpp(int d, CharacterVector seq1, CharacterVector seq2)
{
    int n = seq1.size();
    int m = seq2.size();  
    NumericMatrix S(n, m);
    S = similarity_Rcpp(seq1, seq2);
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
//' @param n length of first sequence
//' @param m length of second sequence
//' @return a similarity matrix
//' @export
// [[Rcpp::export]] //mandatory to export the function

List getAlignement_Rcpp(int d, CharacterVector seq1, CharacterVector seq2)
{

    int n = seq1.size();
    int m = seq2.size();
    NumericMatrix S(n, m);
    NumericMatrix Fij(n+1, m+1);
    S = similarity_Rcpp(seq1, seq2);
    Fij = scoring_Rcpp(d, seq1, seq2);
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
    return List::create(Named("Alignment A") = AlignmentA, Named("Alignment B") = AlignmentB);
}

