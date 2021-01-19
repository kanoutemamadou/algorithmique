
library(microbenchmark)
library(ggplot2)

##  MIT License
## Copyright (c) 2020 Projet6 Algorithmique

# library(microbenchmark)
# library(ggplot2)


#' Similarity  matrix R
#'
#' @description Creating similarity matrix
#' @param A first sequence
#' @param B second sequence
#' @param match matching
#' @param mismatch mismatching
#' @return a similarity matrix

similarity<-function(A,B, match, mismatch)
{
  n<-length(A)
  m<-length(B)
  S<-matrix(0,nrow = n,ncol = m,dimnames = list(A,B))
  for(i in 1:n)
  {
    for(j in 1:m)
    {
      S[i,j]<-ifelse(A[i]==B[j],match,mismatch)
    }
  }
  return(S)
}


#' Score  matrix R
#'
#' @description Creating score matrix
#' @param A first sequence
#' @param B second sequence
#' @param match matching
#' @param mismatch mismatching
#' @param d deletion
#' @return a score matrix

Scoring<-function(A,B, match, mismatch, d=-1)
{
  n<-length(A)
  m<-length(B)
  S<-similarity(A=A,B=B, match, mismatch)
  Fij<-matrix(0,nrow = n+1,ncol = m+1,dimnames = list(c(" ",A),c(" ",B)))
  Fij[1,2:(m+1)]<-seq(1,m)*d
  Fij[2:(n+1),1]<-seq(1,n)*d
  
  for(i in 2:(n+1))
  {
    for(j in 2:(m+1))
    {
      Match<-Fij[(i-1),(j-1)]+S[(i-1),(j-1)]
      Delete<-Fij[(i-1),j]+d
      Insert<-Fij[i,(j-1)]+d
      Fij[i,j]<-max(Match,Insert,Delete)
    }
  }
  return(Fij)
}


#' Alignments R
#'
#' @description Get best alignments of two sequences
#' @param A first sequence
#' @param B second sequence
#' @param match matching
#' @param mismatch mismatching
#' @param d deletion
#' @return alignments of two sequence by score max

NeedlemanWunsch<-function(A,B, match, mismatch, d=-1)
{

  
  if (typeof(A) == "character" && length(A) == 1) {
    A = unlist(strsplit(A, split=""))
  }
  
  if (typeof(B) == "character" && length(B) == 1) {
    B = unlist(strsplit(B, split=""))
  }  
  
  S<-similarity(A=A,B=B, match, mismatch)
  Fij<-Scoring(A=A,B=B, match, mismatch, d)
  n<-length(A)
  m<-length(B)
  AlignmentA<-" "
  AlignmentB<-" "
  i<-n
  j<-m
  
  while(i>0||j>0)
  {
    
    if((i>0)&&(j>0)&&(Fij[i+1,j+1]==Fij[i,j]+S[i,j]))
    {
      AlignmentA<-paste(A[i],AlignmentA)
      AlignmentB<-paste(B[j],AlignmentB)
      i<-i-1
      j<-j-1
    } else if((i>0)&&(Fij[i+1,j+1]==Fij[i,j+1]+d)){
      AlignmentA<-paste(A[i],AlignmentA)
      AlignmentB<-paste("-",AlignmentB)
      i<-i-1
    } else{
      AlignmentA<-paste("-",AlignmentA)
      AlignmentB<-paste(B[j],AlignmentB)
      j<-j-1
    }
  }
  return(list(AlignmentA=AlignmentA,AlignmentB=AlignmentB, Fij=Fij, score = Fij[n+1, m+1]))
  
}



similarity2 = function(S1, S2, match, mismatch) {
  S = match
  if (S1 != S2) {
    S = mismatch
  }
  return (S)
}


#' Score  matrix V2 R
#'
#' @description Creating score matrix
#' @param A first sequence
#' @param B second sequence
#' @param match matching
#' @param mismatch mismatching
#' @param d deletion
#' @return a score matrix V2


ScoringV2<-function(A,B, match, mismatch, d)
{
  
  # A = c(" ", A)
  # B = c(" ", B)
  n = length(A)
  colNames = B
  rowNames = A
  matF<-matrix(-Inf,nrow = n+1,ncol = n+1,dimnames = list(c(" ",A),c(" ",B)))
  # matF = matrix(-Inf, n, n)
  # row.names(matF) = rowNames
  # colnames(matF) = colNames
  
  
  matF[1,1] = 0
  matF[1,2] = d
  matF[2,1] = d
  
  for(i in 2:(n+1)) {
    
    S = similarity2(A[i-1], B[i-1], match, mismatch)
    matF[i,i] = max(matF[i-1,i] + d, matF[i,i-1] + d, matF[i-1,i-1] + S)
    
    if(i+1 <= (n+1)) {
      S = similarity2(A[i], B[i-1], match, mismatch)
      matF[i+1,i] = max(matF[i,i] + d, matF[i,i-1] + S)
      S = similarity2(A[i-1], B[i], match, mismatch)
      matF[i,i+1] = max(matF[i,i] + d, matF[i-1,i] + S)
    }
  }
  return (matF)  
}


#' Alignments V2 R
#'
#' @description Get best alignments of two sequences V2
#' @param A first sequence
#' @param B second sequence
#' @param match matching
#' @param mismatch mismatching
#' @param d deletion
#' @return alignments of two sequence by score max V2


NeedlemanWunschV2<-function(A,B, match, mismatch, d)
{

  if (typeof(A) == "character" && length(A) == 1) {
    A = unlist(strsplit(A, split=""))
  }
  
  if (typeof(B) == "character" && length(B) == 1) {
    B = unlist(strsplit(B, split=""))
  }  
 
  # match = 1
  # mismatch = -1  
   
  Fij = ScoringV2(A, B, match, mismatch, d) 
  n<-length(A)
  m<-length(B)
  AlignmentA<-" "
  AlignmentB<-" "
  i<-n
  j<-m
  while(i>0||j>0)
  {
    if((i>0)&&(j>0)&&(Fij[i+1,j+1]==Fij[i,j]+similarity2(A[i], B[j], match, mismatch)))
    {
      AlignmentA<-paste(A[i],AlignmentA)
      AlignmentB<-paste(B[j],AlignmentB)
      i<-i-1
      j<-j-1
    } else if((i>0)&& (Fij[i+1,j+1]==Fij[i,j+1]+d)){
      AlignmentA<-paste(A[i],AlignmentA)
      AlignmentB<-paste("-",AlignmentB)
      i<-i-1
    } else{
      AlignmentA<-paste("-",AlignmentA)
      AlignmentB<-paste(B[j],AlignmentB)
      j<-j-1
    }
  }
  return(list(AlignmentA=AlignmentA,AlignmentB=AlignmentB, Fij=Fij, score = Fij[n+1, n+1]))
  
}


SimulateSeq = function(n,m) {
  s <- sample(c("A","C","G","T"),size = n, replace = TRUE)
  snew <- s
  for(i in 0:m)
    
  {
    snew <- append(snew, sample(c("A","C","G","T"), 1), sample(length(snew), 1))
    snew <- snew[-sample(length(snew), 1)]
    snew[sample(length(snew),1)] <- sample(c("A","C","G","T"), 1)
    
  }
  s <- paste(s, collapse = "")
  snew <- paste(snew, collapse = "")  
  return (list(A=s, B=snew))  
}

timeByFunction = function(n = 1000, match, mismatch, d, fun="V1_R") {
  v = SimulateSeq(n,1)
  # print("v")
  # print(v)
  if (fun == "V1_R") {
    t = system.time(NeedlemanWunsch(v$A, v$B, match, mismatch, d))[[1]]
  } else if(fun == "V2_R") {
    t = system.time(NeedlemanWunschV2(v$A, v$B, match, mismatch, d))[[1]]
  } else if(fun == "V1_Cpp") {
    t = system.time(NeedlemanWunsch_Rcpp(v$A, v$B, match, mismatch, d))[[1]]
  }  else if(fun == "V2_Cpp") {
    t = system.time(NeedlemanWunschV2_Rcpp(v$A, v$B, match, mismatch, d))[[1]]
  } 
  return (t)
}

# n <- 1000
# 
# res <- microbenchmark(timeByFunction(n, match, mismatch, d, fun="V1_Cpp"), timeByFunction(n, match, mismatch, d, fun="V2_Cpp"), times = 50)
# autoplot(res)

nbSimus <- 20
vector_n <- seq(from = 10000, to = 100000, length.out = nbSimus)
nbRep <- 50
res_Heap <- data.frame(matrix(0, nbSimus, nbRep + 1))
colnames(res_Heap) <- c("n", paste0("Rep",1:nbRep))

j <- 1
for(i in vector_n)
{
  res_Heap[j,] <- c(i, replicate(nbRep, timeByFunction(n, match, mismatch, d, fun="V2_Cpp")))  
  #print(j)
  j <- j + 1
}

res <- rowMeans(res_Heap[,-1])
plot(vector_n, res, type = 'b', xlab = "data length", ylab = "mean time in seconds")


