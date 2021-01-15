similarity<-function(A,B)
{
  n<-length(A)
  m<-length(B)
  S<-matrix(0,nrow = n,ncol = m,dimnames = list(A,B))
  for(i in 1:n)
  {
    for(j in 1:m)
    {
      S[i,j]<-ifelse(A[i]==B[j],1,-1)
    }
  }
  return(S)
}

Scoring<-function(d,A,B)
{
  n<-length(A)
  m<-length(B)
  S<-similarity(A=A,B=B)
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
  return(t(Fij))
}


Needleman.Wunsch<-function(d=-1,A,B)
{

  
  if (typeof(A) == "character" && length(A) == 1) {
    A = unlist(strsplit(A, split=""))
  }
  
  if (typeof(B) == "character" && length(B) == 1) {
    B = unlist(strsplit(B, split=""))
  }  
  
  S<-similarity(A=A,B=B)
  Fij<-Scoring(d=d,A=A,B=B)
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
      AlignmentA<-trimws(paste0(A[i],AlignmentA))
      AlignmentB<-trimws(paste0(B[j],AlignmentB))
      i<-i-1
      j<-j-1
    }else if((i>0)&&(Fij[i+1,j+1]==Fij[i,j+1]+d)){
      AlignmentA<-trimws(paste0(A[i],AlignmentA))
      AlignmentB<-trimws(paste0("-",AlignmentB))
      i<-i-1
    }else{
      AlignmentA<-trimws(paste0("-",AlignmentA))
      AlignmentB<-trimws(paste0(B[j],AlignmentB))
      j<-j-1
    }
  }
  # cat("\n")
  # print("F dans V1")
  # print(Fij)
  # cat("\n")
  return(list(AlignmentA=AlignmentA,AlignmentB=AlignmentB, Fij=Fij, score = Fij[n+1, m+1]))
  
}



similarity2 = function(S1, S2, match, mismatch) {
  S = match
  if (S1 != S2) {
    S = mismatch
  }
  return (S)
}


ScoringV2<-function(gap,A,B, match, mismatch)
{
  
  match = 1
  mismatch = -1
  
  A = c(" ", A)
  B = c(" ", B)
  n = length(A)
  colNames = B
  rowNames = A
  matF = matrix(-Inf, n, n)
  row.names(matF) = rowNames
  colnames(matF) = colNames
  
  
  matF[1,1] = 0
  matF[1,2] = gap
  matF[2,1] = gap
  
  for(i in 2:n) {
    
    S = similarity2(B[i], A[i], match, mismatch)
    matF[i,i] = max(matF[i-1,i] + gap, matF[i,i-1] + gap, matF[i-1,i-1] + S)
    
    if(i+1 <= n) {
      S = similarity2(B[i+1], A[i], match, mismatch)
      matF[i+1,i] = max(matF[i,i] + gap, matF[i,i-1] + S)
      S = similarity2(B[i], A[i+1], match, mismatch)
      matF[i,i+1] = max(matF[i,i] + gap, matF[i-1,i] + S)
    }
  }
  return (matF)  
}



Needleman.WunschV2<-function(gap=-1,A,B)
{

  if (typeof(A) == "character" && length(A) == 1) {
    A = unlist(strsplit(A, split=""))
  }
  
  if (typeof(B) == "character" && length(B) == 1) {
    B = unlist(strsplit(B, split=""))
  }  
 
  match = 1
  mismatch = -1  
   
  Fij = ScoringV2(gap, A, B) 
  n<-length(A)
  m<-length(B)
  AlignmentA<-" "
  AlignmentB<-" "
  i<-n
  j<-m
  
  while(i>0||j>0)
  {
    if((i>0)&&(j>0)&&(Fij[i+1,j+1]==Fij[i,j]+similarity2(B[i], A[j], match, mismatch)))
    {
      AlignmentA<-trimws(paste0(B[i],AlignmentA))
      AlignmentB<-trimws(paste0(A[j],AlignmentB))
      i<-i-1
      j<-j-1
    }else if((i>0)&& (Fij[i+1,j+1]==Fij[i,j+1]+gap)){
      AlignmentA<-trimws(paste0(B[i],AlignmentA))
      AlignmentB<-trimws(paste0("-",AlignmentB))
      i<-i-1
    }else{
      AlignmentA<-trimws(paste0("-",AlignmentA))
      AlignmentB<-trimws(paste0(A[j],AlignmentB))
      j<-j-1
    }
  }
  # cat('\n')
  # print("F dans V2")
  # print(Fij)
  # cat("\n")
  return(list(AlignmentA=AlignmentA,AlignmentB=AlignmentB, score = Fij[n+1, n+1]))
  
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
  return (list(A=s, B=snew))  
}


getTimeExcecution = function(d, A, B, fun="V1") {
  t = 0
  if (fun == "V1") {
    t = system.time(Needleman.Wunsch(d,A,B))[[1]]
    # print(paste("tempps d'excécution première version:", toString(t1)))
  } else {
    t = system.time(Needleman.WunschV2(d,A,B))[[1]]
  }
  # print(paste("temps d'excécution deuxième version:", toString(t2)))
 return (t)
}
