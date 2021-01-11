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
  return(Fij)
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
  return(list(AlignmentA=AlignmentA,AlignmentB=AlignmentB))
  
}
