
LOD_R2<-function(LOD,n)
{
  return(1-(10^((-2/n)*LOD)))
  
}

getP<-function(LOD,n,df,ddf)
{
  ff<-(10^((2/n)*LOD)-1)*((ddf)/df)
  pp<- -(pf(ff,df,ddf,lower.tail=FALSE,log=TRUE)/log(10))
  #pp<- pf(ff,df,ddf,lower.tail=FALSE)
  return(pp)
}


getF<-function(LOD,n,df)
{
  ddf<-n-(df+1)
  ff<-(10^((2/n)*LOD)-1)*((n-df-1)/df)
  pp<- -log10(pf(ff,df,ddf,lower.tail=FALSE))
  return(ff)
}

getLOD<-function(Fscore,n,df)
{
  (n/2)*log10((Fscore*(df/(n-df-1)))+1)
}


P_LOD<-function(P,n,ndf,ddf)
{
  Fscore<-qf(P,ndf,ddf,lower.tail=FALSE)
  return((n/2)*log10((Fscore*(ndf/(n-ndf-1)))+1))
}
