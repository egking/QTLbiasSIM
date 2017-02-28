pfind<-function(pp, cM, th, tol.dist)
{
  
  d1<-c(diff(pp),0)
  d2<- c(0,-d1[-length(d1)])
  
  init.p<- poslist[which(pp > th & d1 <= 0 & d2 <= 0),]
  init.p$LL<-pp[which(pp > th & d1 <= 0 & d2 <= 0)]
  init.p<-init.p[order(-init.p$LL),]
  #if none - return 0
  cc<-1
  for(kk in 1:nrow(init.p))
  {
    if(cc<=nrow(init.p))
    {
      ind.e<-which(init.p[,'chr']==init.p[cc,'chr'] & abs(init.p[,'Gpos'] - init.p[cc, 'Gpos']) < tol.dist)
      ind.e<-ind.e[-which(ind.e==cc)]
      if(length(ind.e)>0)
      {
        init.p<-init.p[-ind.e,,drop=FALSE]
      }
      cc<-cc+1
    }
  }
  return(nrow(init.p))
  
}

source("Functions/LOD_P_F_convert.R")
library(DSPRqtl)
data(positionlist_wgenetic)

#match DGRP
#thresh.lod<-seq(3.6,7,by=0.2)
#pp.th<-getP(thresh.lod, 185,1,183)

#P_LOD<-function(P,n,ndf,ddf)
#P_LOD(10^-pp.th,878,7,870)

th.set<-seq(1,6, by=0.25)

ll<-list.files('Data/DSPR/Perm_LODs/')
L.hit<-matrix(NA, length(ll),length(th.set))
  counter<-1
  for(file.i in 1:4000)
  {
    load(paste('Data/DSPR/Perm_LODs/',ll[file.i],sep=''))
    if(file.i<2001)
    {
      p.set<-getP(LOD.SET$lod,100,7,92)
    }else{
      p.set<-getP(LOD.SET$lod,878,7,870)
    }
    
      for(tt in 1:length(th.set))
    {
      
    L.hit[counter,tt]<- pfind(p.set,cM=poslist[,c('chr','Gpos')] ,th=th.set[tt], tol.dist=2)
    
    }
    counter<-counter+1 
    cat(counter, "\n")
  }
  
  
save(L.hit, file='Data/DSPR/FPR_rates_ALL.rda')


maxP<-numeric(4000)

for(file.i in 1:4000)
{
  load(paste('Data/DSPR/Perm_LODs/',ll[file.i],sep=''))
  if(file.i<2001)
  {
    p.set<-getP(LOD.SET$lod,100,7,92)
  }else{
    p.set<-getP(LOD.SET$lod,878,7,870)
  }
  maxP[file.i]<-max(p.set)
}

save(maxP, file='Data/DSPR/FWER.rda')
