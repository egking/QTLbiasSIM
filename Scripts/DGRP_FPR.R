pfind<-function(pp, poslist, th, tol.dist)
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


convG<-function(chr,pos,genmap)
{
  library(pspline)
  map<-genmap[[chr]]
  Mpos<-as.numeric(pos)/1000000
  Gpos<-as.numeric(predict(map,Mpos))
  return(Gpos) 
  
}
args=(commandArgs(TRUE))

tt<-as.numeric(args[1])

load(file='Data/geneticmapLIST.rda')

load(file="Data/DGRP/imputed_genos.rda")

arms<-c('X','2L','2R','3L','3R')

for(kk in arms)
{
  dgrp.p<- imputed.genos[imputed.genos$chr==kk,c('chr','pos')]
  
  dgrp.p$Gpos<-convG(kk, dgrp.p$pos, genmap)
  if(kk =='X')
  {
  dgrp.pos<-dgrp.p  
  }else{
  dgrp.pos<-rbind(dgrp.pos, dgrp.p)  
  }
}  
  

source("Functions/LOD_P_F_convert.R")

ffs<-list.files("Data/DGRP/Perm_LODs/")

#tt<-seq(3.6,7,by=0.2)

#oo<-matrix(NA,4000,length(tt))
oo<-numeric(length=4000)

for(ii in 1:4000)
{
  #time1<-Sys.time()
  tryCatch(
    load(paste("Data/DGRP/Perm_LODs/",ffs[ii],sep="")),
    error=function(e) cat(ii))
  if(ii<2001)
  {
    pp<-getP(ll, 100,1,98)
  }else{
    pp<-getP(ll, 185,1,183) 
  }
  #for(th in 1:length(tt))
  #{
    oo[ii] <- pfind(pp, dgrp.pos, tt,0.5)
  #}
 
   cat(ii, "\n")
}

save(oo, file=paste('Data/DGRP/FPR_rates',tt,'.rda',sep=''))

#nn<-getP(ll, 185,1,183)




 tt<-seq(3.75,7.75,by=0.25)
 all.n<-matrix(NA,4000,length(tt))
 
 for(ifile in 1:length(tt))
 {
   load(paste("Data/DGRP/Perm_LODs/",ffs[ii],sep="")),
   all.n[,ifile]<-oo
 }
 
 save(all.n, file='Data/DGRP/FPR_rates_ALL.rda')    
#   
# 
# 
# 

 

