#random set
#get random regions
#get max in region
#get FWER
setwd("/home/kingeg/Projects/BeavisProj/")
source("Functions/LOD_P_F_convert.R")

#DGRP
load(file="Data/DGRP/imputed_genos.rda")
ffs<-list.files("Data/DGRP/Perm_LODs/")
maxP<-numeric(length=length(ffs))

for(ii in 1:4000)
{
  #time1<-Sys.time()
  tryCatch(
    load(paste("Data/DGRP/Perm_LODs/",ffs[ii],sep="")),
    error=function(e) cat(ii))
  if(ii<2001)
  {
    pp<-getP(ll, 100,1,98)
    r.ii<-sample(seq(1, length(pp)),1)
    ccc<-imputed.genos$chr[r.ii]
    ppp<-imputed.genos$pos[r.ii]
    set.p<-pp[(imputed.genos$chr == ccc & 
                            imputed.genos$pos > (ppp-5000000) &
                            imputed.genos$pos < (ppp+5000000))]
    maxP[ii]<-max(set.p)
    
  }else{
    pp<-getP(ll, 185,1,183) 
    r.ii<-sample(seq(1, length(pp)),1)
    ccc<-imputed.genos$chr[r.ii]
    ppp<-imputed.genos$pos[r.ii]
    set.p<-pp[(imputed.genos$chr == ccc & 
                 imputed.genos$pos > (ppp-5000000) &
                 imputed.genos$pos < (ppp+5000000))]
    maxP[ii]<-max(set.p)
    
  }
  cat("1\t",ii, "\n")
}

save(maxP, file="Data/DGRP/FWER_region_DGRP.rda")  

#DSPR
ll<-list.files('Data/DSPR/Perm_LODs/')
maxPs<-numeric(length=length(ll))
  
for(file.i in 1:length(ll))
{
  load(paste('Data/DSPR/Perm_LODs/',ll[file.i],sep=''))
  if(file.i<2001)
  {
    p.set<-getP(LOD.SET$lod,100,7,92)
    ii<-sample(seq(1,length(p.set)),1)
    st<-ii-500
    if(st<0){st<-0}
    end<-ii+500
    if(end>length(p.set)){end<-length(p.set)}
    maxPs[file.i]<-max(p.set[st:end])
  }else{
    p.set<-getP(LOD.SET$lod,878,7,870)
    ii<-sample(seq(1,length(p.set)),1)
    st<-ii-500
    if(st<0){st<-0}
    end<-ii+500
    if(end>length(p.set)){end<-length(p.set)}
    maxPs[file.i]<-max(p.set[st:end])
  }
  cat("2\t",file.i, "\n")
}

save(maxPs, file="Data/DSPR/FWER_region_DSPR.rda")  



load(file="Data/DSPR/FWER_region_DSPR.rda")
quantile(maxPs, 0.95)

load(file="Data/DGRP/FWER_region_DGRP.rda")
quantile(maxP, 0.95,na.rm=TRUE)
