source("Functions/LOD_P_F_convert.R")

ffs<-list.files("Data/DGRP/Perm_LODs/")




maxP<-numeric(4000)

for(file.i in 1:4000)
{
  load(paste("Data/DGRP/Perm_LODs/",ffs[file.i],sep=""))
  if(file.i<2001)
  {
    p.set<-getP(ll,100,1,98)
  }else{
    p.set<-getP(ll,185,1,183)
  }
  maxP[file.i]<-max(p.set)
  cat(file.i,"\n")
}

save(maxP, file="Data/DGRP/FWER.rda")


load(file="Data/DGRP/FWER.rda")
quantile(maxP[1:2000],0.95)
quantile(maxP[2001:4000],0.95)
