proc.res<-function(eff,samp.size)
{
library(DSPRqtl)
data(positionlist_wgenetic)

fname<-paste("Data/DSPR/SimPhenos/pheno_biallelic_",100*eff,".rda",sep="")
load(file=fname)
pheno$CHROM<-as.character(pheno$CHROM)

oname<-paste("Data/DSPR/Obs_LODs/LODseed_biallelic_KLOCO_E",100*eff,"_S",samp.size,".rda",sep="")
load(file=oname)

dspr.out<-pheno[c(1:1000),1:3]
dspr.out$samps<-rep(samp.size,nrow(dspr.out))

dspr.out$maxLOD<-unlist(lapply(lod.matrix,function(x) max(x$lod)))
dspr.out$maxPos<-rep(NA,nrow(dspr.out))
dspr.out$QLOD<-rep(NA,nrow(dspr.out))

for(ii in 1:1000)
{
  ind<-which(poslist$chr==pheno$CHROM[ii] & poslist$Ppos==(round(pheno$POS[ii]/1e4,0)*1e4))
  st<-ind-500
  if(st<0){st<-0}
  end<-ind+500
  if(end>nrow(poslist)){end<-nrow(poslist)}
  p.small<-poslist[st:end,]
  dspr.out$maxPos[ii]<-p.small[which.max(lod.matrix[[ii]]$lod),'Ppos']
  dspr.out$QLOD[ii]<-lod.matrix[[ii]]$lod[which(p.small$Ppos==round(pheno$POS[ii]/1e4,0)*1e4)]
  
}
return(dspr.out)
}

proc.res.ex<-function(eff,samp.size)
{
  library(DSPRqtl)
  data(positionlist_wgenetic)
  
  fname<-paste("Data/DSPR/SimPhenos/pheno_biallelic_",100*eff,".rda",sep="")
  load(file=fname)
  pheno$CHROM<-as.character(pheno$CHROM)
  
  oname<-paste("Data/DSPR/Obs_LODs/LODseed_biallelic_KLOCO_E",100*eff,"_S",samp.size,"_2.rda",sep="")
  load(file=oname)
  
  dspr.out<-pheno[c(1001:2000),1:3]
  dspr.out$samps<-rep(samp.size,nrow(dspr.out))
  
  dspr.out$maxLOD<-unlist(lapply(lod.matrix,function(x) max(x$lod)))
  dspr.out$maxPos<-rep(NA,nrow(dspr.out))
  dspr.out$QLOD<-rep(NA,nrow(dspr.out))
  cc<-1
  for(ii in 1001:2000)
  {
    ind<-which(poslist$chr==pheno$CHROM[ii] & poslist$Ppos==(round(pheno$POS[ii]/1e4,0)*1e4))
    st<-ind-500
    if(st<0){st<-0}
    end<-ind+500
    if(end>nrow(poslist)){end<-nrow(poslist)}
    p.small<-poslist[st:end,]
    dspr.out$maxPos[cc]<-p.small[which.max(lod.matrix[[cc]]$lod),'Ppos']
    dspr.out$QLOD[cc]<-lod.matrix[[cc]]$lod[which(p.small$Ppos==round(pheno$POS[ii]/1e4,0)*1e4)]
    cc<-cc+1
  }
  return(dspr.out)
}
effs<-c(0.05,0.1,0.2)
sss<-c(100,185,500,878)

for(j in effs)
{
  for(k in sss)
  {
    oset<-proc.res(j,k)
    if(j==effs[1] & k==sss[1])
    {
      dspr.all<-oset
    }else{
      dspr.all<-rbind(dspr.all,oset)
    }
    cat(j,"\t",k,"\n")
  }
}

dex<-proc.res.ex(0.05,100)
dspr.all<-rbind(dspr.all, dex)


source("Functions/LOD_P_F_convert.R")

dspr.all$R2max <- rep(NA, nrow(dspr.all))
dspr.all$R2Q <- rep(NA, nrow(dspr.all))

for(i in 1:nrow(dspr.all))
{
  dspr.all$R2max[i]<-LOD_R2(dspr.all$maxLOD[i],dspr.all$samps[i])
  dspr.all$R2Q[i]<-LOD_R2(dspr.all$QLOD[i],dspr.all$samps[i])
#  cat(i)
}


save(dspr.all, file="Data/DSPR/Obs_LODs/map_combined_LOCOseed.rda")




