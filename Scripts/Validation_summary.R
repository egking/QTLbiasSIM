library(plyr)

setwd("/home/kingeg/Projects/BeavisProj/")

source("Functions/LOD_P_F_convert.R")

effs<-c(0.05,0.1,0.2)
sss<-c(100,185,500,878)

#DSPR.th<- -log10(0.00001)
#DGRP.th<- -log10(0.00001)
  
DGRP.th<-6.37

DSPR.th<-3.15

#DSPR

load(file="Data/DSPR/SimPhenos/DSPR_validate_biallelic.rda")

load(file="Data/DSPR/Obs_LODs/Validate_LODs_biallelic_DSPR.rda")
phenos$maxLOD<-unlist(lapply(lod.matrix, function(x) max(x$lod)))
phenos$R2max<-LOD_R2(phenos$maxLOD,phenos$samps)
phenos$P<-getP(phenos$maxLOD,phenos$samps,7,phenos$samps-8)
phenos$pow<-rep(0,nrow(phenos))
phenos$pow[phenos$P>=DSPR.th]<-1

dspr.b.v<-phenos[,c('CHROM','POS','effect','samps','maxLOD','R2max','P','pow')]
ddply(dspr.b.v, c("effect","samps"), function(dspr.b.v) c(mean(dspr.b.v$pow),nrow(dspr.b.v[dspr.b.v$pow==1,])))


load(file="Data/DSPR/SimPhenos/DSPR_validate_multi.rda")

load(file="Data/DSPR/Obs_LODs/Validate_LODs_multi_DSPR.rda")
phenos$maxLOD<-unlist(lapply(lod.matrix, function(x) max(x$lod)))
phenos$R2max<-LOD_R2(phenos$maxLOD,phenos$samps)
phenos$P<-getP(phenos$maxLOD,phenos$samps,7,phenos$samps-8)
phenos$pow<-rep(0,nrow(phenos))
phenos$pow[phenos$P>=DSPR.th]<-1

dspr.m.v<-phenos[,c('CHROM','POS1','POS2','POS3','effect','samps','maxLOD','R2max','P','pow')]
ddply(dspr.m.v, c("effect","samps"), function(dspr.m.v) c(mean(dspr.m.v$pow),nrow(dspr.m.v[dspr.m.v$pow==1,])))



#DGRP
load(file="Data/DGRP/Obs_LODs/Validate_LODs_biallelic_DGRP.rda")
set$R2max<-LOD_R2(set$maxLOD,set$samps)
set$P<-getP(set$maxLOD,set$samps,1,set$samps-2)
set$pow<-rep(0,nrow(set))
set$pow[set$P>=DGRP.th]<-1
dgrp.b.v<-set
ddply(dgrp.b.v, c("effect","samps"), function(dgrp.b.v) c(mean(dgrp.b.v$pow),nrow(dgrp.b.v[dgrp.b.v$pow==1,])))


load(file="Data/DGRP/Obs_LODs/Validate_LODs_multi_DGRP.rda")
set$R2max<-LOD_R2(set$maxLOD,set$samps)
set$P<-getP(set$maxLOD,set$samps,1,set$samps-2)
set$pow<-rep(0,nrow(set))
set$pow[set$P>=DGRP.th]<-1
dgrp.m.v<-set
ddply(dgrp.m.v, c("effect","samps"), function(dgrp.m.v) c(mean(dgrp.m.v$pow),nrow(dgrp.m.v[dgrp.m.v$pow==1,])))


########BETWEEN
effs<-c(0.05,0.1)
sss<-c(185,878)

#DSPR.th<- -log10(0.00001)

#DSPR.th<-4.45

#DSPR



pp.l<-vector(mode='list',length=2)

for(samps in 1:length(sss))
{
  load(file="Data/DSPR/SimPhenos/DSPR_validate_btwpop.rda")
  load(file=paste("Data/DSPR/Obs_LODs/Validate_LODs_btw_DSPR_",sss[samps],".rda",sep=""))
  phenos$maxLOD<-unlist(lapply(lod.matrix, function(x) max(x$lod)))
  phenos$R2max<-LOD_R2(phenos$maxLOD,sss[samps])
  phenos$P<-getP(phenos$maxLOD,sss[samps],7,sss[samps]-8)
  phenos$pow<-rep(0,nrow(phenos))
  phenos$pow[phenos$P>=DSPR.th]<-1
  phenos$samps<-sss[samps]
  pp.l[[samps]]<-phenos
}

dspr.btw.v<-rbind(pp.l[[1]][,c('CHROM','POS','effect','samps','maxLOD','R2max','P','pow')],
                  pp.l[[2]][,c('CHROM','POS','effect','samps','maxLOD','R2max','P','pow')])
ddply(dspr.btw.v, c("effect","samps"), function(dspr.btw.v) c(mean(dspr.btw.v$pow),nrow(dspr.btw.v[dspr.btw.v$pow==1,])))



load(file="Data/DGRP/Obs_LODs/Validate_LODs_btw_DGRP_185.rda")
set$samps<-185
set$R2max<-LOD_R2(set$maxLOD,185)
set$P<-getP(set$maxLOD,185,1,185-2)
set$pow<-rep(0,nrow(set))
set$pow[set$P>=DGRP.th]<-1
dgrp.btw.v<-set
ddply(dgrp.btw.v, c("effect","samps"), function(dgrp.btw.v) c(mean(dgrp.btw.v$pow),nrow(dgrp.btw.v[dgrp.btw.v$pow==1,])))



#ALL
load(file="Data/DGRP/Obs_LODs/map_combined_LOCO_biallelic.rda")
out.all$R2max<-LOD_R2(out.all$maxLOD,out.all$samps)
out.all$R2Q<-LOD_R2(out.all$QLOD,out.all$samps)
out.all$P<-getP(out.all$maxLOD,out.all$samps,1,out.all$samps-2)
out.all$pow<-rep(0,nrow(out.all))
out.all$pow[out.all$P>=DGRP.th]<-1
out.all$type<-'b'
dgrp.comb<-out.all
rm(out.all)

load(file="Data/DGRP/Obs_LODs/map_combined_LOCO_multi.rda")
out.all$R2max<-LOD_R2(out.all$maxLOD,out.all$samps)
out.all$R2Q<-LOD_R2(out.all$QLOD,out.all$samps)
out.all$P<-getP(out.all$maxLOD,out.all$samps,1,out.all$samps-2)
out.all$pow<-rep(0,nrow(out.all))
out.all$pow[out.all$P>=DGRP.th]<-1
out.all$type<-'m'
out.all$POS<-rowMeans(out.all[,c('POS1','POS2','POS3')])
out.all<-out.all[,colnames(dgrp.comb)]
dgrp.comb<-rbind(dgrp.comb, out.all)

load(file="Data/DSPR/Obs_LODs/map_combined_LOCOseed.rda")
dspr.all$P<-getP(dspr.all$maxLOD,dspr.all$samps,7,dspr.all$samps-8)
dspr.all$pow<-rep(0,nrow(dspr.all))
dspr.all$pow[dspr.all$P>=DSPR.th]<-1
dspr.all$type<-'b'
dspr.comb<-dspr.all

load(file="Data/DSPR/Obs_LODs/map_combined_multiLOCO.rda")
dspr.all$P<-getP(dspr.all$maxLOD,dspr.all$samps,7,dspr.all$samps-8)
dspr.all$pow<-rep(0,nrow(dspr.all))
dspr.all$pow[dspr.all$P>=DSPR.th]<-1
dspr.all$type<-'m'
dspr.comb<-rbind(dspr.comb,dspr.all)


types<-c('b','m')
effs<-c(0.05,0.1,0.2)
sss<-c(100,185,500,878)

dspr.sum<-data.frame("type"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "eff"=numeric(length=(length(sss)*length(effs)*length(types))), 
                     "samps"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "pow"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "Meff"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "SEeff"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "N"=numeric(length=(length(sss)*length(effs)*length(types))))
counter<-1
for(tt in types)
{
  for(j in effs)
  {
    for(k in sss)
    {
      dd<-subset(dspr.comb, effect==j & samps==k & type==tt)
      dspr.sum$type[counter]<-tt
      dspr.sum$eff[counter]<-j
      dspr.sum$samps[counter]<-k
      dspr.sum$pow[counter]<-mean(dd$pow)
      dspr.sum$Meff[counter]<-mean(dd[dd$pow==1,'R2max'])
      dspr.sum$SEeff[counter]<-sd(dd[dd$pow==1,'R2max'])/length(dd[dd$pow==1,'R2max'])
      dspr.sum$N[counter]<-length(dd[dd$pow==1,'R2max'])
      counter<-counter+1  
    }
  }
}


types<-c('b','m')
effs<-c(0.05,0.1,0.2)
sss<-c(100,185)

dgrp.sum<-data.frame("type"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "eff"=numeric(length=(length(sss)*length(effs)*length(types))), 
                     "samps"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "pow"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "Meff"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "SEeff"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "N"=numeric(length=(length(sss)*length(effs)*length(types))))
counter<-1
for(tt in types)
{
  for(j in effs)
  {
    for(k in sss)
    {
      dd<-subset(dgrp.comb, effect==j & samps==k & type==tt)
      dgrp.sum$type[counter]<-tt
      dgrp.sum$eff[counter]<-j
      dgrp.sum$samps[counter]<-k
      dgrp.sum$pow[counter]<-mean(dd$pow)
      dgrp.sum$Meff[counter]<-mean(dd[dd$pow==1,'R2max'])
      dgrp.sum$SEeff[counter]<-sd(dd[dd$pow==1,'R2max'])/length(dd[dd$pow==1,'R2max'])
      dgrp.sum$N[counter]<-length(dd[dd$pow==1,'R2max'])
      counter<-counter+1  
    }
  }
}




