

library(rrBLUP)
load(file="Data/dgrp2.rda")

dgrp.dat[1:10,1:10]
dgrp.dat[1:10,213:214]

maf<-dgrp.dat$refc/(dgrp.dat$refc+dgrp.dat$altc)
which.maf<-which(maf<=0.95 & maf>0.05)
dgrp.dat<-dgrp.dat[which.maf,]

dgrp.pos<-dgrp.dat[,1:2]
dgrp.dat<-dgrp.dat[,10:214]

dgrp.dat<-data.matrix(dgrp.dat)

nas<-apply(dgrp.dat,1,function(x) length(which(is.na(x))))
nas<-nas/ncol(dgrp.dat)
nas.c<-which(nas>=0.2)

dgrp.dat<-dgrp.dat[-nas.c,]
dgrp.pos<-dgrp.pos[-nas.c,]

dgrp.dat[dgrp.dat==2]<- 1
dgrp.dat[dgrp.dat==0]<- -1

load(file="Data/dropped_lines.rda")
dgrp.dat<-dgrp.dat[,-which(colnames(dgrp.dat) %in% elim.set)]

arms<-c('X','2L','2R','3L','3R')

kinall<-vector(mode='list',length=5)
names(kinall)<-arms


for(arm in arms)
{
  kinall[[arm]]<-A.mat(t(dgrp.dat[-which(dgrp.pos$chr==arm),]))
  cat(arm,"\n")
}

save(kinall,file="Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda")


