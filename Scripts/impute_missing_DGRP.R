library(rrBLUP)
load(file="Data/DGRP/dgrp2.rda")

load(file="Data/dropped_lines.rda")

elim.n<-which(colnames(dgrp.dat) %in% elim.set)
dgrp.dat<-dgrp.dat[,-elim.n]

maf<-dgrp.dat$refc/(dgrp.dat$refc+dgrp.dat$altc)
which.maf<-which(maf<=0.975 & maf>0.025)
dgrp.dat<-dgrp.dat[which.maf,]

dgrp.info<-dgrp.dat[,1:9]
dgrp.dat<-dgrp.dat[,10:ncol(dgrp.dat)]

dgrp.dat<-data.matrix(dgrp.dat)

nas<-apply(dgrp.dat,1,function(x) length(which(is.na(x))))
nas<-nas/ncol(dgrp.dat)
nas.c<-which(nas>=0.2)

dgrp.dat<-dgrp.dat[-nas.c,]
dgrp.info<-dgrp.info[-nas.c,]

dgrp.dat[dgrp.dat==2]<- 1
dgrp.dat[dgrp.dat==0]<- -1

tester<-A.mat(t(dgrp.dat),impute.method="EM",return.imputed=TRUE)

imputed.genos<-data.frame('chr'=dgrp.info$chr, 'pos'=dgrp.info$pos, stringsAsFactors=FALSE)
imputed.genos<-cbind(imputed.genos, t(tester$imputed))

rbind(dgrp.dat[1,],
imputed.genos[1,])

save(imputed.genos, file="Data/DGRP/imputed_genos.rda")

