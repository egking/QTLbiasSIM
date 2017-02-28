set.seed(3221)

setwd("/home/kingeg/Projects/BeavisProj/")

source("Functions/hacked_DOQTL.R")

library(regress)

load(file="Data/DGRP/imputed_genos.rda")


load(file= "Data/DGRP/SimPhenos/DGRP_validate_multi.rda")

load(file="Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda")

#output<-vector(mode='list',length=nrow(phenos))
set<-phenos[,1:6]
set$maxLOD<-rep(NA,nrow(set))
set$maxPos<-rep(NA,nrow(set))
set$QLOD<-rep(NA,nrow(set))
counter<-1

for(i in 1:nrow(phenos))
{
line.set<-sample(colnames(phenos[,7:ncol(phenos)]),phenos$samps[i])
ccc<-phenos[i,'CHROM']
ppp<-mean(as.numeric(phenos[i,c('POS1','POS2','POS3')]))

set.p<-imputed.genos[(imputed.genos$chr == ccc & 
                       imputed.genos$pos > (ppp-5000000) &
                       imputed.genos$pos < (ppp+5000000)),c('chr','pos',line.set)]

#named phenotype numeric
#named kinship
#lines by positions

#all.equal(colnames(all.phenos)[4:ncol(all.phenos)],colnames(kinmat_cut))
kinmat<-kinall[[ccc]]
kinmat<-kinmat/2
kinmat<-kinmat[line.set,line.set]
phen.i<-t(phenos[i,line.set])
mod<-snp.scan(pheno=phen.i, 
                  K=kinmat, 
                  geno=t(set.p[,3:ncol(set.p)]))
out<-set.p[,1:2] 
out$LOD<-mod[,1]
set$maxLOD[counter]<-max(mod[,1])
set$maxPos[counter]<-set.p[which.max(mod[,1]),2]
qpos<-phenos[i,'POS1']

if(length(which(set.p[,2]==qpos))<=1)
{
  set$QLOD[counter]<-mod[which(set.p[,2]==qpos),1]
}else{
  set$QLOD[counter]<-mod[which(set.p[,2]==qpos)[1],1]
}


 counter<-counter+1
cat(i, "\t", Sys.time(),"\n")
}

save(set, file="Data/DGRP/Obs_LODs/Validate_LODs_multi_DGRP.rda")







