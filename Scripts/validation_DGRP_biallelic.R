set.seed(3221)

setwd("/home/kingeg/Projects/BeavisProj/")

source("Functions/hacked_DOQTL.R")

library(regress)

load(file="Data/DGRP/imputed_genos.rda")


load(file= "Data/DGRP/SimPhenos/DGRP_validate_biallelic.rda")

load(file="Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda")

#output<-vector(mode='list',length=nrow(phenos))
set<-phenos[,1:4]
set$maxLOD<-rep(NA,nrow(set))
set$maxPos<-rep(NA,nrow(set))
set$QLOD<-rep(NA,nrow(set))
counter<-1

for(i in 1:nrow(phenos))
{
line.set<-sample(colnames(phenos[,5:ncol(phenos)]),phenos$samps[i])
ccc<-phenos[i,'CHROM']
ppp<-phenos[i,'POS']

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
qpos<-phenos[i,'POS']

if(length(which(set.p[,2]==qpos))<=1)
{
  set$QLOD[counter]<-mod[which(set.p[,2]==qpos),1]
}else{
  set$QLOD[counter]<-mod[which(rownames(set.p)==rownames(phenos)[i]),1]
}


 counter<-counter+1
cat(i, "\t", Sys.time(),"\n")
}

save(set, file="Data/DGRP/Obs_LODs/Validate_LODs_biallelic_DGRP.rda")







