args=(commandArgs(TRUE))

st.ind<-as.numeric(args[1])

en.ind<-st.ind+999

samp.size<-as.numeric(args[2]) 
iseed<-as.numeric(args[3])



source("Functions/hacked_DOQTL.R")

library(regress)

load(file= "Data/randNum.rda")
set.seed(randNum[iseed])


load(file="Data/DGRP/imputed_genos.rda")

load(file= "Data/DGRP/SimPhenos/allphenos_noinv_multi.rda")

load(file="Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda")

output<-vector(mode='list',length=en.ind-(st.ind-1))

set<-all.phenos[st.ind:en.ind,1:5]
set$maxLOD<-rep(NA,nrow(set))
set$maxPos<-rep(NA,nrow(set))
set$QLOD<-rep(NA,nrow(set))

counter<-1
for(i in st.ind:en.ind)
{
line.set<-sample(colnames(all.phenos[,6:ncol(all.phenos)]),samp.size)
ccc<-all.phenos[i,'CHROM']

ppp<-mean(as.numeric(all.phenos[i,c('POS1','POS2','POS3')]))

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
phen.i<-t(all.phenos[i,line.set])
mod<-snp.scan(pheno=phen.i, 
                  K=kinmat, 
                  geno=t(set.p[,3:ncol(set.p)]))
out<-set.p[,1:2] 
out$LOD<-mod[,1]
set$maxLOD[counter]<-max(mod[,1])
set$maxPos[counter]<-set.p[which.max(mod[,1]),2]

qpos<-all.phenos[i,'POS1']

if(length(which(set.p[,2]==qpos))<=1)
{
  set$QLOD[counter]<-mod[which(set.p[,2]==qpos),1]
}else{
  set$QLOD[counter]<-mod[which(set.p[,2]==qpos)[1],1]
}


output[[counter]]<-out
 counter<-counter+1
cat(i, "\t", Sys.time(),"\n")
}



#fname<-paste("Data/DGRP/Obs_LODs/DGRPscansLOCOmulti_",samp.size,"_",st.ind,"_",en.ind,".rda",sep="")
#save(output, file=fname)

fname2<-paste("Data/DGRP/Obs_LODs/DGRPsetLOCOmulti_",samp.size,"_",st.ind,"_",en.ind,".rda",sep="")
save(set, file=fname2)







