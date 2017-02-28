args=(commandArgs(TRUE))

st.ind<-as.numeric(args[1])

samp.size<-as.numeric(args[2]) 
iseed<-as.numeric(args[3])

load(file= "Data/randNum.rda")
set.seed(randNum[iseed])

rm(randNum)

en.ind<-st.ind+199

source("Functions/hacked_DOQTL.R")

library(regress)

load(file="Data/DGRP/imputed_genos.rda")

load(file="Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda")

arms<-c('X','2L','2R','3L','3R')

#all.lods<-matrix(NA,nrow(imputed.genos[imputed.genos$chr %in% arms,]),en.ind-(st.ind-1))
  #named phenotype numeric
  #named kinship
  #lines by positions
counter<-1
for(i in st.ind:en.ind)
{
  time1<-Sys.time()
  line.set<-sample(colnames(imputed.genos[,3:ncol(imputed.genos)]),samp.size)
  for(cc in arms)
  {
  ii<-imputed.genos[imputed.genos$chr==cc,line.set]
  kinmat<-kinall[[cc]]
  kinmat<-kinmat/2
  kinmat<-kinmat[line.set,line.set]
  mod.s<-snp.scan(pheno=rnorm(samp.size), 
                K=kinmat, 
                geno=t(ii))
  if(cc=='X')
  {
    mod<-mod.s
  }else{
    mod<-rbind(mod,mod.s)
  }
  rm(mod.s,kinmat,ii)
  gc()
  }
  
  #all.equal(colnames(all.phenos)[4:ncol(all.phenos)],colnames(kinmat_cut))
 
  #all.lods[,counter]<-mod[,1]
  counter<-counter+1
  time2<-Sys.time()
cat(counter, "\t",time2-time1,"\n")
ll<-mod[,1]
save(ll, file=paste("Data/DGRP/Perm_LODs/All_lods_randLOCO_",samp.size,"_",i,".rda",sep=""))
rm(mod, ll)
gc()
}



# 
# st.inds<-c(1,501,1001,1501)
# for(i in st.inds)
# {
# load(file=paste("OutputFiles/maxlods_rand_",i,".rda",sep=""))
#   max.lods<-max.lods[i:(i+499)]
#   if(i==1)
#   {
#     max.all<-max.lods
#   }else{
#     max.all<-c(max.all,max.lods)
#   }
#   
# }
# quantile(max.all,0.95)
# 
# 
# LODnom<-P_LOD(10^-5,185,1,183)
# length(which(max.all>=LODnom))


  
