library(DSPRqtl)

data(positionlist_wgenetic)

big.list<-vector('list',nrow(poslist))
load(file="Data/DSPR/SimPhenos/pheno_biallelic_1.rda")


for (i in 1:nrow(poslist)) 
  
{  
  
  filename<-paste("../GenomeCache/Release2/A/",poslist[i,1],"_",poslist[i,2],".rda",sep="")
  
  load(filename)
  
  
  genotypes <-Ahmm_genos[,c("ril","AA1","AA2","AA3","AA4","AA5","AA6","AA7","AA8")]
  
  genos<-genotypes[order(genotypes$ril),c('ril','AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')]
  rr<-genos$ril
  genos<-as.matrix(genotypes[order(genotypes$ril),c('AA1','AA2','AA3','AA4','AA5','AA6','AA7','AA8')])
  
  rownames(genos)<-rr
  big.list[[i]]<-genos
  
  cat(i,"\n")
  
}

save(big.list,file="Data/DSPR/pA_biglist.rda")

all.equal(as.numeric(colnames(pheno)[4:ncol(pheno)]),as.numeric(genotypes$ril))


source("/Functions/convert_list_array.R")
big.array<-ListTo3dArray(big.list)
save(big.array,file="Data/DSPR/pA_bigarray.rda")

