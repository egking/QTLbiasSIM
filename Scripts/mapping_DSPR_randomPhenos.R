#run this before below
#set.seed(90)
#info.mat<-cbind(seq(1,4000),c(rep(878,2000),rep(100,2000)),as.integer(runif(4000,1,100000)))
#save(info.mat,file='InputFiles/perms_DSPR_info.rda')



args=(commandArgs(TRUE))

ii<-as.numeric(args[1])

load(file='Data/DSPR/perms_DSPR_info.rda')

set.seed(info.mat[ii,3])
samp.size<-info.mat[ii,2]

source("Functions/hacked_DOQTL.R")

library(regress)

library(DSPRqtl)
data(positionlist_wgenetic)
load(file="Data/DSPR/pA_bigarray.rda")

load("Data/DSPR/Kin_LOCO_A.rda")

    
for(kk in c('X','2L','2R','3L','3R'))
    {
    
    rilset<-sample(rownames(big.array),samp.size)
    KK<-kinall[[kk]][rilset,rilset]
    pos.set<-poslist[poslist$chr==kk,c('Ppos','chr','Ppos','Gpos')]
    names(pos.set)<-c('SNP_ID','CHROM','Mb','cM')
    pos.set$Mb<-pos.set$Mb/1e6
    small.array<-big.array[rilset,,which(poslist$chr==kk)]
    lod.set<-fast.qtlrel.hap(pheno=rnorm(samp.size), probs=small.array, K=KK,snps=pos.set)
    if(kk=='X')
    {
    LOD.SET<-lod.set$lod  
    }else{
    LOD.SET<-rbind(LOD.SET,lod.set$lod)
    }
    
    }

 
save(LOD.SET, file=paste("Data/DSPR/Perm_LODs/LODSET_",samp.size,"_",ii,".rda",sep=""))






