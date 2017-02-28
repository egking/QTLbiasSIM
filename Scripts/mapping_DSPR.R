# Arguments: pheno: numeric vector with no NAs containing the phenotype.
#            probs: 3D numeric array of haplotype probabilities.
#            K: numeric matrix.
#fast.qtlrel.hap = function(pheno, probs, K, addcovar, snps) {


#single kinship A

args=(commandArgs(TRUE))

eff<-as.numeric(args[1])
samp.size<-as.numeric(args[2])
iseed<-as.numeric(args[3])
#eff<-0.2
#samp.size<-878

setwd("/home/kingeg/Projects/BeavisProj/")

source("Functions/hacked_DOQTL.R")

library(regress)

library(DSPRqtl)
data(positionlist_wgenetic)
load(file="Data/DSPR/pA_bigarray.rda")

load("Data/DSPR/Kin_LOCO_A.rda")

#all.equal(rownames(Phys.Out),dimnames(big.array)[[1]])



#mapall<-function(eff,samp.size)
#{
  
  fname<-paste("Data/DSPR/SimPhenos/pheno_biallelic_",100*eff,".rda",sep="")
  
  load(file=fname)
  
  pheno$CHROM<-as.character(pheno$CHROM)
  
  pheno<-pheno[1:1000,]
  
  lod.matrix<-vector('list',nrow(pheno))

    
for (j in 1:nrow(pheno))
  {
    #sample
    rilset<-sample(colnames(pheno[,4:ncol(pheno)]),samp.size)
    phenotype<-as.numeric(pheno[j,rilset])
    KK<-kinall[[pheno$CHROM[j]]][rilset,rilset]
    pp<-pheno$POS[j]
    ind<-which(poslist$chr==pheno$CHROM[j] & poslist$Ppos==(round(pp/1e4,0)*1e4))
    st<-ind-500
    if(st<0){st<-0}
    end<-ind+500
    if(end>nrow(poslist)){end<-nrow(poslist)}
    small.array<-big.array[rilset,,st:end]
    pos.set<-poslist[st:end,c('Ppos','chr','Ppos','Gpos')]
    names(pos.set)<-c('SNP_ID','CHROM','Mb','cM')
    pos.set$Mb<-pos.set$Mb/1e6
    lod.set<-fast.qtlrel.hap(pheno=phenotype, probs=small.array, K=KK,snps=pos.set)
    lod.matrix[[j]]<-lod.set[['lod']][,c('perc.var','lod','neg.log10.p')]
    cat(j, Sys.time(),"\n")
  }
  oname<-paste("Data/DSPR/Obs_LODs/LODseed_biallelic_KLOCO_E",100*eff,"_S",samp.size,".rda",sep="")
  
  save(lod.matrix,file=oname)

#}



