# Arguments: pheno: numeric vector with no NAs containing the phenotype.
#            probs: 3D numeric array of haplotype probabilities.
#            K: numeric matrix.
#fast.qtlrel.hap = function(pheno, probs, K, addcovar, snps) {


#single kinship A

set.seed(4215)

setwd("/home/kingeg/Projects/BeavisProj/")

source("Functions/hacked_DOQTL.R")

library(regress)

library(DSPRqtl)
data(positionlist_wgenetic)

load(file="Data/DSPR/pA_bigarray.rda")

load("Data/DSPR/Kin_LOCO_A.rda")

load(file="Data/DSPR/SimPhenos/DSPR_validate_multi.rda")


  lod.matrix<-vector('list',nrow(phenos))

    
for (j in 1:nrow(phenos))
  {
    #sample
    rilset<-sample(colnames(phenos[,7:ncol(phenos)]),phenos$samps[j])
    phenotype<-as.numeric(phenos[j,rilset])
    KK<-kinall[[phenos$CHROM[j]]][rilset,rilset]
    pp<-as.integer(mean(as.numeric(phenos[j,2:4])))
    ind<-which(poslist$chr==phenos$CHROM[j] & poslist$Ppos==(round(pp/1e4,0)*1e4))
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
  save(lod.matrix, file="Data/DSPR/Obs_LODs/Validate_LODs_multi_DSPR.rda")
  