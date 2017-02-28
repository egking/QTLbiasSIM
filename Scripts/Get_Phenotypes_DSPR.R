#little function to get order for SNPs
decoder<-function(x,rils)
{
  which(rils==x)
}

getPhenos<-function(chr, mafs, ids,rils, All.loc,eff)
{
  
  exp.s<-subset(All.loc, CHROM==chr)
  exp.s<-exp.s[order(exp.s$POS),]
  
  
  mafs.exp<-mafs[which(mafs[,2] %in% exp.s$POS),1:3]
  all.equal(mafs.exp[,2],exp.s$POS)
  
  Aprobma<-do.call(rbind, strsplit(as.character(mafs.exp[,3]),"/"))
  
  Aprobma<-apply(Aprobma,2,function(x) as.numeric(x)/100)
  
  decodeA<-sapply(rils, function(x) decoder(x,ids[,1]))
  
  Aprobma<-Aprobma[,decodeA]
  
  genos<-apply(Aprobma, 2,function(x) rbinom(length(x),1,x))
  
  envs<-matrix(,nrow(genos),ncol(genos))
  
  
  for(i in 1:length(genos[,1]))
  {
    envs[i,]<-rnorm(length(genos[i,]),0,sqrt(((1/eff)-1)*var(genos[i,])))
  }
  
  phenos<-genos+envs
  phenos<-data.frame('CHROM'=exp.s$CHROM, 'POS'=exp.s$POS, 'effect'=exp.s$effect,phenos)
  colnames(phenos)[4:ncol(phenos)]<-rils
  
  return(phenos)
}

set.seed(351109)

library(DSPRqtl)
data(positionlist_wgenetic)

chromosomes<-c('X','2L','2R','3L','3R')

#get list of A rils
#filename<-paste("../GenomeCache/Release2/A/",poslist[1,1],"_",poslist[1,2],".rda",sep="")

#load(filename)


#pARILs <-as.integer(Ahmm_genos[,'ril'])

#save(pARILs, file="Data/DSPR_Release4/pARILs.rda")


#number of positions
npos<-10000

#get effects

#pull from dist
#exp.effects <- rexp(npos,10)
#unif.effects <- runif(npos,0.001,1)

#get snp table


load(file="Data/DSPR_Release4/pARILs.rda")

rilA<-read.table(file= "Data/DSPR_Release4/snpfreqs_RILS_A_X.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)

chromosomes<-c('2L','2R','3L','3R')

for(k in chromosomes)
{
  
  fname<-paste(file= "Data/DSPR_Release4/snpfreqs_RILS_A_",k,".txt",sep="")  
  ss<-read.table(fname,sep="\t",stringsAsFactors=FALSE,header=TRUE)  
  rilA<-rbind(rilA,ss) 
  
}

colnames(rilA)<-c("CHROM","POS","freqARILs")


for(k in chromosomes)
{
  minp<-min(poslist[poslist$chr==k,'Ppos'])
  maxp<-max(poslist[poslist$chr==k,'Ppos'])
  rilA<-rilA[-which(rilA$CHROM==k & rilA$POS>maxp),]
  rilA<-rilA[-which(rilA$CHROM==k & rilA$POS<minp),]
  
}

rilA<-rilA[-which(rilA$freqARILs>0.975 | rilA$freqARILs<0.025),]

######################set of effects
eff.all<-c(0.01,seq(0.025,0.25,by=0.025),0.5)

for(zz in 1:length(eff.all))
{
eff<-eff.all[zz]

#get random set

All.loc<-rilA[sample(seq(1,nrow(rilA)),npos),c('CHROM','POS')]



All.loc$effect<-eff
All.loc$pid<-seq(1,nrow(All.loc))

#generate phenos
#
chromosomes<-c('2L','2R','3L','3R')

idname<-paste('Data/DSPR_Release4/A_RILLIST_release4_','X','.txt',sep='')
ids<-read.table(idname,sep="\t",stringsAsFactors=FALSE,header=FALSE)

mafname<-paste('Data/DSPR_Release4/SNPtable_A_RILs_inferred_Release4_','X','.txt',sep='')
mafs<-read.table(mafname,sep="\t",stringsAsFactors=FALSE,header=FALSE)

pheno<-getPhenos("X",mafs,ids, pARILs,All.loc,eff)


for (k in chromosomes)
{
  
  idname<-paste('Data/DSPR_Release4/A_RILLIST_release4_',k,'.txt',sep='')
  ids<-read.table(idname,sep="\t",stringsAsFactors=FALSE,header=FALSE)
  
  mafname<-paste('Data/DSPR_Release4/SNPtable_A_RILs_inferred_Release4_',k,'.txt',sep='')
  mafs<-read.table(mafname,sep="\t",stringsAsFactors=FALSE,header=FALSE)
  
  pheno<-rbind(pheno,getPhenos(k,mafs,ids, pARILs,All.loc,eff))
  
}

pheno<-pheno[sample(seq(1,npos)),]

outf<-paste("Data/DSPR/SimPhenos/pheno_biallelic_",100*eff,".rda",sep="")
save(pheno,file=outf)
cat(eff, "\n")

}