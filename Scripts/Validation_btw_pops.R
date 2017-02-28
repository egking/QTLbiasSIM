setwd("/home/kingeg/Projects/BeavisProj/")

#############FUNCTIONS

#little function to get order for SNPs
decoder<-function(x,rils)
{
  which(rils==x)
}

getGvar<-function(chr, pos, eff, mafs,ids,rils)
{
  mafs.pos<-mafs[[chr]][which(mafs[[chr]][,2]== pos),1:3]
  
  Aprobma<-do.call(rbind, strsplit(as.character(mafs.pos[,3]),"/"))
  
  Aprobma<-apply(Aprobma,2,function(x) as.numeric(x)/100)
  
  decodeA<-sapply(rils, function(x) decoder(x,ids[[chr]][,1]))
  
  Aprobma<-Aprobma[decodeA]
  
  genos<- rbinom(length(Aprobma),1,Aprobma)
  
  
  return(var(genos))
}

getPhenos_dgrp<-function(chr, pos, eff,gvar, mafs,ids,rils)
{
  mafs.pos<-mafs[[chr]][which(mafs[[chr]][,2]== pos),1:3]
  
  Aprobma<-do.call(rbind, strsplit(as.character(mafs.pos[,3]),"/"))
  
  Aprobma<-apply(Aprobma,2,function(x) as.numeric(x)/100)
  
  decodeA<-sapply(rils, function(x) decoder(x,ids[[chr]][,1]))
  
  Aprobma<-Aprobma[decodeA]
  
  genos<- rbinom(length(Aprobma),1,Aprobma)
  
  envs<-rnorm(length(genos),0,sqrt(((1/eff)-1)*gvar))
  
  phenos<-genos+envs
  
  names(phenos)<-rils
  
  return(phenos)
}



#############FUNCTIONS
#############FUNCTIONS
#############FUNCTIONS
source("Functions/LOD_P_F_convert.R")



library(DSPRqtl)

data(positionlist_wgenetic)



#make new phenotypes for all
set.seed(56329)


load(file= 'Data/DSPR/Obs_LODs/All_Obs_LODs.rda')

#biallelic
dspr.comb<-dspr.comb[dspr.comb$type=='b',]


effs<-c(0.05,0.1)
sss<-c(185,878)

dspr.pos<-dspr.comb[0,c('CHROM','POS','effect','samps')]

for(ee in effs)
{
  for(ss in sss)
  {
      
      dspr.pos<-rbind(dspr.pos, dspr.comb[dspr.comb$pow==1 & dspr.comb$effect==ee & dspr.comb$samps==ss,c('CHROM','POS','effect','samps')])
  
  }
}



load(file="Data/DGRP/imputed_genos.rda")

ff<-as.matrix(imputed.genos[,3:187])
ff[ff< -1]<- -1
ff[ff> 1]<- 1

ff<-(ff+1)/2

imputed.genos<-imputed.genos[,c(1,2)]
imputed.genos<-cbind(ii, ff)


#get gen var in dspr
chromosomes<-c('X','2L','2R','3L','3R')
load(file="Data/DSPR_Release4/pARILs.rda")

ids<-vector(mode='list', length=5)
mafs<-vector(mode='list', length=5)
for(ii in 1:5)
{
  idname<-paste('Data/DSPR_Release4/A_RILLIST_release4_',chromosomes[ii],'.txt',sep='')
  ids[[ii]]<-read.table(idname,sep="\t",stringsAsFactors=FALSE,header=FALSE)
  
  mafname<-paste('Data/DSPR_Release4/SNPtable_A_RILs_inferred_Release4_',chromosomes[ii],'.txt',sep='')
  mafs[[ii]]<-read.table(mafname,sep="\t",stringsAsFactors=FALSE,header=FALSE)  
  
}
names(mafs)<-chromosomes
names(ids)<-chromosomes

dspr.gvar<-t(apply(dspr.pos, 1, function(x) getGvar(chr=as.character(x[1]), pos=as.numeric(x[2]), eff=as.numeric(x[3]),mafs=mafs, ids=ids, rils=pARILs)))
dspr.pos$gvar<-as.numeric(dspr.gvar)

pos.all<-merge(imputed.genos, dspr.pos, by.x=c('chr','pos'),by.y=c('CHROM','POS'))

(nrow(dspr.pos)-nrow(pos.all))/nrow(dspr.pos)

genos<-data.matrix(pos.all[,3:187])
envs<-matrix(,nrow(genos),ncol(genos))
for(i in 1:length(genos[,1]))
{
  envs[i,]<-rnorm(length(genos[i,]),0,sqrt(((1/pos.all[i,'effect'])-1)*pos.all[i,'gvar']))
}
phenos<-genos+envs

phenos<-data.frame('CHROM'=pos.all$chr, 'POS'=pos.all$pos, 'effect'=pos.all$effect, 'samps'=pos.all$samps,phenos,stringsAsFactors=FALSE)

save(phenos, file="Data/DGRP/SimPhenos/DGRP_validate_btwpop.rda")


###DGRP hits
load(file= 'Data/DGRP/Obs_LODs/All_Obs_LODs.rda')
#biallelic
dgrp.comb<-dgrp.comb[dgrp.comb$type=='b',]

effs<-c(0.05,0.1)
sss<-c(185)
dgrp.pos<-dgrp.comb[0,c('CHROM','POS','effect','samps')]

for(ee in effs)
{
  for(ss in sss)
  {
    
      dgrp.pos<-rbind(dgrp.pos, dgrp.comb[dgrp.comb$pow==1 & dgrp.comb$effect==ee & dgrp.comb$samps==ss,c('CHROM','POS','effect','samps')])
    
    
  }
}

pp<-rbind(mafs[['X']][,1:2],mafs[['2L']][,1:2],mafs[['2R']][,1:2],mafs[['3L']][,1:2],mafs[['3R']][,1:2])
colnames(pp)<-c('CHROM','POS')
pp.all<-merge(pp, dgrp.pos, by=c('CHROM','POS'))

(nrow(dgrp.pos)-nrow(pp.all))/nrow(dgrp.pos)

pp.dgrp<-merge(pp.all, imputed.genos, by.x=c('CHROM','POS'), by.y=c('chr','pos'))
gvar<-apply(pp.dgrp[,5:189], 1, var)

pp.dgrp$gvar<-gvar

pp.dspr<-pp.dgrp[,c('CHROM','POS','effect','samps','gvar')]

dspr.phenos<-t(apply(pp.dspr, 1, function(x) getPhenos_dgrp(chr=as.character(x[1]), pos=as.numeric(x[2]), eff=as.numeric(x[3]),gvar=as.numeric(x[5]),mafs=mafs, ids=ids, rils=pARILs)))


phenos<-data.frame('CHROM'=pp.dspr$CHROM, 'POS'=pp.dspr$POS,
                   'effect'=pp.dspr$effect, 'samps'=pp.dspr$samps,dspr.phenos,stringsAsFactors=FALSE)
colnames(phenos)[5:ncol(phenos)]<-colnames(dspr.phenos)
save(phenos, file="Data/DSPR/SimPhenos/DSPR_validate_btwpop.rda")



