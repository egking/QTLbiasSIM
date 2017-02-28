#############FUNCTIONS

#little function to get order for SNPs
decoder<-function(x,rils)
{
  which(rils==x)
}

getPhenos<-function(chr, pos, eff, mafs,ids,rils)
{
  mafs.pos<-mafs[[chr]][which(mafs[[chr]][,2]== pos),1:3]
  
  Aprobma<-do.call(rbind, strsplit(as.character(mafs.pos[,3]),"/"))
  
  Aprobma<-apply(Aprobma,2,function(x) as.numeric(x)/100)
  
  decodeA<-sapply(rils, function(x) decoder(x,ids[[chr]][,1]))
  
  Aprobma<-Aprobma[decodeA]
  
  genos<- rbinom(length(Aprobma),1,Aprobma)
  
  envs<-rnorm(length(genos),0,sqrt(((1/eff)-1)*var(genos)))
  
  phenos<-genos+envs
  
  names(phenos)<-rils
  
  return(phenos)
}

getPhenosMulti<-function(chr, pos, eff, mafs,ids,rils)
{
  gg<-matrix(NA,length(rils),length(pos))
  for(ii in 1:length(pos))
  {
  mafs.pos<-mafs[[chr]][which(mafs[[chr]][,2]== as.integer(pos[ii])),1:3]
  
  Aprobma<-do.call(rbind, strsplit(as.character(mafs.pos[,3]),"/"))
  
  Aprobma<-apply(Aprobma,2,function(x) as.numeric(x)/100)
  
  decodeA<-sapply(rils, function(x) decoder(x,ids[[chr]][,1]))
  
  Aprobma<-Aprobma[decodeA]
  
  gg[,ii]<- rbinom(length(Aprobma),1,Aprobma)
  }
  
  genos<-rowMeans(gg)
  
  envs<-rnorm(length(genos),0,sqrt(((1/eff)-1)*var(genos)))
  
  phenos<-genos+envs
  
  names(phenos)<-rils
  
  return(phenos)
}


#############FUNCTIONS
#############FUNCTIONS
#############FUNCTIONS
source("Functions/LOD_P_F_convert.R")






#make new phenotypes for all
set.seed(422001456)


load(file= 'Data/DSPR/Obs_LODs/All_Obs_LODs.rda')

#biallelic
dspr.comb<-dspr.comb[dspr.comb$type=='b',]


effs<-c(0.05,0.1,0.2)
sss<-c(100,185,500,878)

dspr.pos<-dspr.comb[0,c('CHROM','POS','effect','samps')]

for(ee in effs)
{
  for(ss in sss)
  {
    if(nrow(dspr.comb[dspr.comb$pow==1 & dspr.comb$effect==ee & dspr.comb$samps==ss,])<100)
    {
      
      dspr.pos<-rbind(dspr.pos, dspr.comb[dspr.comb$pow==1 & dspr.comb$effect==ee & dspr.comb$samps==ss,c('CHROM','POS','effect','samps')])
      
    }else{
      
      dd<-dspr.comb[dspr.comb$pow==1 & dspr.comb$effect==ee & dspr.comb$samps==ss,c('CHROM','POS','effect','samps')]
      dspr.pos<-rbind(dspr.pos, dd[sample(seq(1,nrow(dd)),100),])
      
    }
    
    
  }
}

#make phenos

#DSPR

library(DSPRqtl)
data(positionlist_wgenetic)

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
  
dspr.phenos<-t(apply(dspr.pos, 1, function(x) getPhenos(chr=as.character(x[1]), pos=as.numeric(x[2]), eff=as.numeric(x[3]),mafs=mafs, ids=ids, rils=pARILs)))
phenos<-data.frame('CHROM'=dspr.pos$CHROM, 'POS'=dspr.pos$POS,
                   'effect'=dspr.pos$effect, 'samps'=dspr.pos$samps,dspr.phenos,stringsAsFactors=FALSE)
colnames(phenos)[5:ncol(phenos)]<-colnames(dspr.phenos)
save(phenos, file="Data/DSPR/SimPhenos/DSPR_validate_biallelic.rda")


###########DSPR MULTI

load(file="Data/DSPR/SimPhenos/dspr.multi.pos.rda")

load(file= 'Data/DSPR/Obs_LODs/All_Obs_LODs.rda')

dspr.comb<-dspr.comb[dspr.comb$type=='m',]
dspr.comb$POS1<-dspr.multi.pos$POS1
dspr.comb$POS2<-dspr.multi.pos$POS2
dspr.comb$POS3<-dspr.multi.pos$POS3


effs<-c(0.05,0.1,0.2)
sss<-c(100,185,500,878)

dspr.pos<-dspr.comb[0,c('CHROM','POS1','POS2','POS3','effect','samps')]

for(ee in effs)
{
  for(ss in sss)
  {
    if(nrow(dspr.comb[dspr.comb$pow==1 & dspr.comb$effect==ee & dspr.comb$samps==ss,])<100)
    {
      
      dspr.pos<-rbind(dspr.pos, dspr.comb[dspr.comb$pow==1 & dspr.comb$effect==ee & dspr.comb$samps==ss,c('CHROM','POS1','POS2','POS3','effect','samps')])
      
    }else{
      
      dd<-dspr.comb[dspr.comb$pow==1 & dspr.comb$effect==ee & dspr.comb$samps==ss,c('CHROM','POS1','POS2','POS3','effect','samps')]
      dspr.pos<-rbind(dspr.pos, dd[sample(seq(1,nrow(dd)),100),])
      
    }
    
    
  }
}
dspr.phenos<-t(apply(dspr.pos, 1, function(x) getPhenosMulti(chr=as.character(x[1]), pos=as.numeric(x[2:4]), eff=as.numeric(x[5]),mafs=mafs, ids=ids, rils=pARILs)))
phenos<-data.frame('CHROM'=dspr.pos$CHROM, 'POS1'=dspr.pos$POS1,'POS2'=dspr.pos$POS2,'POS3'=dspr.pos$POS3,
                   'effect'=dspr.pos$effect, 'samps'=dspr.pos$samps,dspr.phenos,stringsAsFactors=FALSE)
colnames(phenos)[7:ncol(phenos)]<-colnames(dspr.phenos)
save(phenos, file="Data/DSPR/SimPhenos/DSPR_validate_multi.rda")

###DGRP

load(file= 'Data/DGRP/Obs_LODs/All_Obs_LODs.rda')
#biallelic
dgrp.comb<-dgrp.comb[dgrp.comb$type=='b',]

effs<-c(0.05,0.1,0.2)
sss<-c(100,185)
dgrp.pos<-dgrp.comb[0,c('CHROM','POS','effect','samps')]

for(ee in effs)
{
  for(ss in sss)
  {
    if(nrow(dgrp.comb[dgrp.comb$pow==1 & dgrp.comb$effect==ee & dgrp.comb$samps==ss,])<100)
    {
      
      dgrp.pos<-rbind(dgrp.pos, dgrp.comb[dgrp.comb$pow==1 & dgrp.comb$effect==ee & dgrp.comb$samps==ss,c('CHROM','POS','effect','samps')])
      
    }else{
      
      dd<-dgrp.comb[dgrp.comb$pow==1 & dgrp.comb$effect==ee & dgrp.comb$samps==ss,c('CHROM','POS','effect','samps')]
      dgrp.pos<-rbind(dgrp.pos, dd[sample(seq(1,nrow(dd)),100),])
      
    }
    
    
  }
}



load(file="Data/DGRP/imputed_genos.rda")

ss<-merge(imputed.genos, dgrp.pos, by.x=c('chr','pos'),by.y=c('CHROM','POS'))

genos<-data.matrix(ss[,3:187])
envs<-matrix(,nrow(genos),ncol(genos))
for(i in 1:length(genos[,1]))
{
  envs[i,]<-rnorm(length(genos[i,]),0,sqrt(((1/ss[i,188])-1)*var(genos[i,])))
}
phenos<-genos+envs

phenos<-data.frame('CHROM'=ss$chr, 'POS'=ss$pos, 'effect'=ss$effect, 'samps'=ss$samps,phenos,stringsAsFactors=FALSE)

save(phenos, file="Data/DGRP/SimPhenos/DGRP_validate_biallelic.rda")

#multi
DGRP.th<-7.4
load(file="Data/DGRP/Obs_LODs/map_combined_LOCO_multi.rda")
out.all$R2max<-LOD_R2(out.all$maxLOD,out.all$samps)
out.all$R2Q<-LOD_R2(out.all$QLOD,out.all$samps)
out.all$P<-getP(out.all$maxLOD,out.all$samps,1,out.all$samps-2)
out.all$pow<-rep(0,nrow(out.all))
out.all$pow[out.all$P>=DGRP.th]<-1
out.all$type<-'m'
out.all$POS<-rowMeans(out.all[,c('POS1','POS2','POS3')])

effs<-c(0.05,0.1,0.2)
sss<-c(100,185)
dgrp.comb<-out.all
dgrp.pos<-dgrp.comb[0,c('CHROM','POS1','POS2','POS3','effect','samps')]

for(ee in effs)
{
  for(ss in sss)
  {
    if(nrow(dgrp.comb[dgrp.comb$pow==1 & dgrp.comb$effect==ee & dgrp.comb$samps==ss,])<100)
    {
      
      dgrp.pos<-rbind(dgrp.pos, dgrp.comb[dgrp.comb$pow==1 & dgrp.comb$effect==ee & dgrp.comb$samps==ss,c('CHROM','POS1','POS2','POS3','effect','samps')])
      
    }else{
      
      dd<-dgrp.comb[dgrp.comb$pow==1 & dgrp.comb$effect==ee & dgrp.comb$samps==ss,c('CHROM','POS1','POS2','POS3','effect','samps')]
      dgrp.pos<-rbind(dgrp.pos, dd[sample(seq(1,nrow(dd)),100),])
      
    }
    
    
  }
}

dgrp.pos$pid<-seq(1,nrow(dgrp.pos))

ss1<-merge(imputed.genos, dgrp.pos, by.x=c('chr','pos'),by.y=c('CHROM','POS1'))
ss2<-merge(imputed.genos, dgrp.pos, by.x=c('chr','pos'),by.y=c('CHROM','POS2'))
ss3<-merge(imputed.genos, dgrp.pos, by.x=c('chr','pos'),by.y=c('CHROM','POS3'))

ss1<-ss1[order(ss1$pid),]
ss2<-ss2[order(ss2$pid),]
ss3<-ss3[order(ss3$pid),]

genos<-as.matrix(ss1[,3:187]+ss2[,3:187]+ss3[,3:187])

envs<-matrix(,nrow(genos),ncol(genos))
for(i in 1:length(genos[,1]))
{
  envs[i,]<-rnorm(length(genos[i,]),0,sqrt(((1/ss1$effect[i])-1)*var(genos[i,])))
}
phenos<-genos+envs

phenos<-data.frame('CHROM'=dgrp.pos$CHROM, 'POS1'=dgrp.pos$POS1,'POS2'=dgrp.pos$POS2,'POS3'=dgrp.pos$POS3, 'effect'=dgrp.pos$effect, 'samps'=dgrp.pos$samps,phenos,stringsAsFactors=FALSE)

save(phenos, file="Data/DGRP/SimPhenos/DGRP_validate_multi.rda")


