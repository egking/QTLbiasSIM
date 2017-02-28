#little function to get order for SNPs
decoder<-function(x,rils)
{
  which(rils==x)
}



getPhenos<-function(chrom, mafs, ids,rils, All.loc,eff)
{
  
  exp.s<-subset(All.loc, chr==chrom)
  exp.s$pid<-seq(1,nrow(exp.s))

  
  mafs1<-merge(exp.s[,c('chr','p1','pid')], mafs,by.y='V2',by.x="p1")
  mafs1<-mafs1[order(mafs1$pid),]
  mafs2<-merge(exp.s[,c('chr','p2','pid')], mafs,by.y='V2',by.x="p2")
  mafs2<-mafs2[order(mafs2$pid),]
  mafs3<-merge(exp.s[,c('chr','p3','pid')], mafs,by.y='V2',by.x="p3")
  mafs3<-mafs3[order(mafs3$pid),]
  
  Aprobma<-do.call(rbind, strsplit(as.character(mafs1[,5 ]),"/"))
  Aprobma<-apply(Aprobma,2,function(x) as.numeric(x)/100)
  decodeA<-sapply(rils, function(x) decoder(x,ids[,1]))
  Aprobma<-Aprobma[,decodeA]
  genos1<-apply(Aprobma, 2,function(x) rbinom(length(x),1,x))
  
  Aprobma<-do.call(rbind, strsplit(as.character(mafs2[,5 ]),"/"))
  Aprobma<-apply(Aprobma,2,function(x) as.numeric(x)/100)
  decodeA<-sapply(rils, function(x) decoder(x,ids[,1]))
  Aprobma<-Aprobma[,decodeA]
  genos2<-apply(Aprobma, 2,function(x) rbinom(length(x),1,x))
  
  Aprobma<-do.call(rbind, strsplit(as.character(mafs3[,5 ]),"/"))
  Aprobma<-apply(Aprobma,2,function(x) as.numeric(x)/100)
  decodeA<-sapply(rils, function(x) decoder(x,ids[,1]))
  Aprobma<-Aprobma[,decodeA]
  genos3<-apply(Aprobma, 2,function(x) rbinom(length(x),1,x))
  
  genos.all<-genos1+genos2+genos3
  
  envs<-matrix(,nrow(genos.all),ncol(genos.all))
  
  
  for(i in 1:length(genos.all[,1]))
  {
    envs[i,]<-rnorm(length(genos.all[i,]),0,sqrt(((1/eff)-1)*var(genos.all[i,])))
  }
  
  phenos<-genos.all+envs
  phenos<-data.frame('CHROM'=exp.s$chr, 'POS1'=exp.s$p1,'POS2'=exp.s$p2,'POS3'=exp.s$p3,
                     'effect'=eff,phenos)
  colnames(phenos)[6:ncol(phenos)]<-rils
  
  return(phenos)
}

set.seed(351109)

library(DSPRqtl)
data(positionlist_wgenetic)

chromosomes<-c('X','2L','2R','3L','3R')

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

#number of positions
npos<-10000

gg<-read.table(file="Data/genes_Dmel_5.txt", sep="\t",stringsAsFactors=FALSE,header=TRUE)

min(rilA[rilA$CHROM=='X','POS'])
min(rilA[rilA$CHROM=='2L','POS'])
min(rilA[rilA$CHROM=='2R','POS'])
min(rilA[rilA$CHROM=='3L','POS'])
min(rilA[rilA$CHROM=='3R','POS'])

max(rilA[rilA$CHROM=='X','POS'])
max(rilA[rilA$CHROM=='2L','POS'])
max(rilA[rilA$CHROM=='2R','POS'])
max(rilA[rilA$CHROM=='3L','POS'])
max(rilA[rilA$CHROM=='3R','POS'])

#gg<-subset(gg, posS > 192400 & posE < 21109758)

chromosomes<-c('X','2L','2R','3L','3R')
all.chr<-character(length=0)
all.pos<-matrix(NA,0,3)

for(cc in chromosomes)
{
  rilAX<-as.numeric(rilA[rilA$CHROM==cc,'POS'])
  gg.x<-as.matrix(gg[gg$chr==cc,c('posS','posE')])
  test<-matrix(NA,nrow(gg.x),3)
  for(i in 1:nrow(gg.x))
  {
  set<-rilAX[rilAX>=(gg.x[i,1]-1000) & rilAX<=(gg.x[i,2]+1000)]
  if(length(set)>=3)
  {
  test[i,]<-sample(set,3)
  }
  }
  all.pos<-rbind(all.pos,test)
  all.chr<-c(all.chr,rep(cc, nrow(gg.x)))
  cat(cc, "\n")
}

rand.pos<-data.frame('chr'=all.chr, 'p1'=all.pos[,1],'p2'=all.pos[,2],'p3'=all.pos[,3])
rand.pos<-rand.pos[-which(is.na(all.pos[,1])),]



######################set of effects
eff.all<-c(0.01,seq(0.025,0.25,by=0.025),0.5)

for(zz in 1:length(eff.all))
{
eff<-eff.all[zz]
All.loc<-rand.pos[sample(seq(1, nrow(rand.pos)),npos),]

All.loc$effect<-eff

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

outf<-paste("Data/DSPR/SimPhenos/pheno_multi3_gene_",100*eff,".rda",sep="")
save(pheno,file=outf)
cat(eff, "\n")

}