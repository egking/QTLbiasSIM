
set.seed(7444209)

npos<-10000

load(file="Data/DGRP/imputed_genos.rda")
ww1<-which(imputed.genos$chr=='2L' & imputed.genos$pos>400000 & imputed.genos$pos<14900000)
ww2<-which(imputed.genos$chr=='2R' & imputed.genos$pos>9000000 & imputed.genos$pos<18000000)
ww3<-which(imputed.genos$chr=='3R' & imputed.genos$pos>6000000 & imputed.genos$pos<27000000)
ww4<-which(imputed.genos$chr=='4')

imputed.genos<-imputed.genos[-c(ww1,ww2,ww3,ww4),]
#exclude regions near inversions????
#excludes 1090771
#left with 1276049 total snps

######GENES
gg<-read.table(file="Data/genes_Dmel_5.txt", sep="\t",stringsAsFactors=FALSE,header=TRUE)

chromosomes<-c('X','2L','2R','3L','3R')
all.chr<-character(length=0)
all.pos<-matrix(NA,0,3)
genos1<-matrix(NA, nrow=nrow(gg),ncol=185)
genos2<-matrix(NA, nrow=nrow(gg),ncol=185)
genos3<-matrix(NA, nrow=nrow(gg),ncol=185)
counter<-1
for(cc in chromosomes)
{
  ii<-data.matrix(imputed.genos[imputed.genos$chr==cc,2:ncol(imputed.genos)])
  posX<-as.numeric(imputed.genos[imputed.genos$chr==cc,'pos'])
  gg.x<-as.matrix(gg[gg$chr==cc,c('posS','posE')])
  test<-matrix(NA,nrow(gg.x),3)
  for(i in 1:nrow(gg.x))
  {
    set<-posX[posX>=(gg.x[i,1]-1000) & posX<=(gg.x[i,2]+1000)]
    id.p<-which(posX>=(gg.x[i,1]-1000) & posX<=(gg.x[i,2]+1000))
    if(length(set)>=3)
    {
      ind.s<-sample(seq(1:length(set)),3)
      test[i,]<-set[ind.s]
      genos1[counter,]<-ii[id.p[ind.s[1]],2:ncol(ii)]
      genos2[counter,]<-ii[id.p[ind.s[2]],2:ncol(ii)]
      genos3[counter,]<-ii[id.p[ind.s[3]],2:ncol(ii)]
    }
    counter<-counter+1
  }
  all.pos<-rbind(all.pos,test)
  all.chr<-c(all.chr,rep(cc, nrow(gg.x)))
  cat(cc, "\n")
}

rand.pos<-data.frame('chr'=all.chr, 'p1'=all.pos[,1],'p2'=all.pos[,2],'p3'=all.pos[,3],stringsAsFactors=FALSE)
rand.pos<-rand.pos[-which(is.na(all.pos[,1])),]
genos1<-genos1[-which(is.na(all.pos[,1])),]
genos2<-genos2[-which(is.na(all.pos[,1])),]
genos3<-genos3[-which(is.na(all.pos[,1])),]
colnames(genos1)<-colnames(ii)[2:ncol(ii)]
colnames(genos2)<-colnames(ii)[2:ncol(ii)]
colnames(genos3)<-colnames(ii)[2:ncol(ii)]
#CHECK
#genos1[10001,1:10]
#tester<-subset(imputed.genos, chr==rand.pos[10001,1] & pos==rand.pos[10001,2])
#tester[1,1:10]


######################set of effects
eff.all<-c(0.01,seq(0.025,0.25,by=0.025),0.5)

for(k in eff.all)
{
  randSamp<-sample(seq(1,nrow(rand.pos)),npos)
  sub.pos<-rand.pos[randSamp,]

  genos1s<-genos1[randSamp,]
  genos2s<-genos2[randSamp,]
  genos3s<-genos3[randSamp,]
  
  genos<-genos1s+genos2s+genos3s
  
  envs<-matrix(,nrow(genos),ncol(genos))
  for(i in 1:length(genos[,1]))
  {
    envs[i,]<-rnorm(length(genos[i,]),0,sqrt(((1/k)-1)*var(genos[i,])))
  }
  phenos<-genos+envs
  #mm<-apply(genos,1,var)
  #ww<-apply(phenos[4:ncol(phenos)],1,var)
  #hist(mm/ww)
  phenos<-data.frame('CHROM'=sub.pos$chr, 'POS1'=sub.pos$p1,'POS2'=sub.pos$p2,'POS3'=sub.pos$p3, 'effect'=k,phenos,stringsAsFactors=FALSE)
  if(k==0.01)
  {
    all.phenos<-phenos
  }else{
    all.phenos<-rbind(all.phenos, phenos)
  }
}


save(all.phenos, file= "Data/DGRP/SimPhenos/allphenos_noinv_multi.rda")




