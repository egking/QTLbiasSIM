
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

######################set of effects
eff.all<-c(0.01,seq(0.025,0.25,by=0.025),0.5)

for(k in eff.all)
{
  sub.pos<-imputed.genos[sample(seq(1,nrow(imputed.genos)),npos),]
  genos<-data.matrix(sub.pos[3:ncol(sub.pos)])
  envs<-matrix(,nrow(genos),ncol(genos))
  for(i in 1:length(genos[,1]))
  {
    envs[i,]<-rnorm(length(genos[i,]),0,sqrt(((1/k)-1)*var(genos[i,])))
  }
  phenos<-genos+envs
  #mm<-apply(genos,1,var)
  #ww<-apply(phenos[4:ncol(phenos)],1,var)
  #hist(mm/ww)
  phenos<-data.frame('CHROM'=sub.pos$chr, 'POS'=sub.pos$pos, 'effect'=k,phenos,stringsAsFactors=FALSE)
  if(k==0.01)
  {
    all.phenos<-phenos
  }else{
    all.phenos<-rbind(all.phenos, phenos)
  }
}


save(all.phenos, file= "Data/DGRP/SimPhenos/allphenos_noinv.rda")



