
setwd("/home/kingeg/Projects/BeavisProj/")
library(DSPRqtl)
library(ggplot2)
library(cowplot)

font_size <- 10
my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size),
  legend.text = element_text(size = font_size),
  plot.title = element_text(size = font_size + 1))



data(positionlist_wgenetic)

#get SNPs in DGRP with freqs
load(file="Data/DGRP/imputed_genos.rda")

ff<-as.matrix(imputed.genos[,3:187])
ff[ff< -1]<- -1
ff[ff> 1]<- 1

ff<-(ff+1)/2

ff<-rowMeans(ff) #each line separate. different than read count qualification

dgrp.f<-imputed.genos[,1:2]
dgrp.f$freq<-ff

dgrp.f<-dgrp.f[-which(dgrp.f$freq>0.975 | dgrp.f$freq<0.025),]
#drop those too.

dgrp.f$freq[dgrp.f$freq>0.5]<-1-dgrp.f$freq[dgrp.f$freq>0.5]


#get SNPs in DSPR with freqs

rilA<-read.table(file= "Data/DSPR_Release4/snpfreqs_RILS_A_X.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)

chromosomes<-c('2L','2R','3L','3R')

for(k in chromosomes)
{
  
  fname<-paste(file= "Data/DSPR_Release4/snpfreqs_RILS_A_",k,".txt",sep="")  
  ss<-read.table(fname,sep="\t",stringsAsFactors=FALSE,header=TRUE)  
  rilA<-rbind(rilA,ss) 
  
}

colnames(rilA)<-c("CHROM","POS","freqARILs")

rilA<-rilA[-which(rilA$freqARILs>0.975 | rilA$freqARILs<0.025),]

rilA$freqARILs[rilA$freqARILs>0.5]<-1-rilA$freqARILs[rilA$freqARILs>0.5]



#do comparisons

sh.p<-merge(rilA, dgrp.f, by.x=c('CHROM','POS'), by.y=c('chr','pos'), all=TRUE)

dgrp.un<-sh.p[is.na(sh.p$freqARILs),]
dspr.un<-sh.p[is.na(sh.p$freq),]

sh.p.set<-merge(rilA, dgrp.f, by.x=c('CHROM','POS'), by.y=c('chr','pos'), all=FALSE)

sh.p.set$dd<-sh.p.set$freqARILs-sh.p.set$freq


mean(dgrp.un$freq)
mean(dspr.un$freqARILs)


dgrp.p<-sh.p.set
dgrp.p$SNPs<-'DGRP shared'
dgrp.p<-dgrp.p[,c('CHROM','POS','freq','SNPs')]

dg1<-dgrp.un
dg1$SNPs<-'DGRP unique'
dg1<-dg1[,c('CHROM','POS','freq','SNPs')]

ds1<-dspr.un
ds1$SNPs<-'DSPR unique'
ds1<-ds1[,c('CHROM','POS','freqARILs','SNPs')]
colnames(ds1)<-c('CHROM','POS','freq','SNPs')

dspr.p<-sh.p.set
dspr.p$SNPs<-'DSPR shared'
dspr.p<-dspr.p[,c('CHROM','POS','freqARILs','SNPs')]
colnames(dspr.p)<-c('CHROM','POS','freq','SNPs')

all.p<-rbind(dgrp.p, ds1, dg1, dspr.p)


g1<-ggplot(all.p, aes(freq, fill=SNPs, color=SNPs)) +
  geom_density(alpha=0.1) +
  xlab('minor allele frequency') +
  my_theme

#scale_fill_manual(values=c('a'='red','b'='royalblue','c'='indianred', 'd'='cyan4'))

cor.test(sh.p.set$freqARILs, sh.p.set$freq)

nrow(all.p)
nrow(all.p[all.p$SNPs=='DGRP shared',])
nrow(all.p[all.p$SNPs=='DGRP unique',])
nrow(all.p[all.p$SNPs=='DSPR shared',])
nrow(all.p[all.p$SNPs=='DSPR unique',])

g2<-ggplot(sh.p.set, aes(dd)) +
  geom_histogram(bins=30, fill='gray80') +
  geom_vline(xintercept=0)+
  xlab('difference in MAF') +
  annotate('text', x=-0.3, y = 50000, label='DSPR\nMAF lower') +
  annotate('text', x=0.3, y = 50000, label='DGRP\nMAF lower') +
  my_theme

plots <- plot_grid(g1, g2,
                   ncol = 2, nrow = 1,
                   labels = c("a.", "b."),
                   label_size = 9)
                   
ggsave(plots, file="Plots/Freq_compare.pdf", width=8, height=3)

