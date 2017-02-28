setwd("/home/kingeg/Projects/BeavisProj/")
library(ggplot2)
library(cowplot)

library(dplyr)
library(magrittr)
source("Functions/LOD_P_F_convert.R")

effs<-c(0.05,0.1,0.2)
sss<-c(100,185,500,878)

#DGRP.th<-7.4

#DSPR.th<-4.45


load(file="Data/DGRP/Obs_LODs/map_combined_LOCO_biallelic.rda")
out.all$R2max<-LOD_R2(out.all$maxLOD,out.all$samps)
out.all$R2Q<-LOD_R2(out.all$QLOD,out.all$samps)
out.all$P<-getP(out.all$maxLOD,out.all$samps,1,out.all$samps-2)
out.all$pow<-rep(0,nrow(out.all))
out.all$pow[out.all$P>=DGRP.th]<-1
out.all$type<-'b'
dgrp.comb<-out.all
rm(out.all)

load(file="Data/DGRP/Obs_LODs/map_combined_LOCO_multi.rda")
out.all$R2max<-LOD_R2(out.all$maxLOD,out.all$samps)
out.all$R2Q<-LOD_R2(out.all$QLOD,out.all$samps)
out.all$P<-getP(out.all$maxLOD,out.all$samps,1,out.all$samps-2)
out.all$pow<-rep(0,nrow(out.all))
out.all$pow[out.all$P>=DGRP.th]<-1
out.all$type<-'m'
out.all$POS<-rowMeans(out.all[,c('POS1','POS2','POS3')])
out.all<-out.all[,colnames(dgrp.comb)]
dgrp.comb<-rbind(dgrp.comb, out.all)

save(dgrp.comb,file= 'Data/DGRP/Obs_LODs/All_Obs_LODs.rda')


load(file="Data/DSPR/Obs_LODs/map_combined_LOCOseed.rda")
dspr.all$P<-getP(dspr.all$maxLOD,dspr.all$samps,7,dspr.all$samps-8)
dspr.all$pow<-rep(0,nrow(dspr.all))
dspr.all$pow[dspr.all$P>=DSPR.th]<-1
dspr.all$type<-'b'
dspr.comb<-dspr.all

load(file="Data/DSPR/Obs_LODs/map_combined_multiLOCO.rda")
dspr.all$P<-getP(dspr.all$maxLOD,dspr.all$samps,7,dspr.all$samps-8)
dspr.all$pow<-rep(0,nrow(dspr.all))
dspr.all$pow[dspr.all$P>=DSPR.th]<-1
dspr.all$type<-'m'
dspr.comb<-rbind(dspr.comb,dspr.all)

save(dspr.comb,file= 'Data/DSPR/Obs_LODs/All_Obs_LODs.rda')

font_size <- 10
my_theme <- theme(
  axis.text = element_text(size = font_size),
  axis.title = element_text(size = font_size),
  legend.text = element_text(size = font_size),
  plot.title = element_text(size = font_size + 1))

p1 <- dspr.comb %>% 
  filter(effect == effs[1] & pow == 1) %>% 
  ggplot(aes(x=factor(samps), y=R2max, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0.01,0.6)) +
  geom_hline(yintercept = effs[1]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  my_theme

p2 <- dspr.comb %>% 
  filter(effect == effs[2] & pow == 1) %>% 
  ggplot(aes(x=factor(samps), y=R2max, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0.01,0.6)) +
  geom_hline(yintercept = effs[2]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  my_theme

p3 <- dspr.comb %>% 
  filter(effect == effs[3] & pow == 1) %>% 
  ggplot(aes(x=factor(samps), y=R2max, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0.01,0.6)) +
  geom_hline(yintercept = effs[3]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = font_size)) +
  my_theme


p4 <- dgrp.comb %>% 
    filter(effect == effs[1] & pow == 1) %>% 
    ggplot(aes(x=factor(samps), y=R2max, fill = factor(type))) +
    geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
    geom_boxplot(outlier.size = 0, alpha=0.6) + 
    ylim(c(0.01,0.6)) +
    geom_hline(yintercept = effs[1]) +
    xlab("Sample Size") + ylab("Percent Variance Explained") +
    scale_fill_discrete(name="QTL Type",
                        breaks=c("b", "m"),
                        labels=c("biallelic", "multiallelic")) +
    theme(axis.title.x = element_blank(),
          legend.position = "none") +
    my_theme

p5 <- dgrp.comb %>% 
  filter(effect == effs[2] & pow == 1) %>% 
  ggplot(aes(x=factor(samps), y=R2max, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0.01,0.6)) +
  geom_hline(yintercept = effs[2]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.y = element_blank(),
        legend.position = "none") +
  my_theme

p6 <- dgrp.comb %>% 
  filter(effect == effs[3] & pow == 1) %>% 
  ggplot(aes(x=factor(samps), y=R2max, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0.01,0.6)) +
  geom_hline(yintercept = effs[3]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = font_size)) +
  my_theme

plots <- plot_grid(p1, p2, p3, p4, p5, p6,
                   ncol = 3, nrow = 2,
                   labels = c("a.", "b.", "c.","d.", "e.", "f."),
                   label_size = 9,
                   hjust = 0,
                   rel_widths=c(1, 0.9, 1.4),
                   align = "h")
                   
ggsave(plots, file = "Plots/DSPR_DGRP_beavis.pdf",width=8,height=6)






############################ ALL EFFECTS
############################ ALL EFFECTS
############################ ALL EFFECTS
p1 <- dspr.comb %>% 
  filter(effect == effs[1]) %>% 
  ggplot(aes(x=factor(samps), y=R2Q, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0,0.5)) +
  geom_hline(yintercept = effs[1]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  my_theme

p2 <- dspr.comb %>% 
  filter(effect == effs[2]) %>% 
  ggplot(aes(x=factor(samps), y=R2Q, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0,0.5)) +
  geom_hline(yintercept = effs[2]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  my_theme

p3 <- dspr.comb %>% 
  filter(effect == effs[3]) %>% 
  ggplot(aes(x=factor(samps), y=R2Q, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0,0.5)) +
  geom_hline(yintercept = effs[3]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = font_size)) +
  my_theme


p4 <- dgrp.comb %>% 
  filter(effect == effs[1]) %>% 
  ggplot(aes(x=factor(samps), y=R2Q, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0,0.5)) +
  geom_hline(yintercept = effs[1]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        legend.position = "none") +
  my_theme

p5 <- dgrp.comb %>% 
  filter(effect == effs[2]) %>% 
  ggplot(aes(x=factor(samps), y=R2Q, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0,0.5)) +
  geom_hline(yintercept = effs[2]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.y = element_blank(),
        legend.position = "none") +
  my_theme

p6 <- dgrp.comb %>% 
  filter(effect == effs[3]) %>% 
  ggplot(aes(x=factor(samps), y=R2Q, fill = factor(type))) +
  geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
  geom_boxplot(outlier.size = 0, alpha=0.6) + 
  ylim(c(0,0.5)) +
  geom_hline(yintercept = effs[3]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  scale_fill_discrete(name="QTL Type",
                      breaks=c("b", "m"),
                      labels=c("biallelic", "multiallelic")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = font_size)) +
  my_theme

plots <- plot_grid(p1, p2, p3, p4, p5, p6,
                   ncol = 3, nrow = 2,
                   labels = c("a.", "b.", "c.","d.", "e.", "f."),
                   label_size = 9,
                   hjust = 0,
                   rel_widths=c(1, 0.9, 1.4),
                   align = "h")

ggsave(plots, file = "Plots/DGRP_ALLeff_AtQTL.pdf",width=8,height=6)




############################ ALL EFFECTS
############################ ALL EFFECTS
############################ ALL EFFECTS

dspr.1<-subset(dspr.comb, samps %in% c(100,185) & pow==1 & type=='b')
dspr.1$pop<-'DSPR'

dgrp.1<-subset(dgrp.comb,type=='b' & pow==1)
dgrp.1$pop<-'DGRP'
dgrp.1<-dgrp.1[,colnames(dspr.1)]

dda<-rbind(dspr.1,dgrp.1)

#dda<-subset(dda, samps==100 & effect==0.1)
#summary(aov(log(dda$R2max)~dda$pop))

summary(aov(log(dspr.comb$R2max)~dspr.comb$samps+dspr.comb$effect+dspr.comb$type))

pl.list<-vector(mode='list',length=length(effs))

for(pp in 1:length(effs))
{
  tt<-subset(dda, effect==effs[pp])
  
  pl.list[[pp]] <- ggplot(tt, aes(x=factor(samps), y=R2max, fill = factor(pop))) +
    
    geom_point(pch=21,position=position_jitterdodge(), alpha=0.3) +
    geom_boxplot(outlier.size = 0, alpha=0.6) + 
    ylim(c(0,max(dda$R2max))) +
    geom_hline(yintercept = effs[pp]) +
    xlab("Sample Size") + ylab("Percent Variance Explained") +
    annotate("text", x = 2, y = 0.5, label = effs[pp]) +
    scale_fill_discrete(name="Population")
}

#pdf(file="Plots/DSPR_DGRP_v1.pdf",width=12,height=3.25)
#multiplot(pl.list[[1]],pl.list[[2]],pl.list[[3]],cols=3)
#dev.off()
  

dgrp.comb$pop<-'DGRP'
dspr.comb$pop<-'DSPR'
dd.all<-rbind(dgrp.comb, dspr.comb)

dd.all.p<-subset(dd.all, pow==1)

pall<-ggplot(dd.all, aes(x=factor(effect), y=R2max)) +
  # geom_point(aes(color=factor(samps), pch=factor(pop)),
  #            stat = "identity",
  #            size = 1,
  #            alpha = 1/5) +
  stat_summary(aes(color=factor(samps), pch=factor(pop)), 
    fun.y = mean,
               position = position_dodge(width = 0.5),
               geom = "point",
               size = 1) +
  stat_summary(aes(color=factor(samps), pch=factor(pop)),
               fun.data = mean_se, geom = "errorbar", width = 0.1)
pall







###############
types<-c('b','m')
effs<-c(0.05,0.1,0.2)
sss<-c(100,185,500,878)

dspr.sum<-data.frame("type"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "eff"=numeric(length=(length(sss)*length(effs)*length(types))), 
                     "samps"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "pow"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "Meff"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "SEeff"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "N"=numeric(length=(length(sss)*length(effs)*length(types))))
counter<-1
for(tt in types)
{
for(j in effs)
{
  for(k in sss)
  {
    dd<-subset(dspr.comb, effect==j & samps==k & type==tt)
    dspr.sum$type[counter]<-tt
    dspr.sum$eff[counter]<-j
    dspr.sum$samps[counter]<-k
    dspr.sum$pow[counter]<-mean(dd$pow)
    dspr.sum$Meff[counter]<-mean(dd[dd$pow==1,'R2max'])
    dspr.sum$SEeff[counter]<-sd(dd[dd$pow==1,'R2max'])/length(dd[dd$pow==1,'R2max'])
    dspr.sum$N[counter]<-length(dd[dd$pow==1,'R2max'])
    counter<-counter+1  
  }
}
}


types<-c('b','m')
effs<-c(0.05,0.1,0.2)
sss<-c(100,185)

dgrp.sum<-data.frame("type"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "eff"=numeric(length=(length(sss)*length(effs)*length(types))), 
                     "samps"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "pow"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "Meff"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "SEeff"=numeric(length=(length(sss)*length(effs)*length(types))),
                     "N"=numeric(length=(length(sss)*length(effs)*length(types))))
counter<-1
for(tt in types)
{
  for(j in effs)
  {
    for(k in sss)
    {
      dd<-subset(dgrp.comb, effect==j & samps==k & type==tt)
      dgrp.sum$type[counter]<-tt
      dgrp.sum$eff[counter]<-j
      dgrp.sum$samps[counter]<-k
      dgrp.sum$pow[counter]<-mean(dd$pow)
      dgrp.sum$Meff[counter]<-mean(dd[dd$pow==1,'R2max'])
      dgrp.sum$SEeff[counter]<-sd(dd[dd$pow==1,'R2max'])/length(dd[dd$pow==1,'R2max'])
      dgrp.sum$N[counter]<-length(dd[dd$pow==1,'R2max'])
      counter<-counter+1  
    }
  }
}

table(dgrp.comb$effect,dgrp.comb$samps,dgrp.comb$type)
table(dspr.comb$effect,dspr.comb$samps,dspr.comb$type)


dgrp.sum$pop<-'DGRP'
dspr.sum$pop<-'DSPR'
dd.all<-rbind(dgrp.sum, dspr.sum)
dd.all<-subset(dd.all, type=='b')

pall<-ggplot(dd.all, aes(x=eff, y=Meff)) +
  #geom_hline(yintercept = c(0.05,0.1,0.2), col='grey70')+
  geom_abline(slope=1,intercept=0,col='grey80')+
  geom_point(aes(color=factor(samps), pch=factor(pop)),
              stat = "identity",
              size = 4,
              alpha=1/2)+
  ylim(c(0.03,0.35))+
  scale_color_discrete(name = "Sample Size") +
  scale_shape_discrete(name = "Population") +
  ylab('Estimated % Variance Explained')+
  xlab('True % Variance Explained')
  
#pdf("Plots/True_Est_Biallelic.pdf",width=6,height=4)
#pall
#dev.off()





#####################

Ns<-seq(10,1500, by=10)
LLdgrp<-LOD_R2(6.72, Ns)

#pdf("Plots/R2_samplesize.pdf", width=4, height=4)
#plot(Ns, LLdgrp, xlab="Sample Size", ylab = "Minimal % Variance Detectable", pch=16)
#abline(v=185)
#dev.off()

LOD_R2(6.72, 100)
LOD_R2(6.8, 100)



#####POWER



p1 <- ggplot(dgrp.comb, aes(x = effect, y = pow, fill=type)) +
  stat_summary(fun.y = mean,
               #position = position_dodge(width = 0.5),
               pch=21,
               geom = "point",
               size = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  theme(legend.position = "none") +
  xlab("Effect") +
  ylab("Power") +
  ylim(c(0,1))+
  theme(axis.text.x = element_blank())
p1

