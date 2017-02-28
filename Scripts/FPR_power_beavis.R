multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


library(ggplot2)
library(cowplot)
library(tidyverse)

source("Functions/LOD_P_F_convert.R")

effs<-c(0.05,0.1,0.2)
sss<-c(100,185,500,878)

DGRP.th<-seq(4.25,7.75,by=0.25)

DSPR.th<-seq(1,6, by=0.25)

load(file="Data/DGRP/Obs_LODs/map_combined_LOCO_biallelic.rda")
out.all$R2max<-LOD_R2(out.all$maxLOD,out.all$samps)
out.all$R2Q<-LOD_R2(out.all$QLOD,out.all$samps)
out.all$P<-getP(out.all$maxLOD,out.all$samps,1,out.all$samps-2)
dgrp.comb<-out.all
rm(out.all)


load(file="Data/DSPR/Obs_LODs/map_combined_LOCOseed.rda")
dspr.all$P<-getP(dspr.all$maxLOD,dspr.all$samps,7,dspr.all$samps-8)
dspr.comb<-dspr.all

load(file='Data/DGRP/FPR_rates_ALL.rda')    
load(file='Data/DSPR/FPR_rates_ALL.rda')

all.n<-all.n[,3:17]
all.n<-all.n[,3:17]
dgrp.mm.185<-colMeans(all.n[2001:4000,])
dgrp.mm.100<-colMeans(all.n[1:2000,])


dspr.mm.878<-colMeans(L.hit[2001:4000,])
dspr.mm.100<-colMeans(L.hit[1:2000,])

eff.s<-0.1
samp.s<-185

dgrp.s<-subset(dgrp.comb, samps==samp.s & effect==eff.s)
dspr.s<-subset(dspr.comb, samps==samp.s & effect==eff.s)

pp.dgrp<-data.frame('P'=DGRP.th, 'pow'=numeric(length(DGRP.th)),
                    'powsd'=numeric(length(DGRP.th)),
                    'beav'=numeric(length(DGRP.th)),
                    'beavsd'=numeric(length(DGRP.th)),
                    'fp1'=numeric(length(DGRP.th)),
                    'fp1sd'=numeric(length(DGRP.th)),
                    'fp2'=numeric(length(DGRP.th)),
                    'fp2sd'=numeric(length(DGRP.th)))

pp.dspr<-data.frame('P'=DSPR.th, 'pow'=numeric(length(DSPR.th)),
                    'powsd'=numeric(length(DSPR.th)),
                    'beav'=numeric(length(DSPR.th)),
                    'beavsd'=numeric(length(DSPR.th)),
                    'fp1'=numeric(length(DSPR.th)),
                    'fp1sd'=numeric(length(DSPR.th)),
                    'fp2'=numeric(length(DSPR.th)),
                    'fp2sd'=numeric(length(DSPR.th)))


for(tt in 1:length(DSPR.th))
{
  pow<-rep(0,nrow(dspr.s))
  pow[dspr.s$P>=DSPR.th[tt]]<-1
  pp.dspr[tt,'pow']<-mean(pow)
  pp.dspr[tt,'powsd']<-sd(pow)/sqrt(nrow(dspr.s))
  pp.dspr[tt,'beav']<-mean(dspr.s[pow==1,'R2max'])
  pp.dspr[tt,'beavsd']<-sd(dspr.s[pow==1,'R2max'])/sqrt(length(pow[pow==1]))
  pp.dspr[tt,'fp1']<-mean(L.hit[1:2000,tt])
  pp.dspr[tt,'fp1sd']<-sd(L.hit[1:2000,tt])/sqrt(2000)
  pp.dspr[tt,'fp2']<-mean(L.hit[2001:4000,tt])
  pp.dspr[tt,'fp2sd']<-sd(L.hit[2001:4000,tt])/sqrt(2000)
  
}


for(tt in 1:length(DGRP.th))
{
  pow<-rep(0,nrow(dgrp.s))
  pow[dgrp.s$P>=DGRP.th[tt]]<-1
  pp.dgrp[tt,'pow']<-mean(pow)
  pp.dgrp[tt,'powsd']<-sd(pow)/sqrt(nrow(dgrp.s))
  pp.dgrp[tt,'beav']<-mean(dgrp.s[pow==1,'R2max'])
  pp.dgrp[tt,'beavsd']<-sd(dgrp.s[pow==1,'R2max'])/sqrt(length(pow[pow==1]))
  pp.dgrp[tt,'fp1']<-mean(all.n[1:2000,tt])
  pp.dgrp[tt,'fp1sd']<-sd(all.n[1:2000,tt])/sqrt(2000)
  pp.dgrp[tt,'fp2']<-mean(all.n[2001:4000,tt])
  pp.dgrp[tt,'fp2sd']<-sd(all.n[2001:4000,tt])/sqrt(2000)
  
}

pp.dspr$population<-'DSPR'
pp.dgrp$population<-'DGRP'

pp.all<-rbind(pp.dspr, pp.dgrp)

p1<- ggplot() +
  geom_point(data = pp.dspr, aes(x = P, y = fp1), color = "pink") +
  geom_point(data = pp.dspr, aes(x = P, y = fp2), color = "red") +
  geom_point(data = pp.dgrp, aes(x = P, y = fp1), color = "gray70") +
  geom_point(data = pp.dgrp, aes(x = P, y = fp2), color = "black") +
  geom_point(aes(x = 4.45, y = -3), pch = 2, color = "red") +
  geom_point(aes(x = 7.4, y = -3), pch = 2, color = "black") +
  xlab(expression("-log"[10]*"(P-value)")) + 
  ylab("# of False Positives") +
  scale_x_continuous(limits = c(1, 8))
p1

p2<- ggplot() +
  geom_point(data = pp.dspr, aes(x = P, y = beav), color = "red") +
  geom_point(data = pp.dgrp, aes(x = P, y = beav), color = "black") +
  xlab(expression("-log"[10]*"(P-value)")) + 
  ylab("Estimated PVE") +
  scale_y_continuous(limits=c(0.1,0.23))+
  scale_x_continuous(limits = c(1, 8))
p2

p3<- ggplot() +
  geom_point(data = pp.dspr, aes(x = P, y = pow), color = "red") +
  geom_point(data = pp.dgrp, aes(x = P, y = pow), color = "black") +
  xlab(expression("-log"[10]*"(P-value)")) + 
  ylab("Power") +
  scale_x_continuous(limits = c(1, 8))
p3


plots <- plot_grid(p1, p2, p3, ncol = 1, nrow = 3,
                   labels = c("a.", "b.", "c."))
ggsave(plots, file = "Plots/FPR_pow_beav.pdf",width=3,height=8)










