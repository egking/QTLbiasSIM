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

library(plyr)
source("Functions/LOD_P_F_convert.R")
eff.all<-c(0.01,seq(0.025,0.25,by=0.025),0.5)

effs<-c(0.05,0.2)
sss<-c(200,878)

DSPR.th<-6.8


load(file="Data/DSPR/Obs_LODs/CompareK/map_combined_v1.rda")

dspr.all$pow<-rep(0,nrow(dspr.all))
dspr.all$pow[dspr.all$maxLOD>=DSPR.th]<-1
dspr.all<-subset(dspr.all, effect==0.05 | effect==0.2)
dspr.all<-subset(dspr.all, samps==200 | samps==878)
dspr.all<-dspr.all[c(1:1000,10001:11000,20001:21000,30001:31000),]
dspr.all$type<-'No K'

dat.comb<-dspr.all
rm(dspr.all)
load(file="Data/DSPR/Obs_LODs/CompareK/map_combined_K.rda")
dspr.all$pow<-rep(0,nrow(dspr.all))
dspr.all$pow[dspr.all$maxLOD>=DSPR.th]<-1
dspr.all$type<-'Global K'
dat.comb<-rbind(dat.comb,dspr.all)
rm(dspr.all)

load(file="Data/DSPR/Obs_LODs/CompareK/map_combined_KLOCO.rda")
dspr.all$pow<-rep(0,nrow(dspr.all))
dspr.all$pow[dspr.all$maxLOD>=DSPR.th]<-1
dspr.all$type<-'LOCO'
dat.comb<-rbind(dat.comb,dspr.all)
rm(dspr.all)


pl.list<-vector(mode='list',length=length(effs))

for(pp in 1:length(effs))
{
tt<-subset(dat.comb, effect==effs[pp] & pow==1)

pl.list[[pp]] <- ggplot(tt, aes(factor(samps), R2max)) +
  geom_boxplot(aes(fill = factor(type))) + 
  geom_hline(yintercept = effs[pp]) +
  xlab("Sample Size") + ylab("Percent Variance Explained") +
  annotate("text", x = 2, y = max(tt$R2max), label = effs[pp]) +
  scale_fill_discrete(name="Model")
}

pdf(file="Plots/DSPR_K_LOCO.pdf",width=8,height=4)
multiplot(pl.list[[1]],pl.list[[2]],cols=2)
dev.off()


pl.list<-vector(mode='list',length=length(effs))

for(pp in 1:length(effs))
{
  tt<-subset(dat.comb, effect==effs[pp])
  
  pl.list[[pp]] <- ggplot(tt, aes(factor(samps), R2Q)) +
    geom_boxplot(aes(fill = factor(type))) + 
    geom_hline(yintercept = effs[pp]) +
    xlab("Sample Size") + ylab("Percent Variance Explained") +
    annotate("text", x = 2, y = max(tt$R2Q), label = effs[pp]) +
    scale_fill_discrete(name="Model")
}

pdf(file="Plots/DSPR_K_LOCO_AllPOS.pdf",width=8,height=4)
multiplot(pl.list[[1]],pl.list[[2]],cols=2)
dev.off()

tts<-c("No K", "Global K", "LOCO")
ss<-c(200, 878)
effs<-c(0.05,0.2)

oo<-data.frame('eff'=numeric(length=12),'ss'=numeric(length=12),
               'type'=character(length=12), 'pow'=numeric(length=12),stringsAsFactors=FALSE)
cc<-1
for(ii in ss)
{
  for(jj in effs)
  {
    for(zz in tts)
    {
      dd<-subset(dat.comb, samps==ii & effect==jj & type==zz)
      oo[cc,]<-c(jj,ii,zz,mean(dd$pow))
      cc<-cc+1
    }
  }
  
}

oo$ss_f <- as.factor(oo$ss)
oo$eff_f <- as.factor(oo$eff)

ggplot(oo, aes(x = eff_f, y = pow, color = type, shape = ss_f,
               group = ss_f)) +
  geom_point(position = position_jitter(width = 0.1),
             size = 3)

ggplot(dat.comb, aes(x = effect, y = pow, color = type, 
                     group = as.factor(samps),
                     linetype = as.factor(samps))) +
  geom_point()




