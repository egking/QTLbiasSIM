library(rrBLUP)
load(file="Data/DGRP/dgrp2.rda")

dgrp.dat[1:10,1:10]
dgrp.dat[1:10,213:214]

maf<-dgrp.dat$refc/(dgrp.dat$refc+dgrp.dat$altc)
which.maf<-which(maf<=0.95 & maf>0.05)
dgrp.dat<-dgrp.dat[which.maf,]

dgrp.dat<-dgrp.dat[,10:214]

dgrp.dat<-data.matrix(dgrp.dat)

nas<-apply(dgrp.dat,1,function(x) length(which(is.na(x))))
nas<-nas/ncol(dgrp.dat)
nas.c<-which(nas>=0.2)

dgrp.dat<-dgrp.dat[-nas.c,]

dgrp.dat[dgrp.dat==2]<- 1
dgrp.dat[dgrp.dat==0]<- -1

kinmat<-A.mat(t(dgrp.dat))

kinmat_master<-kinmat

diag(kinmat)<-0

all.lines<-character(0)

for(i in 1:nrow(kinmat))
{
  tt<-which(kinmat[i,]>0.5)
  zz<-colnames(kinmat)[tt]
  all.lines<-c(all.lines,zz)
}
all.lines<-unique(all.lines)


ii<-numeric(0)
for(k in all.lines)
{
  ss<-which(kinmat[k,]>0.5)
  if(length(ss)>1)
  {cat(k, "\t", ss, "\n")}
  
}

min_lines<-which(rownames(kinmat) %in% c('line_306','line_303',
                                         'line_358','line_348',
                                         'line_392','line_395',
                                         'line_385','line_850'))

kinmat1<-kinmat[-min_lines,-min_lines]

all.lines<-character(0)

for(i in 1:nrow(kinmat1))
{
  tt<-which(kinmat1[i,]>0.5)
  zz<-colnames(kinmat1)[tt]
  all.lines<-c(all.lines,zz)
  cat(i, zz,"\n")
}
all.lines<-unique(all.lines)


dd<-matrix(,length(all.lines),2)
for(k in 1:length(all.lines))
{
  if(!all.lines[k] %in% dd[,2])
  {
  dd[k,1]<-all.lines[k]
  dd[k,2]<-rownames(kinmat1)[which(kinmat1[all.lines[k],]>0.5)]
  }
  
}

dd<-dd[-which(is.na(dd[,1])),]

elim.set<-character(0)
for(ii in 1:nrow(dd))
{
  elim.set<-c(elim.set, sample(dd[ii,],1))  
}

elim.set<-c(elim.set,'line_306','line_303',
                       'line_358','line_348',
                       'line_392','line_395',
                       'line_385','line_850')

elim.n<-which(rownames(kinmat_master) %in% elim.set)

save(elim.set, file="Data/DGRP/dropped_lines.rda")

