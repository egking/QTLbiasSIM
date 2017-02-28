

st.inds<-c(20001,21001,22001,23001,24001,25001, 40001,41001,42001, 80001)
           
en.inds<-st.inds+999

s<-100

counter<-1
for(k in 1:length(st.inds))
{
  st.ind<-st.inds[k]
  en.ind<-en.inds[k]
fname2<-paste("Data/DGRP/Obs_LODs/DGRPsetLOCO_",s,"_",st.ind,"_",en.ind,".rda",sep="")
load(file=fname2)
set$samps<-s
if(counter==1)
{
  out.all<-set
}else{
out.all<-rbind(out.all,set)
}
counter<-counter+1
  
}

st.inds<-c(20001,21001,22001, 40001,80001)

en.inds<-st.inds+999

s<-185

for(k in 1:length(st.inds))
{
  st.ind<-st.inds[k]
  en.ind<-en.inds[k]
  fname2<-paste("Data/DGRP/Obs_LODs/DGRPsetLOCO_",s,"_",st.ind,"_",en.ind,".rda",sep="")
  load(file=fname2)
  set$samps<-s
  if(counter==1)
  {
    out.all<-set
  }else{
    out.all<-rbind(out.all,set)
  }
  counter<-counter+1
  
}


save(out.all, file="Data/DGRP/Obs_LODs/map_combined_LOCO_biallelic.rda")


st.inds<-c(20001,21001,22001,23001,24001,25001, 40001,41001,42001, 80001)
s<-100
en.inds<-st.inds+999

counter<-1
for(k in 1:length(st.inds))
{
  
    st.ind<-st.inds[k]
    en.ind<-en.inds[k]
    fname2<-paste("Data/DGRP/Obs_LODs/DGRPsetLOCOmulti_",s,"_",st.ind,"_",en.ind,".rda",sep="")
    load(file=fname2)
    set$samps<-s
    
    if(counter==1)
    {
      out.all<-set
    }else{
      out.all<-rbind(out.all,set)
    }
    counter<-counter+1
  
}

st.inds<-c(20001,21001,22001,23001,24001, 40001, 80001)
s<-185
en.inds<-st.inds+999

for(k in 1:length(st.inds))
{
  
  st.ind<-st.inds[k]
  en.ind<-en.inds[k]
  fname2<-paste("Data/DGRP/Obs_LODs/DGRPsetLOCOmulti_",s,"_",st.ind,"_",en.ind,".rda",sep="")
  load(file=fname2)
  set$samps<-s
  
  if(counter==1)
  {
    out.all<-set
  }else{
    out.all<-rbind(out.all,set)
  }
  counter<-counter+1
  
}

save(out.all, file="Data/DGRP/Obs_LODs/map_combined_LOCO_multi.rda")

