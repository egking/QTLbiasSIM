
######LIBRARIRES AND CONNECTIONS
library(RMySQL)

con <- dbConnect(MySQL(), user="eking",password="markov",dbname = "HMM")


Agenos<-dbGetQuery(con, "SELECT DISTINCT(A.ril) FROM HMMgenosA2 A, stocks B WHERE B.gstatus=1 AND A.ril=B.ril")

big.list<-vector('list',nrow(Agenos))

for (z in 1:nrow(Agenos))
{

  ril<-Agenos[z,]
  query<-paste("SELECT * FROM HMMregGeneticA2 WHERE ril=\"",ril,"\"",sep="")
  ril.dat<-dbGetQuery(con, query)   
  big.list[[z]]<-as.matrix(ril.dat[,40:47])
  if(z==1){pp<-ril.dat[,2:3]}
  cat(z, "\n")
}

arms<-c('X','2L','2R','3L','3R')

for(arm in arms)
{


#PHYSICAL
output<-matrix(NA,nrow(Agenos),nrow(Agenos))
rownames(output)<-Agenos[,1]
colnames(output)<-Agenos[,1]

diag(output)<-rep(1,nrow(output))

for (i in 1:(nrow(Agenos)-1))
{
  for (j in (i+1):nrow(Agenos)) 
    {  
    
      #print(c(i,j))
      ril.dat1<-big.list[[i]][pp$chr!=arm,]
      ril.dat2<-big.list[[j]][pp$chr!=arm,]
      big<-ril.dat1*ril.dat2
      bvec<-rowSums(big)
      r<-mean(bvec)
      output[i,j]<-r
      output[j,i]<-r
    
  }
}
Phys.Out<-output

oname<-paste("/home/eking/GWAS_DSPR/Kinship/kinship_Genetic_",arm,".rda",sep="")
save(Phys.Out,file=oname)
cat(arm)
}

dbDisconnect(con)



arms<-c('X','2L','2R','3L','3R')
kinall<-vector(mode='list',length=5)
names(kinall)<-arms

for(arm in arms)
{
  oname<-paste("/home/kingeg/Projects/BeavisProj/DSPR/InputFiles/kinship_Genetic_",arm,".rda",sep="")
  load(file=oname)
  kinall[[arm]]<-Phys.Out
}  

save(kinall, file="/home/kingeg/Projects/BeavisProj/DSPR/InputFiles/Kin_LOCO_A.rda")














