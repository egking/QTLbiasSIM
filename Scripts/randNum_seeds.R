set.seed(9)
randNum<-as.integer(runif(100000,1,6e6))
save(randNum,file= "Data/randNum.rda")
