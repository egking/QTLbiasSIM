
ListTo3dArray<-function(oldlist)
{

  newarray <- array(dim=c(dim(oldlist[[1]]), length(oldlist)))
  for(i in seq(along=oldlist)) {newarray[,,i] <- as.matrix(oldlist[[i]])}

  dimnames(newarray)[[1]]<-rownames(oldlist[[1]])
  dimnames(newarray)[[2]]<-colnames(oldlist[[1]])
  if(is.null(names(oldlist))==FALSE)
     {
    dimnames(newarray)[[3]]<-names(oldlist)
  }
  return(newarray)
}

#big.array<-ListTo3dArray(big.list1)
