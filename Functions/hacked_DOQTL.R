################################################################################
# A faster impelmentation of the QTLRel algorithm for additive covariates with
# a kinship matrix.  For interactive covariates, use the main QTLRel function.
# This is only about twice as fast as QTLRel.
# Code adapted from Riyan Cheng's QTLRel package on CRAN.
# NOTE: No NAs are allowed anywhere! Filter your data before calling this.
# Daniel Gatti
# Dan.Gatti@jax.org
# Aug, 2, 2013

#MODIFIED BY KING
#fast.qtlrel.hap = haplotype based mapping
#fast.qtlrel.snp = snp based mapping
#for snp, give single column of snp probabilities instead of haplotype probs
################################################################################
# Arguments: pheno: numeric vector with no NAs containing the phenotype.
#            probs: 3D numeric array of haplotype probabilities.
#            K: numeric matrix.
fast.qtlrel.hap = function(pheno, probs, K, addcovar, snps) {
  
  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) < 200)
  
  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3], snp = snps[,1])
  class(prdat) = c(class(prdat), "addEff")
  err.cov = NULL
  
  if(!missing(K)) {
    
    K = as.matrix(K)
    vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
                EE = diag(nrow(K)))
    vc = NULL
    if(missing(addcovar)) {
      vc = estVC(y = pheno, v = vTmp)
    } else {
      vc = estVC(y = pheno, x = addcovar, v = vTmp)
    } # else
    
    err.cov = matrix(0, nrow(K), ncol(K))
    for(j in which(vc$nnl)) {
      err.cov = err.cov + vTmp[[j]] * vc$par[names(vTmp)[j]]
    } # for(j)
    
    rm(vTmp)
    
    # Invert the covariance matrix.
    eW = eigen(err.cov, symmetric = TRUE)
    if (min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)) {
      stop("fast.qtlrel: W is not positive definite")
    } else {
      eW$values[eW$values <= 0] = Inf
    } # else
    err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
    rm(eW)
    
  } # if(!missing(K))
  
  # Add the intercept.
  if(missing(addcovar)) {
    addcovar = matrix(1, nrow = length(pheno), ncol = 1, dimnames =
                        list(rownames(pheno), "Intercept"))
  } else {
    addcovar = as.matrix(cbind(rep(1, length(pheno)), addcovar))
    colnames(addcovar)[1] = "Intercept"
  } # else
  
  # Remove A as the basis.
  probs = probs[,-1,]
  
  # Null model.
  ss.null = 0
  if(!is.null(err.cov)) {
    
    ytmp = err.cov %*% pheno
    xtmp = err.cov %*% addcovar
    qr.null = qr(xtmp)
    ss.null = sum(qr.resid(qr.null, ytmp)^2)
    
  } else {
    
    qr.null = qr(addcovar)
    ss.null = sum(qr.resid(qr.null, pheno)^2)
    
  } # else
  
  # Additive model for all SNPs.
  addx = cbind(addcovar, probs[,,1])
  ss = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), ncol(addx), dimnames = list(snps[,1],
                                                           colnames(addx)))
  rng = (ncol(addcovar)+1):ncol(addx)
  
  perc.var = 0
  lrs = 0
  lod = 0
  
  if(!is.null(err.cov)) {
    
    for(s in 1:nrow(snps)) {
      addx[,rng] = probs[,,s]
      xtmp = err.cov %*% addx
      qr.add = qr(xtmp)
      ss[s] = sum(qr.resid(qr.add, ytmp)^2)
      coef[s,] = qr.coef(qr.add, ytmp)
    } # for(s)
    
  } else {
    
    for(s in 1:nrow(snps)) {
      
      addx[,rng] = probs[,,s]
      qr.add = qr(addx)
      ss[s] = sum(qr.resid(qr.add, pheno)^2)
      coef[s,] = qr.coef(qr.add, pheno)
      
    } # for(s)
    
  } # else
  
  perc.var = 100 * (1.0 - (ss / ss.null))
  lrs = -length(pheno) * log(ss / ss.null)
  lod = lrs / (2 * log(10))
  
  # Get the p-value from the LRS.
  p = pchisq(q = -length(pheno) * log(ss / ss.null), 
             df = ncol(addx) - ncol(addcovar), lower.tail = FALSE)
  
  return(list(lod = cbind(snps[,1:4], perc.var  = perc.var, lrs = lrs, 
                          lod = lod, p = p, neg.log10.p = -log(p, 10)), coef = coef))
  
} # fast.qtlrel()


fast.qtlrel.hapAB = function(pheno, probs, K, addcovar, snps) {
  
  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) < 200)
  
  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3], snp = snps[,1])
  class(prdat) = c(class(prdat), "addEff")
  err.cov = NULL
  
  if(!missing(K)) {
    
    K = as.matrix(K)
    vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
                EE = diag(nrow(K)))
    vc = NULL
    if(missing(addcovar)) {
      vc = estVC(y = pheno, v = vTmp)
    } else {
      vc = estVC(y = pheno, x = addcovar, v = vTmp)
    } # else
    
    err.cov = matrix(0, nrow(K), ncol(K))
    for(j in which(vc$nnl)) {
      err.cov = err.cov + vTmp[[j]] * vc$par[names(vTmp)[j]]
    } # for(j)
    
    rm(vTmp)
    
    # Invert the covariance matrix.
    eW = eigen(err.cov, symmetric = TRUE)
    if (min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)) {
      stop("fast.qtlrel: W is not positive definite")
    } else {
      eW$values[eW$values <= 0] = Inf
    } # else
    err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
    rm(eW)
    
  } # if(!missing(K))
  
  # Add the intercept.
  if(missing(addcovar)) {
    addcovar = matrix(1, nrow = length(pheno), ncol = 1, dimnames =
                        list(rownames(pheno), "Intercept"))
  } else {
    addcovar = as.matrix(cbind(rep(1, length(pheno)), addcovar))
    colnames(addcovar)[1] = "Intercept"
  } # else
  
  # Remove A as the basis.
  probs = probs[,-c(1,9),]
  
  # Null model.
  ss.null = 0
  if(!is.null(err.cov)) {
    
    ytmp = err.cov %*% pheno
    xtmp = err.cov %*% addcovar
    qr.null = qr(xtmp)
    ss.null = sum(qr.resid(qr.null, ytmp)^2)
    
  } else {
    
    qr.null = qr(addcovar)
    ss.null = sum(qr.resid(qr.null, pheno)^2)
    
  } # else
  
  # Additive model for all SNPs.
  addx = cbind(addcovar, probs[,,1])
  ss = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), ncol(addx), dimnames = list(snps[,1],
                                                           colnames(addx)))
  rng = (ncol(addcovar)+1):ncol(addx)
  
  perc.var = 0
  lrs = 0
  lod = 0
  
  if(!is.null(err.cov)) {
    
    for(s in 1:nrow(snps)) {
      addx[,rng] = probs[,,s]
      xtmp = err.cov %*% addx
      qr.add = qr(xtmp)
      ss[s] = sum(qr.resid(qr.add, ytmp)^2)
      coef[s,] = qr.coef(qr.add, ytmp)
    } # for(s)
    
  } else {
    
    for(s in 1:nrow(snps)) {
      
      addx[,rng] = probs[,,s]
      qr.add = qr(addx)
      ss[s] = sum(qr.resid(qr.add, pheno)^2)
      coef[s,] = qr.coef(qr.add, pheno)
      
    } # for(s)
    
  } # else
  
  perc.var = 100 * (1.0 - (ss / ss.null))
  lrs = -length(pheno) * log(ss / ss.null)
  lod = lrs / (2 * log(10))
  
  # Get the p-value from the LRS.
  p = pchisq(q = -length(pheno) * log(ss / ss.null), 
             df = ncol(addx) - ncol(addcovar), lower.tail = FALSE)
  
  return(list(lod = cbind(snps[,1:4], perc.var  = perc.var, lrs = lrs, 
                          lod = lod, p = p, neg.log10.p = -log(p, 10)), coef = coef))
  
} # fast.qtlrel()




fast.qtlrel.snp = function(pheno, probs, K, addcovar, snps) {
  
  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) < 200)
  
  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3], snp = snps[,1])
  class(prdat) = c(class(prdat), "addEff")
  err.cov = NULL
  
  if(!missing(K)) {
    
    K = as.matrix(K)
    vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
                EE = diag(nrow(K)))
    vc = NULL
    if(missing(addcovar)) {
      vc = estVC(y = pheno, v = vTmp)
    } else {
      vc = estVC(y = pheno, x = addcovar, v = vTmp)
    } # else
    
    err.cov = matrix(0, nrow(K), ncol(K))
    for(j in which(vc$nnl)) {
      err.cov = err.cov + vTmp[[j]] * vc$par[names(vTmp)[j]]
    } # for(j)
    
    rm(vTmp)
    
    # Invert the covariance matrix.
    eW = eigen(err.cov, symmetric = TRUE)
    if (min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)) {
      stop("fast.qtlrel: W is not positive definite")
    } else {
      eW$values[eW$values <= 0] = Inf
    } # else
    err.cov = eW$vector %*% diag(eW$values^-0.5) %*% t(eW$vector)
    rm(eW)
    
  } # if(!missing(K))
  
  # Add the intercept.
  if(missing(addcovar)) {
    addcovar = matrix(1, nrow = length(pheno), ncol = 1, dimnames =
                        list(rownames(pheno), "Intercept"))
  } else {
    addcovar = as.matrix(cbind(rep(1, length(pheno)), addcovar))
    colnames(addcovar)[1] = "Intercept"
  } # else
  
  # Remove A as the basis. DO NOT DO THIS FOR SNP ANALYSIS (ONLY ONE COLUMN)
  #probs = probs[,-1,]
  
  # Null model.
  ss.null = 0
  if(!is.null(err.cov)) {
    
    ytmp = err.cov %*% pheno
    xtmp = err.cov %*% addcovar
    qr.null = qr(xtmp)
    ss.null = sum(qr.resid(qr.null, ytmp)^2)
    
  } else {
    
    qr.null = qr(addcovar)
    ss.null = sum(qr.resid(qr.null, pheno)^2)
    
  } # else
  
  # Additive model for all SNPs.
  addx = cbind(addcovar, probs[,,1])
  ss = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), ncol(addx), dimnames = list(snps[,1],
                                                           colnames(addx)))
  rng = (ncol(addcovar)+1):ncol(addx)
  
  perc.var = 0
  lrs = 0
  lod = 0
  
  if(!is.null(err.cov)) {
    
    for(s in 1:nrow(snps)) {
      addx[,rng] = probs[,,s]
      xtmp = err.cov %*% addx
      qr.add = qr(xtmp)
      ss[s] = sum(qr.resid(qr.add, ytmp)^2)
      coef[s,] = qr.coef(qr.add, ytmp)
    } # for(s)
    
  } else {
    
    for(s in 1:nrow(snps)) {
      
      addx[,rng] = probs[,,s]
      qr.add = qr(addx)
      ss[s] = sum(qr.resid(qr.add, pheno)^2)
      coef[s,] = qr.coef(qr.add, pheno)
      
    } # for(s)
    
  } # else
  
  perc.var = 100 * (1.0 - (ss / ss.null))
  lrs = -length(pheno) * log(ss / ss.null)
  lod = lrs / (2 * log(10))
  
  # Get the p-value from the LRS.
  p = pchisq(q = -length(pheno) * log(ss / ss.null), 
             df = ncol(addx) - ncol(addcovar), lower.tail = FALSE)
  
  return(list(lod = cbind(snps[,1:4], perc.var  = perc.var, lrs = lrs, 
                          lod = lod, p = p, neg.log10.p = -log(p, 10)), coef = coef))
  
} # fast.qtlrel()

################################################################################
# Use the methods in Matrix eQTL to perform fast eQTL mapping for expression
# data. Method is from:
# Shabalin AA., Matrix eQTL: ultra fast eQTL analysis via large matrix
# operations, Bioinformatics. 2012 May 15;28(10):1353-8. PMID: 22492648
# Andrey Shabalin collaborated and wrote the code.
# Daniel Gatti
# Dan.Gatti@jax.org
# May 22, 2012

#KING NOTES
#this function is used in assoc.map. gives identical results with snp.scan (see below) as fast.qtlrel.snp

################################################################################
# Helper function for use in merge.analysis(). This uses the matrix eqtl 
# algorithm on a matrix of SNPs that are coded as 0, 0.5 and 1, with the 
# possibility of values in between those.
################################################################################
matrixeqtl.snps = function(pheno, geno, K, addcovar) {
  pheno = as.matrix(t(pheno))
  geno  = as.matrix(t(geno))
  # Create an error covariance matrix.
  if(!missing(K)) {
    eig = eigen(K, symmetric = TRUE)
    if(any(eig$values <= 0)) {
      stop("The covariance matrix is not positive definite")
    } # if(any(eig$values <= 0))
    correctionMatrix = eig$vectors %*% diag(1.0 / sqrt(eig$values)) %*% 
      t(eig$vectors)
    rm(eig)
  } else {
    correctionMatrix = numeric()
  } # else
  # Add an intercept and rotate the covariates.
  cvrt = matrix(1, nrow = 1, ncol = ncol(pheno))
  if(!missing(addcovar)) {
    cvrt = rbind(cvrt, t(addcovar))
  } # if(!missing(addcovar))
  if(length(correctionMatrix) > 0) {
    cvrt = cvrt %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  q = qr(t(cvrt))
  cvrt = t(qr.Q(q))
  # Rotate and center the genes.
  if(length(correctionMatrix) > 0) {
    pheno = pheno %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  pheno = pheno - tcrossprod(pheno, cvrt) %*% cvrt
  div = sqrt(rowSums(pheno^2))
  div[div == 0] = 1
  pheno = pheno / div
  rm(div)
  # Rotate and center the SNPs. 
  if(length(correctionMatrix) > 0) {
    geno = geno %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  geno = geno - tcrossprod(geno, cvrt) %*% cvrt
  div = sqrt(rowSums(geno^2))
  drop = div < 5 * sqrt(ncol(geno) * .Machine$double.eps)
  div[drop] = 1
  geno = geno / div
  geno[drop,] = 0
  # Note: we return the R^2.
  return(tcrossprod(geno, pheno)^2)
} # matrixeqtl.snps()


###KING FUNCTION FROM ASSOC.MAP AND OTHERS FROM DOQTL
#depends on regress library
#pheno is single column of phenotypes
#geno is single matrix of genotypes (rils by position)
#names for both and row and col names of K must be ril ids and must match
snp.scan<-function(pheno, K, geno)
{
#stopifnot(all.equal(names(pheno),rownames(geno)))
#stopifnot(all.equal(names(pheno),rownames(K)))
#stopifnot(all.equal(rownames(K),rownames(geno)))

err.cov = diag(length(pheno))
mod = regress(pheno ~ 1, ~K, pos = c(TRUE, TRUE))
err.cov = mod$sigma[1] * K + mod$sigma[2] * diag(length(pheno))

mm2<-matrixeqtl.snps(pheno=pheno, geno=geno, K=err.cov)
lrs2 = -length(pheno) * log(1.0 - mm2)
lod<-(lrs2 / (2 * log(10)))
return(lod)
}




############################
#CORE ESTVC FUNCTIONS FROM QTLREL (NOT MODIFIED BY KING)
#used in fast.qtlrel

machineEps<- (.Machine$double.eps)^(2/4)
inf<- max(1e+38,sqrt(.Machine$double.xmax))

# extract info from specified variance components
fv <- function(vv){
  # vv: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
  nms<- c("AA","DD","HH","AD","MH","EE")
  if(any(!is.element(nms,names(vv)))){
    cat("Assume components are in the order AA,DD,HH,AD,MH,EE...\n")
    names(vv)[1:6]<- nms
  }
  vvTmp<- vv
  vvTmp$AA<- vvTmp$DD<- vvTmp$HH<- vvTmp$AD<- vvTmp$MH<- vvTmp$EE<- NULL
  vv<- list(AA=vv$AA,
            DD=vv$DD,
            HH=vv$HH,
            AD=vv$AD,
            MH=vv$MH,
            EE=vv$EE)
  vv<- c(vv,vvTmp)
  if(is.null(vv[[nms[2]]])){
    if(!is.null(vv[[nms[3]]])){
      vv[nms[3]]<- list(NULL)
      cat(nms[3], "is set to null because", nms[2], "is null.\n")
    }
    if(!is.null(vv[[nms[5]]])){
      vv[nms[5]]<- list(NULL)
      cat(nms[5], "is set to null because", nms[2], "is null.\n")
    }
  }
  if(is.null(vv[[nms[1]]])){
    if(!is.null(vv[[nms[4]]])){
      vv[nms[4]]<- list(NULL)
      cat(nms[4], "is set to null because", nms[1], "is null.\n")
    }
  }
  if(is.null(vv[[nms[3]]])){
    if(!is.null(vv[[nms[4]]])){
      vv[nms[4]]<- list(NULL)
      cat(nms[4], "is set to null because", nms[3], "is null.\n")
    }
  }
  nv<- length(vv)
  if(nv>0){
    nnl<- NULL
    for(i in 1:nv){
      nnl<- c(nnl,!is.null(vv[[i]]))
      if(!is.null(vv[[i]])){
        if(!all(is.finite(vv[[i]])))
          stop("Only finite numerical elements are allowed in variance matrices!")
      }
    }
    nn<- cumsum(nnl)
  }else stop("At least the environmental variance component should be included.\n")
  
  list(v=vv,nv=nv,nnl=nnl,nn=nn)
}

# one of the main functions, Nelder-Mead method
estVC <-
  function(y,
           x,
           v = vector("list",6),
           initpar,
           nit = 25,
           method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
           control = list(),
           hessian = FALSE)
  {
    UseMethod("estVC")
  }

estVC.default <-
  function(y,
           x,
           v = vector("list",6),
           initpar,
           nit = 25,
           method = c("Nelder-Mead", "BFGS", "CG", "SANN"),
           control = list(),
           hessian = FALSE)
  {
    if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.")
    if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
        stop("x: missing or infinite data points not allowed.")
    
    if(!is.null(dim(y))){
      if(length(dim(y))>2) stop("y: may be wrong.\n")
      if(dim(y)[2]>1)
        warning("y: only the fisrt column will be analyzed.")
      y<- y[,1]
    }
    if(!missing(x)){
      oTmp<- data.frame(y=y,x)
    }else oTmp<- data.frame(y=y)
    oTmp<- model.frame(y~.,oTmp)
    y<- model.response(oTmp)
    x<- model.matrix(y~.,oTmp)
    method<-  match.arg(method)
    
    estVC.4(y = y,
            x = x,
            v = v,
            initpar = initpar,
            nit = nit,
            method = method,
            control = control,
            hessian = hessian)
  }

estVC.1 <-
  function(y,
           x,
           v,
           initpar,
           nit,
           method,
           control,
           hessian)
  {
    # estimate all background genetic variance (bgv)
    # y: vector, response
    # x: covariates
    # v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
    # initpar: initial parameters, will be initilized automatically if missing
    # nit: number of iterations to call optim()
    # method: the method to be used
    # control: A list of control parameters
    # hessian: logical. should a numerically differentiated Hessian matrix be returned?
    fs<- control$fnscale
    if(!is.null(fs) && fs<0) fs<- -1 else fs<- 1
    
    ny<- length(y)
    nb<- ncol(x)
    ov<- fv(v)
    if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))
    for(i in 1:ov$nv){
      if(ov$nnl[i]){
        initpar[nb+ov$nn[i]]<- initpar[nb+ov$nn[i]]/mean(diag(ov$v[[i]]))
        if(i!=4) initpar[nb+ov$nn[i]]<- log(initpar[nb+ov$nn[i]])
      }
    }
    
    optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
        tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
        if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }
      
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(fs*inf)
      
      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*inf)
      
      fs*tmp
    }
    
    oo<- list(par=initpar)
    val1<- val2<- inf
    while(nit>0){
      oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method=method,control=control,hessian=FALSE)
      if(ov$nnl[4]){
        tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
        if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
      }
      
      nit<- nit-1
      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
        break
    }
    
    oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
    if(ov$nnl[4]){
      tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
      if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
    }
    for(i in 1:ov$nv){
      if(ov$nnl[i]){
        if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
      }
    }
    oo$value<- as.numeric(oo$value)
    
    oo$value<- -oo$value #-fs*oo$value
    attributes(oo$value)<- NULL
    names(oo$par)[1:nb]<- colnames(x)
    names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
    oo$y<- as.matrix(y)
    oo$x<- as.matrix(x)
    oo$v<- v
    oo$nv<- ov$nv
    oo$nnl<- ov$nnl
    oo$nn<- ov$nn
    class(oo)<- "bgv"
    oo
  }

# NOTES: as accurate as but about 2.5 times as fast as estVC.1
estVC.2 <-
  function(y,
           x,
           v,
           initpar,
           nit,
           method,
           control,
           hessian)
  {
    # estimate all background genetic variance (bgv)
    # y: vector, response
    # x: covariates
    # v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
    # initpar: initial parameters, will be initilized automatically if missing
    # nit: number of iterations to call optim()
    # method: the method to be used
    # control: A list of control parameters
    # hessian: logical. should a numerically differentiated Hessian matrix be returned?
    fs<- control$fnscale
    if(!is.null(fs) && fs<0) fs<- -1 else fs<- 1
    
    ny<- length(y)
    nb<- ncol(x)
    ov<- fv(v)
    if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))
    for(i in 1:ov$nv){
      if(ov$nnl[i]){
        initpar[nb+ov$nn[i]]<- initpar[nb+ov$nn[i]]/mean(diag(ov$v[[i]]))
        if(i!=4) initpar[nb+ov$nn[i]]<- log(initpar[nb+ov$nn[i]])
      }
    }
    
    optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
        tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
        if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }
      
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(fs*inf)
      
      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*inf)
      
      fs*tmp
    }
    optfct.b<- function(a=list(par=initpar,nb=nb,ny=ny,ov=ov)){
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(a$par[a$nb+a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*a$par[a$nb+a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(a$par[1:a$nb])
      
      dd<- eigen(S,symmetric=T)
      uu<- dd$vec
      dd<- abs(dd$val)
      dd[dd<machineEps]<- machineEps
      yy<- t(uu)%*%y
      yy<- as.matrix(yy)
      xx<- t(uu)%*%x
      xx<- as.matrix(xx)
      b<- lm.wfit(x=xx,y=yy,w=1/dd)$coef
      b[!is.finite(b)]<- 0
      b
    }
    optfct.v<- function(par,a=list(par=initpar,nb=nb,ny=ny,ov=ov)){
      b<- a$par[1:a$nb]
      if(a$ov$nnl[4]){
        tmp<- sqrt(exp(par[a$ov$nn[1]]+par[a$ov$nn[3]])/2)
        if(abs(par[a$ov$nn[4]])>tmp) par[a$ov$nn[4]]<- sign(par[a$ov$nn[4]])*tmp
      }
      
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(par[a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*par[a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(fs*inf)
      
      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*inf)
      
      fs*tmp
    }
    
    oo<- list(par=initpar)
    val1<- val2<- inf
    pparTmp<- initpar
    if(length(initpar) < nb+2){
      while(nit>0){
        pparTmp[1:nb]<- optfct.b(a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
        oo<- optimize(optfct.v,interval=c(-1000,1000),a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),maximum=ifelse(fs<0,TRUE,FALSE))
        pparTmp[-c(1:nb)]<- oo$objective
        oo<- list(par=pparTmp, value=oo$minimum)
        if(ov$nnl[4]){
          tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
          if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
        }
        pparTmp<- oo$par
        
        val1<- val2
        val2<- oo$value
        if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
          break
        nit<- nit-1
      }
    }else{
      while(nit>0){
        pparTmp[1:nb]<- optfct.b(a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
        oo<- optim(pparTmp[-c(1:nb)],optfct.v,gr=NULL,a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),
                   method=method,control=control,hessian=FALSE)
        pparTmp[-c(1:nb)]<- oo$par
        oo$par<- pparTmp
        if(ov$nnl[4]){
          tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
          if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
        }
        pparTmp<- oo$par
        
        val1<- val2
        val2<- oo$value
        if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
          break
        nit<- nit-1
      }
    }
    oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
    if(ov$nnl[4]){
      tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
      if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
    }
    for(i in 1:ov$nv){
      if(ov$nnl[i]){
        if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
      }
    }
    
    oo$value<- as.numeric(oo$value)
    oo$value<- -oo$value
    attributes(oo$value)<- NULL
    names(oo$par)[1:nb]<- colnames(x)
    names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
    oo$y<- as.matrix(y)
    oo$x<- as.matrix(x)
    oo$v<- ov$v
    oo$nv<- ov$nv
    oo$nnl<- ov$nnl
    oo$nn<- ov$nn
    class(oo)<- "bgv"
    oo
  }

# NOTES: as accurate as but about 5 times as fast as estVC.1
# as accurate as but about 2 times as fast as estVC.2; may Not stable!!!
estVC.3 <-
  function(y,
           x,
           v,
           initpar,
           nit,
           method,
           control,
           hessian)
  {
    # estimate all background genetic variance (bgv)
    # y: vector, response
    # x: desig matrix including overall mean !!!
    # v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
    # initpar: initial parameters, will be initilized automatically if missing
    # nit: number of iterations
    control$fnscale<- 1
    ny<- length(y)
    nb<- ncol(x)
    ov<- fv(v)
    if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))
    for(i in 1:ov$nv){
      if(ov$nnl[i]){
        initpar[nb+ov$nn[i]]<- initpar[nb+ov$nn[i]]/mean(diag(ov$v[[i]]))
        if(i!=4) initpar[nb+ov$nn[i]]<- log(initpar[nb+ov$nn[i]])
      }
    }
    
    optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
        tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
        if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }
      
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(inf)
      
      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)
      
      tmp
    }
    optfct.b<- function(a=list(par=par,nb=nb,ny=ny,ov=ov)){
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(a$par[a$nb+a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*a$par[a$nb+a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(a$par[1:a$nb])
      
      dd<- eigen(S,symmetric=T)
      uu<- dd$vec
      dd<- abs(dd$val)
      dd[dd<machineEps]<- machineEps
      yy<- t(uu)%*%y
      yy<- as.matrix(yy)
      xx<- t(uu)%*%x
      xx<- as.matrix(xx)
      b<- lm.wfit(x=xx,y=yy,w=1/dd)$coef
      b[!is.finite(b)]<- 0
      b
    }
    optfct.v<- function(par,a=list(par=par,nb=nb,ny=ny,ov=ov)){
      b<- a$par[1:a$nb]
      if(a$ov$nnl[4]){
        tmp<- sqrt(exp(par[a$ov$nn[1]]+par[a$ov$nn[3]])/2)
        if(abs(par[a$ov$nn[4]])>tmp) par[a$ov$nn[4]]<- sign(par[a$ov$nn[4]])*tmp
      }
      
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(par[a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*par[a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(inf)
      
      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)
      
      tmp
    }
    
    val1<- val2<- inf
    pparTmp<- initpar
    while(nit>0){
      pparTmp[1:nb]<- optfct.b(a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
      oo<- nlm(optfct.v,pparTmp[-c(1:nb)],a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
      pparTmp[-c(1:nb)]<- oo$estimate; oo$value<- oo$minimum
      oo$par<- pparTmp
      if(ov$nnl[4]){
        tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
        if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
      }
      pparTmp<- oo$par
      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
        break
      nit<- nit-1
    }
    oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),method="Nelder-Mead",control=control,hessian=hessian)
    if(ov$nnl[4]){
      tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
      if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
    }
    for(i in 1:ov$nv){
      if(ov$nnl[i]){
        if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
      }
    }
    
    oo$value<- as.numeric(oo$value)
    oo$value<- -oo$value
    attributes(oo$value)<- NULL
    names(oo$par)[1:nb]<- colnames(x)
    names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
    oo$y<- as.matrix(y)
    oo$x<- as.matrix(x)
    oo$v<- ov$v
    oo$nv<- ov$nv
    oo$nnl<- ov$nnl
    oo$nn<- ov$nn
    class(oo)<- "bgv"
    oo
  }

# the above ones give more optional but less sensible solutions
estVC.4 <-
  function(y,
           x,
           v,
           initpar,
           nit,
           method,
           control,
           hessian)
  {
    # estimate all background genetic variance (bgv)
    # y: vector, response
    # x: covariates
    # v: list of variance components -- list(AA,DD,HH,AD,MH,EE,...) (Abney 200 pp635)
    # initpar: initial parameters, will be initilized automatically if missing
    # nit: number of iterations to call optim()
    # method: the method to be used
    # control: A list of control parameters
    # hessian: logical. should a numerically differentiated Hessian matrix be returned?
    fs<- control$fnscale
    if(!is.null(fs) && fs<0) fs<- -1 else fs<- 1
    
    ny<- length(y)
    nb<- ncol(x)
    ov<- fv(v)
    if(missing(initpar)) initpar<- c(rep(mean(y),ncol(x)),rep(var(y),sum(ov$nnl)))
    for(i in 1:ov$nv){
      if(ov$nnl[i]){
        initpar[nb+ov$nn[i]]<- initpar[nb+ov$nn[i]]/mean(diag(ov$v[[i]]))
        if(i!=4) initpar[nb+ov$nn[i]]<- log(initpar[nb+ov$nn[i]])
      }
    }
    
    optfct<- function(par,a=list(nb=nb,ny=ny,ov=ov)){
      b<- par[1:a$nb]
      if(a$ov$nnl[4]){
        tmp<- sqrt(exp(par[a$nb+a$ov$nn[1]]+par[a$nb+a$ov$nn[3]])/2)
        if(abs(par[a$nb+a$ov$nn[4]])>tmp) par[a$nb+a$ov$nn[4]]<- sign(par[a$nb+a$ov$nn[4]])*tmp
      }
      
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(par[a$nb+a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*par[a$nb+a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(fs*inf)
      
      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(fs*inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(fs*inf)
      
      fs*tmp
    }
    
    optfct.b<- function(a=list(par=par,nb=nb,ny=ny,ov=ov)){
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(a$par[a$nb+a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*a$par[a$nb+a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(a$par[1:a$nb])
      
      dd<- eigen(S,symmetric=T)
      uu<- dd$vec
      dd<- abs(dd$val)
      dd[dd<machineEps]<- machineEps
      yy<- t(uu)%*%y
      yy<- as.matrix(yy)
      xx<- t(uu)%*%x
      xx<- as.matrix(xx)
      b<- lm.wfit(x=xx,y=yy,w=1/dd)$coef
      b[!is.finite(b)]<- 0
      b
    }
    optfct.v<- function(par,a=list(par=par,nb=nb,ny=ny,ov=ov)){
      b<- a$par[1:a$nb]
      if(a$ov$nnl[4]){
        tmp<- sqrt(exp(par[a$ov$nn[1]]+par[a$ov$nn[3]])/2)
        if(abs(par[a$ov$nn[4]])>tmp) par[a$ov$nn[4]]<- sign(par[a$ov$nn[4]])*tmp
      }
      
      S<- matrix(0,nrow=a$ny,ncol=a$ny)
      for(i in 1:a$ov$nv){
        if(a$ov$nnl[i]){
          if(i!=4){
            S<- S + a$ov$v[[i]]*exp(par[a$ov$nn[i]])
          }else S<- S + a$ov$v[[i]]*par[a$ov$nn[i]]
        }
      }
      if(!all(is.finite(S))) return(inf)
      
      u<- x%*%b
      tmp<- qr(S)
      if(tmp$rank<ncol(tmp$qr)) return(inf)
      ddtmp<- abs(diag(tmp$qr))
      tmp<- log(2*pi)*(a$ny/2) + sum(log(ddtmp))/2 + t(y-u)%*%solve(tmp,y-u,tol=1e-15)*1/2
      if(!is.finite(tmp)) return(inf)
      
      tmp
    }
    
    upper<- rep(25,length(ov$v))
    names(upper)<- names(ov$v)
    upper["AD"]<- Inf
    upper<- upper[ov$nnl]
    upper<- c(rep(Inf,nb),upper)
    val1<- val2<- inf
    pparTmp<- initpar
    while(nit>0){
      pparTmp[1:nb]<- optfct.b(a=list(par=pparTmp,nb=nb,ny=ny,ov=ov))
      oo<- nlminb(pparTmp[-c(1:nb)],optfct.v,a=list(par=pparTmp,nb=nb,ny=ny,ov=ov),upper=upper[-c(1:nb)])
      pparTmp[-c(1:nb)]<- oo$par; oo$value<- oo$objective
      oo$par<- pparTmp
      if(ov$nnl[4]){
        tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
        if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
      }
      #         for(i in 1:ov$nv){
      #            if(ov$nnl[i]){
      #               if(i!=4) oo$par[nb+ov$nn[i]]<- min(10,exp(oo$par[nb+ov$nn[i]]))
      #            }
      #         }
      pparTmp<- oo$par
      val1<- val2
      val2<- oo$value
      if(abs(val2-val1)<machineEps && abs(val2-val1)<abs(val1)*machineEps)
        break
      nit<- nit-1
    }
    
    oo<- nlminb(oo$par,optfct,gradient=NULL,a=list(nb=nb,ny=ny,ov=ov),upper=upper)
    oo$value<- oo$objective
    if(ov$nnl[4]){
      tmp<- sqrt(exp(oo$par[nb+ov$nn[1]]+oo$par[nb+ov$nn[3]])/2)
      if(abs(oo$par[nb+ov$nn[4]])>tmp) oo$par[nb+ov$nn[4]]<- sign(oo$par[nb+ov$nn[4]])*tmp
    }
    oo<- optim(oo$par,optfct,gr=NULL,a=list(nb=nb,ny=ny,ov=ov),
               method="Nelder-Mead",control=list(maxit=1),hessian=hessian)
    for(i in 1:ov$nv){
      if(ov$nnl[i]){
        if(i!=4) oo$par[nb+ov$nn[i]]<- exp(oo$par[nb+ov$nn[i]])
      }
    }
    
    oo$value<- -fs*oo$value
    attributes(oo$value)<- NULL
    names(oo$par)[1:nb]<- colnames(x)
    names(oo$par)[nb+1:sum(ov$nnl)]<- names(ov$v[ov$nnl])
    oo$y<- as.matrix(y)
    oo$x<- as.matrix(x)
    oo$v<- ov$v
    oo$nv<- ov$nv
    oo$nnl<- ov$nnl
    oo$nn<- ov$nn
    class(oo)<- "bgv"
    oo
  }





####scan.R from QTLREL


W.inv<- function(W, symmetric=TRUE,inverse=TRUE){
  eW <- eigen(W, symmetric=symmetric)
  d <- eW$values
  if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
    stop("'W' is not positive definite")
  else d[d<=0]<- ifelse(inverse, Inf, 0)
  A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
  A # t(A)%*%A = W^{-1}
}

# adapted from lm.gls in MASS
lmGls<- function (formula, data, A, ...) {
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$A <- NULL
  m[[1L]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  yy <- model.response(m)
  y<- A%*%yy
  xx <- model.matrix(Terms, m, contrasts)
  x<- A%*%xx; colnames(x)[1]<- "(Intercept)"
  dtf<- data.frame(y=y,x)
  fit<- lm(y~.-1, data=dtf, ...)
  
  fit
}

# generalized least squares test
scanOne.0 <-
  function(y,
           x,
           prdat,
           cov,
           intcovar = NULL,
           test = c("None","F","Chisq"))
  {
    # prdat$pr: n by ? by ? matrix, allele probabilities
    # vc: object from estVC or aicVC
    # test: "Chisq", "F" or "Cp"
    diag.cov<- diag(cov)
    if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
      if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
        weights<- NULL
      }else weights<- 1/diag.cov
    }else weights<- NA
    gcv<- W.inv(cov)
    test<- match.arg(test)
    
    nsnp<- dim(prdat$pr)[3]
    if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
    model.par<- vector("list",nsnp)
    names(model.par)<- prdat$snp
    Pa<- rep(Inf,nsnp)
    names(Pa)<- prdat$snp
    Va<- Pa
    Pf<- rep(Inf,nsnp)
    names(Pf)<- prdat$snp
    Vf<- Pf
#     if(is.null(intcovar)){
#       if(!missing(x)){
#         oTmp<- data.frame(y=y,x)
#       }else{
#         oTmp<- data.frame(y=y)
#       }
#       if( !is.null(weights[1]) && is.na(weights[1]) ){
#         g0<- lmGls(y~.,data=oTmp,A=gcv)
#       }else{
#         g0<- lm(y~.,data=oTmp,weights=weights)
#       }
#       if(test=="None"){
#         P0<- logLik(g0)
#         for(k in 1:nsnp){
#           if(!missing(x)){
#             #oTmp<- data.frame(y=y,x,prdat$pr[,-1,k])
#             oTmp<- data.frame(y=y,x,prdat$pr[,,k])
#           }else{
#             #oTmp<- data.frame(y=y,a=prdat$pr[,-1,k])
#             oTmp<- data.frame(y=y,a=prdat$pr[,,k])
#           }
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(y~.,data=oTmp,A=gcv)
#           }else{
#             g<- lm(y~.,data=oTmp,weights=weights)
#           }
#           model.par[[k]]<- g$coef
#           P[k]<- logLik(g)
#           V[k]<- sum(g$res^2)
#         }
#         P<- 2*(P-P0)
#       }else{
#         for(k in 1:nsnp){
#           if(!missing(x)){
#             #oTmp<- data.frame(y=y,x,prdat$pr[,-1,k])
#             oTmp<- data.frame(y=y,x,prdat$pr[,,k])
#           }else{
#             #oTmp<- data.frame(y=y,prdat$pr[,-1,k])
#             oTmp<- data.frame(y=y,prdat$pr[,,k])
#           }
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(y~.,data=oTmp,A=gcv)
#           }else{
#             g<- lm(y~.,data=oTmp,weights=weights)
#           }
#           model.par[[k]]<- g$coef
#           P[k]<- anova(g0,g,test=test)$P[2]
#           V[k]<- sum(g$res^2)
#         }
#       }
#       V<- sum(g0$res^2) - V
#       V<- V/sum(anova(g0)[,"Sum Sq"])
#     }else{
    
        oTmp0<- data.frame(y=y)
        oTmp<- data.frame(y=y,intcovar)
  
      if( !is.null(weights[1]) && is.na(weights[1]) ){
        g00<- lmGls(y~.,data=oTmp0,A=gcv)
        g0<- lmGls(y~.,data=oTmp,A=gcv)
      }else{
        g00<- lm(y~.,data=oTmp0,weights=weights)
        g0<- lm(y~.,data=oTmp,weights=weights)
      }
     
        P00<-logLik(g00)
        P0<- logLik(g0)
         for(k in 1:nsnp){

           # oTmp<- data.frame(y=y,intcovar,prdat$pr[,-1,k])
            oTmp<- data.frame(y=y,intcovar,prdat$pr[,,k])

          #nc<- ncol(oTmp); nq<- ncol(prdat$pr[,-1,k])-1
          nc<- ncol(oTmp); nq<- 1
          str<- paste(paste("(",paste(colnames(oTmp)[2],collapse="+"),")",sep=""),
                      paste("(",paste(colnames(oTmp)[3],collapse="+"),")",sep=""),
                      sep=":")
          str<- paste("y~.+",str,sep="")
          
          if( !is.null(weights[1]) && is.na(weights[1]) ){
            ga<- lmGls(y~.,data=oTmp,A=gcv)
            gf<- lmGls(formula(str),data=oTmp,A=gcv)
          }else{
            ga<- lm(y~.,data=oTmp,weights=weights)
            gf<- lm(formula(str),data=oTmp,weights=weights)
            
          }
          
#          model.par[[k]]<- g$coef
          Pa[k]<- logLik(ga)
          Va[k]<- sum(ga$res^2)
          Pf[k]<- logLik(gf)
          Vf[k]<- sum(gf$res^2)
        }

      #P<- 2*(P-P0)
      
     # V<- sum(g0$res^2) - V
    #  V<- V/sum(anova(g0)[,"Sum Sq"])
    #}
    list('L0'=P00, 'Lac'=P0, 'Lag'=Pa, 'Lf'=Pf)
   # list(snp=prdat$snp,
    #     chr=prdat$chr,
    #     dist=prdat$dist,
    #     p=P,
    #     v=V*100,
    #     parameters=model.par)
  }

# scanOne.1 <-
#   function(y,
#            x,
#            prdat,
#            cov,
#            intcovar = NULL,
#            test = c("None","F","Chisq"))
#   {
#     # prdat$pr: n by 3 by ? matrix, conditional probabilities
#     # vc: object from estVC or aicVC
#     # test: "Chisq", "F" or "Cp"
#     diag.cov<- diag(cov)
#     if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
#       if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
#         weights<- NULL
#       }else weights<- 1/diag.cov
#     }else weights<- NA
#     gcv<- W.inv(cov)
#     test<- match.arg(test)
#     
#     nsnp<- dim(prdat$pr)[3]
#     if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
#     model.par<- vector("list",nsnp)
#     names(model.par)<- prdat$snp
#     P<- rep(Inf,nsnp)
#     names(P)<- prdat$snp
#     V<- P
#     if(is.null(intcovar)){
#       if(!missing(x)){
#         oTmp<- data.frame(y=y,x)
#       }else{
#         oTmp<- data.frame(y=y)
#       }
#       if( !is.null(weights[1]) && is.na(weights[1]) ){
#         g0<- lmGls(y~.,data=oTmp,A=gcv)
#       }else{
#         g0<- lm(y~.,data=oTmp,weights=weights)
#       }
#       if(test=="None"){
#         P0<- logLik(g0)
#         for(k in 1:nsnp){
#           if(!missing(x)){
#             oTmp<- data.frame(y=y,x,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
#           }else{
#             oTmp<- data.frame(y=y,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
#           }
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(y~.,data=oTmp,A=gcv)
#           }else{
#             g<- lm(y~.,data=oTmp,weights=weights)
#           }
#           model.par[[k]]<- g$coef
#           P[k]<- logLik(g)
#           V[k]<- sum(g$res^2)
#         }
#         P<- 2*(P-P0)
#       }else{
#         for(k in 1:nsnp){
#           if(!missing(x)){
#             oTmp<- data.frame(y=y,x,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
#           }else{
#             oTmp<- data.frame(y=y,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
#           }
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(y~.,data=oTmp,A=gcv)
#           }else{
#             g<- lm(y~.,data=oTmp,weights=weights)
#           }
#           model.par[[k]]<- g$coef
#           P[k]<- anova(g0,g,test=test)$P[2]
#           V[k]<- sum(g$res^2)
#         }
#       }
#       V<- sum(g0$res^2) - V
#       V<- V/sum(anova(g0)[,"Sum Sq"])
#     }else{
#       if(!missing(x)){
#         oTmp<- data.frame(y=y,x,intcovar)
#       }else{
#         oTmp<- data.frame(y=y,intcovar)
#       }
#       if( !is.null(weights[1]) && is.na(weights[1]) ){
#         g0<- lmGls(y~.,data=oTmp,A=gcv)
#       }else{
#         g0<- lm(y~.,data=oTmp,weights=weights)
#       }
#       if(test=="None"){
#         P0<- logLik(g0)
#         for(k in 1:nsnp){
#           if(!missing(x)){
#             oTmp<- data.frame(y=y,x,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
#           }else{
#             oTmp<- data.frame(y=y,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
#           }
#           #nc<- ncol(oTmp); nq<- ncol(prdat$pr[,-1,k])-1
#           nc<- ncol(oTmp); nq<- ncol(prdat$pr[,,k])-1
#           str<- paste(paste("(",paste(colnames(oTmp)[nc-nq-(nint:1)],collapse="+"),")",sep=""),
#                       paste("(",paste(colnames(oTmp)[(nc-nq):nc],collapse="+"),")",sep=""),
#                       sep=":")
#           str<- paste("y~.+",str,sep="")
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(formula(str),data=oTmp,A=gcv)
#           }else{
#             g<- lm(formula(str),data=oTmp,weights=weights)
#           }
#           model.par[[k]]<- g$coef
#           P[k]<- logLik(g)
#           V[k]<- sum(g$res^2)
#         }
#         P<- 2*(P-P0)
#       }else{
#         for(k in 1:nsnp){
#           if(!missing(x)){
#             oTmp<- data.frame(y=y,x,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
#           }else{
#             oTmp<- data.frame(y=y,intcovar,a=prdat$pr[,1,k]-prdat$pr[,3,k],d=prdat$pr[,2,k])
#           }
#           #nc<- ncol(oTmp); nq<- ncol(prdat$pr[,-1,k])-1
#           nc<- ncol(oTmp); nq<- ncol(prdat$pr[,,k])-1
#           str<- paste(paste("(",paste(colnames(oTmp)[nc-nq-(nint:1)],collapse="+"),")",sep=""),
#                       paste("(",paste(colnames(oTmp)[(nc-nq):nc],collapse="+"),")",sep=""),
#                       sep=":")
#           str<- paste("y~.+",str,sep="")
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(formula(str),data=oTmp,A=gcv)
#           }else{
#             g<- lm(formula(str),data=oTmp,weights=weights)
#           }
#           model.par[[k]]<- g$coef
#           P[k]<- anova(g0,g,test=test)$P[2]
#           V[k]<- sum(g$res^2)
#         }
#       }
#       V<- sum(g0$res^2) - V
#       V<- V/sum(anova(g0)[,"Sum Sq"])
#     }
#     
#     list(snp=prdat$snp,
#          chr=prdat$chr,
#          dist=prdat$dist,
#          p=P,
#          v=V*100,
#          parameters=model.par)
#   }
# 
# scanOne.2 <-
#   function(y,
#            x,
#            gdat,
#            cov,
#            intcovar = NULL,
#            numGeno = FALSE,
#            test = c("None","F","Chisq"))
#   {
#     # gdat: n by ? matrix, marker data. Markers in columes!!!
#     # vc: object from estVC or aicVC
#     # intcover: covariates that interact with QTL
#     # test: "Chisq", "F" or "Cp"
#     diag.cov<- diag(cov)
#     if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
#       if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
#         weights<- NULL
#       }else weights<- 1/diag.cov
#     }else weights<- NA
#     gcv<- W.inv(cov)
#     test<- match.arg(test)
#     if(numGeno){
#       num.geno<- I
#     }else num.geno<- as.factor
#     
#     nsnp<- dim(gdat)[2]
#     if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
#     model.par<- vector("list",nsnp)
#     names(model.par)<- colnames(gdat)
#     P<- rep(Inf,nsnp)
#     names(P)<- colnames(gdat)
#     V<- P
#     if(is.null(intcovar)){
#       if(!missing(x)){
#         oTmp<- data.frame(y=y,x)
#       }else{
#         oTmp<- data.frame(y=y)
#       }
#       if( !is.null(weights[1]) && is.na(weights[1]) ){
#         g0<- lmGls(y~.,data=oTmp,A=gcv)
#       }else{
#         g0<- lm(y~.,data=oTmp,weights=weights)
#       }
#       if(test=="None"){
#         P0<- logLik(g0)
#         for(j in 1:nsnp){
#           if(!missing(x)){
#             oTmp<- data.frame(y=y,x,snp=num.geno(gdat[,j]))
#           }else{
#             oTmp<- data.frame(y=y,snp=num.geno(gdat[,j]))
#           }
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(y~.,data=oTmp,A=gcv)
#           }else{
#             g<- lm(y~.,data=oTmp,weights=weights)
#           }
#           model.par[[j]]<- g$coef
#           P[j]<- logLik(g)
#           V[j]<- sum(g$res^2)
#         }
#         P<- 2*(P - P0)
#       }else{
#         for(j in 1:nsnp){
#           if(!missing(x)){
#             oTmp<- data.frame(y=y,x,snp=num.geno(gdat[,j]))
#           }else{
#             oTmp<- data.frame(y=y,snp=num.geno(gdat[,j]))
#           }
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(y~.,data=oTmp,A=gcv)
#           }else{
#             g<- lm(y~.,data=oTmp,weights=weights)
#           }
#           model.par[[j]]<- g$coef
#           P[j]<- anova(g0,g,test=test)$P[2]
#           V[j]<- sum(g$res^2)
#         }
#       }
#       V<- sum(g0$res^2) - V
#       V<- V/sum(anova(g0)[,"Sum Sq"])
#     }else{
#       if(!missing(x)){
#         oTmp<- data.frame(y=y,x,intcovar)
#       }else{
#         oTmp<- data.frame(y=y,intcovar)
#       }
#       if( !is.null(weights[1]) && is.na(weights[1]) ){
#         g0<- lmGls(y~.,data=oTmp,A=gcv)
#       }else{
#         g0<- lm(y~.,data=oTmp,weights=weights)
#       }
#       if(test=="None"){
#         P0<- logLik(g0)
#         for(j in 1:nsnp){
#           if(!missing(x)){
#             oTmp<- data.frame(y=y,x,intcovar,snp=num.geno(gdat[,j]))
#           }else{
#             oTmp<- data.frame(y=y,intcovar,snp=num.geno(gdat[,j]))
#           }
#           
#           nc<- ncol(oTmp)
#           str<- paste(colnames(oTmp)[nc-(nint:1)],colnames(oTmp)[nc],collapse="+",sep=":")
#           str<- paste("y~.+",str,sep="")
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(formula(str),data=oTmp,A=gcv)
#           }else{
#             g<- lm(formula(str),data=oTmp,weights=weights)
#           }
#           model.par[[j]]<- g$coef
#           P[j]<- logLik(g)
#           V[j]<- sum(g$res^2)
#         }
#         P<- 2*(P-P0)
#       }else{
#         for(j in 1:nsnp){
#           if(!missing(x)){
#             oTmp<- data.frame(y=y,x,intcovar,snp=num.geno(gdat[,j]))
#           }else{
#             oTmp<- data.frame(y=y,intcovar,snp=num.geno(gdat[,j]))
#           }
#           
#           nc<- ncol(oTmp)
#           str<- paste(colnames(oTmp)[nc-(nint:1)],colnames(oTmp)[nc],collapse="+",sep=":")
#           str<- paste("y~.+",str,sep="")
#           
#           if( !is.null(weights[1]) && is.na(weights[1]) ){
#             g<- lmGls(formula(str),data=oTmp,A=gcv)
#           }else{
#             g<- lm(formula(str),data=oTmp,weights=weights)
#           }
#           model.par[[j]]<- g$coef
#           P[j]<- anova(g0,g,test=test)$P[2]
#           V[j]<- sum(g$res^2)
#         }
#       }
#       V<- sum(g0$res^2) - V
#       V<- V/sum(anova(g0)[,"Sum Sq"])
#     }
#     
#     list(p=P,
#          v=V*100,
#          parameters=model.par)
#   }

scanOne<- 
  function(y,
           x,
           gdat,
           prdat = NULL,
           vc = NULL,
           intcovar = NULL,
           numGeno = FALSE,
           test = c("None","F","Chisq"),
           minorGenoFreq = 0,
           rmv = TRUE)
    
  {
    if(!all(is.finite(y)))
      stop("y: non-numeric or infinite data points not allowed.")
    if(!missing(x))
      if(any(sapply(x,is.infinite) | sapply(x,is.na)))
        stop("x: missing or infinite data points not allowed.")
    UseMethod("scanOne")
  }

scanOne.default<- 
  function(y,
           x,
           gdat,
           prdat = NULL,
           vc = NULL,
           intcovar = NULL,
           numGeno = FALSE,
           test = c("None","F","Chisq"),
           minorGenoFreq = 0,
           rmv = TRUE)
  {
    if(!is.null(vc)){
      if(is.element("bgv",attr(vc,"class"))){
        nb<- length(vc$par) - sum(vc$nnl)
        nr<- nrow(vc$y)
        cov<- matrix(0,nrow=nr,ncol=nr)
        for(i in 1:vc$nv)
          if(vc$nnl[i]) cov<- cov + vc$v[[i]]*vc$par[nb+vc$nn[i]]
      }else{
        if(is.data.frame(vc)) vc<- as.matrix(vc)
        if(!is.matrix(vc)) stop("vc should be a matrix.")
        if(!is.numeric(vc)) stop("vc should be a numeric matrix.")
        cov<- vc
      }
    }else cov<- diag(nrow(as.matrix(y)))
    if(!is.null(prdat)){
      if(is.element("addEff",class(prdat))){
        pv<- scanOne.0(y=y,x=x,prdat=prdat,cov=cov,intcovar=intcovar,test=test)
      }else{
        pv<- scanOne.1(y=y,x=x,prdat=prdat,cov=cov,intcovar=intcovar,test=test)
      }
    }else{
      if(any(is.na(gdat)))
        stop("There are missing genotypes...")
      tb<- sort(union(as.matrix(gdat),NULL))
      tbf<- NULL
      for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
      if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
      tbf<- apply(tbf,2,min)
      idx<- (tbf < nrow(gdat)*minorGenoFreq)
      if(sum(idx)>0){
        if(rmv){
          gdat<- as.matrix(gdat)
          tb<- sort(union(as.matrix(gdat),NULL))
          tbf<- NULL
          for(ii in tb) tbf<- rbind(tbf,colSums(gdat==ii))
          if(sum(tbf)!=nrow(gdat)*ncol(gdat)) stop("Error occurred.\n")
          tbf<- apply(tbf,2,min)
          idx<- (tbf < nrow(gdat)*minorGenoFreq)
          gdat<- gdat[,!idx]
          rm(tb,tbf,ii,idx)
        }else{
          cat("minor genotype frequency is too small at one or more SNPs.\n")
          return(NULL)
        }
      }
      
      gdat<- as.data.frame(gdat)
      pv<- scanOne.2(y=y,x=x,gdat=gdat,cov=cov,intcovar=intcovar,numGeno=numGeno, test=test)
    }
    
    class(pv)<- c("scanOne",test)
    pv
  }

print.scanOne <-
  function(x,...)
  {
    tt<- x; class(tt)<- NULL
    tt$parameters<- NULL
    tt<- as.data.frame(tt)
    if(length(tt$p)>5){
      cat("Test statistic:\n")
      print(tt[1:5,])
      cat("... ...\n\n")
      
      cat("Coefficients:\n")
      print(x$par[1:5])
      cat("... ...\n\n")
    }else{
      cat("Test statistic:\n")
      print(tt)
      cat("\n")
      
      cat("Coefficients:\n")
      print(x$par)
    }
  }

# generalized least squares estimates
gls<- function(formula,data=NULL,vc=NULL){
  if(is.null(data)){
    xx<- model.matrix(formula)
  }else xx<- model.matrix(formula,data)
  yy<- model.response(model.frame(formula,data))
  
  nr<- nrow(xx)
  if(!is.null(vc)){
    if(is.element("bgv",attr(vc,"class"))){
      nb<- length(vc$par) - sum(vc$nnl)
      nr<- nrow(vc$y)
      cov<- matrix(0,nrow=nr,ncol=nr)
      for(i in 1:vc$nv)
        if(vc$nnl[i]) cov<- cov + vc$v[[i]]*vc$par[nb+vc$nn[i]]
    }else{
      if(is.data.frame(vc)) vc<- as.matrix(vc)
      if(!is.matrix(vc)) stop("vc should be a matrix.")
      if(!is.numeric(vc)) stop("vc should be a numeric matrix.")
      cov<- vc
    }
  }else cov<- diag(nrow(as.matrix(yy)))
  A<- W.inv(cov)
  
  x<- A%*%xx; colnames(x)[1]<- "(Intercept)"
  y<- A%*%yy
  dtf<- data.frame(y=y,x)
  mdl<- lm(y~.-1, data=dtf)
  mdl$data<- dtf
  
  #   print(logLik(mdl))
  summary(mdl)$coeff
}

### main qtlrel function

qtl.qtlrel = function(pheno, probs, K, addcovar, intcovar, snps) {
  pheno = as.matrix(pheno)
  if(!missing(K)) {
    K = as.matrix(K)
  } # if(!missing(K))
  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) > 200)
  
  # Return value.
  retval = NULL
  
  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3],
               snp = snps[,1])
  vTmp = list(AA = NULL, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
              EE = diag(length(pheno)))
  
  if(!missing(K)) {
    vTmp$AA = 2 * K
  } # if(!missing(K))
  
  # This tells QTLRel to fit the additive model.
  class(prdat) = c(class(prdat), "addEff") 
  vc = NULL
  #if(missing(addcovar)) {
    # No covariates.
   # vc = estVC(y = pheno, v = vTmp)
    #res = scanOne(y = pheno, prdat = prdat, vc = vc, test = "None", numGeno = TRUE)
  #} else {
    vc = estVC(y = pheno, x = intcovar, v = vTmp)
   # if(missing(intcovar)) {
      # Additive covariates only.
    #  res = scanOne(y = pheno, x = addcovar, prdat = prdat, vc = vc, 
     #               numGeno = TRUE, test = "None")
    #} else {
      # Additive and interactive covariates.
      res = scanOne(y = pheno,  prdat = prdat, vc = vc, 
                    intcovar = intcovar, numGeno = TRUE, test = "None")
    #} # else
  #} # else
  
  # Convert the model coefficients to a matrix.
 # coef = matrix(unlist(res$parameters), length(res$parameters),
  #              length(res$parameters[[1]]), dimnames = list(res$snp,
    #                                                         names(res$parameters[[1]])), byrow = TRUE)
  
  # Return the LRS, LOD, p-value and -log10(p-value).
  #p = pchisq(q = res$p, df = dim(probs)[[2]], lower.tail = FALSE)
  #return(list(lod = cbind(snps[,1:4], perc.var = res$v, lrs = res$p,
   #                       lod = res$p / (2 * log(10)), p = p,
    #                      neg.log10.p = -log(p, 10)), coef = coef))
  return(res)
} # qtl.qtlrel()


