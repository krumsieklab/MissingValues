# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# R functions to generate "metabolomics" data and to introduce missingness
# according to different mechanisms
#
# Parts of this are experimental
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


### 1. logit function
logit <- function(x) log(x/(1-x))

### 2. logistic function
logistic <- function(x) exp(x)/(1+exp(x))

### 3. gen.complete.data: generate complete data (global mean=0,var=1)
gen.complete.data <- function(n,                          # sample size
                              p=1,                        # number of variables to be created
                              Cor=diag(1,p),              # Correlation among covariates
                              plus.pheno.last=T,          # should a phenotype be added as last column (as x2 without runday effects)
                              n.rundays=1,                # number of rundays, defaults to one: no runday effects
                              strength.runday.vary=rep(0,p)){ # strength with which runday means vary, i.e. ratio of variance of random intercepts (d2) versus variance of x (==1)
# current version lets runday levels vary in a similar way for all vars, i.e. a runday tends to have consistently higher/larger values for every met. 
    if(n.rundays==1) rundays <- rep(1,n) else rundays <- sample(cut(seq(1,n), breaks=n.rundays,labels=F))
    SD <- sqrt(strength.runday.vary)
    if(p==1) {
      means <- rnorm(n.rundays,sd=SD)
      means <- data.frame(means[rundays]) } else {
      means <- rmvnorm(n.rundays, sigma= diag(SD) %*% (diag(0.2,p)+0.8) %*% diag(SD)) # 0.8 experimental; so far only strength.runday.vary=rep(0,p) used
      means <- means[rundays,] }                                                                             
    if(p==1) x <- rnorm(n,sd=1) else x <- rmvnorm(n,sigma=Cor)
    if(plus.pheno.last) {
      means <- data.frame(means,0) 
      x <- cbind(x,x[,2]) }
    x <- x+means
  list(x=x,rundays=rundays)
  }
  
### 4. gen.miss.lod: generate missingness below a certain threshold (lod), i.e. truncate data at lod
gen.miss.lod <- function(datx,           # Data vector/matrix (function currently assumes the absence of ties)
                         incoms,         # (Vector of) percentage(s) of incomplete entries to be introduced into columns of datx; if last variable is pheno, incom=0
                         ret="x",...){   # what to return; x: data vector with missings; lod: LOD
  
  if(is.vector(datx)) {
    n <- length(datx)
    p <- 1 } else {
    n <- nrow(datx)
    p <- ncol(datx) }

  # # # # # # # # # # # # # # # # # # # # # 
  one.fun <- function(j) {
    incom <- incoms[j]
    if(p>1) x <- datx[,j]  else x <- datx
    if(incom>0) {      
      so <- sort(x)
      thres <- incom*n
      ganz <- floor(thres)
      rest <- thres-ganz 
      if(rest>0) newthres <- ganz + sample(c(0,1), 1, prob=c(1-rest,rest)) else newthres <- ganz
      lod <- (so[newthres+1]+so[newthres])/2    
      if(ret=="x") {
        x[x<lod] <- NA
        re <- x}
      if(ret=="lod") {
        re <- lod }   } else {
      if(ret=="x") re <- x
      if(ret=="lod") re <- NA}
    re}
  # # # # # # # # # # # # # # # # # # # # # 

  if(p==1)  res <- one.fun(1) else     res <- sapply(1:p,function(i) one.fun(i)) 
  res}
  

### 5. gen.miss.trend: generate missingness with increasing probability with decreasing concentration
gen.miss.trend <- function(datx,           # Data matrix/vector
                           incoms,       # (Vector of) percentage(s) of incomplete entries to be introduced intocolumns of datx
                           betas1,...) { # (Vector of) coefficient(s) determining strength of dependency of missingness from concentration (negativ: lower conc--> larger prob to be missing)

  if(is.vector(datx)) {
    n <- length(datx)
    p <- 1 } else {
    n <- nrow(datx)
    p <- ncol(datx) }
  whichmiss <- which(incoms!=0)
  pmiss <- length(whichmiss)

  # # # # # # # # # # # # # # # # # # # # # 
  one.fun <- function(j) {
    incom <- incoms[j]
    if(p>1) x <- datx[,j]  else x <- datx
    if(incom>0) {
      beta1 <- betas1[j]
      so <- sort(x)
      thres <- incom*n
      ganz <- floor(thres)
      rest <- thres-ganz 
      if(rest>0) newthres <- ganz + sample(c(0,1), 1, prob=c(1-rest,rest)) else newthres <- ganz
    
      f <- function(beta0) mean(logistic(beta0+beta1*x))-incom
      beta0 <- tryCatch(uniroot(f,c(-500,500),tol=1E-7)$root , error=function(e) NULL)
      if(is.null(beta0)) beta0 <- mean(logit(incom)-beta1*x)
      lp <- beta0 + beta1*x
      Prob  <- as.numeric(lp>0)-sign(lp)*logistic(-abs(lp))     
      # actually, correct choice of beta0 is not so important here, since in each multinomial experiment, exactly 1 of the n observations is drawn, where their respective probabilty should have only a relative influence on the probability to be drawn

      ids.to.miss <- which(rmultinom(1,1,prob=Prob)[,1]==1)
      for (i in 2:newthres) {
        which.in <- setdiff(1:n,ids.to.miss)
        ids.to.miss <- c(ids.to.miss, which.in[which(rmultinom(1,1,prob=Prob[which.in])[,1]==1)]) }  
      x[ids.to.miss] <- NA}
    x }
  # # # # # # # # # # # # # # # # # # # # # 

  if(is.vector(datx))  res <- one.fun(1) else   res <- sapply(1:p,function(i) one.fun(i))    
  res}
 

### 6.  gen.miss.trend.runday: generate missingness with increasing probability with decreasing concentration, with incom varying between rundays
gen.miss.trend.runday <- function(datx,                        # Data matrix/vector
                                  rundays=NULL,                # Runday assignment (if already available)
                                  n.rundays=1,                 # Number of rundays (if data are not already assigned to rundays) 
                                  incoms,                      # (Vector of) percentage(s) of incomplete entries to be introduced into x (overall)
                                  betas1,                      # (Vector of) coefficient(s) determining strength of dependency of missingness from concentration
                                  cors.incom=1,
                                  cors.beta1=cors.incom,
                                  strengths.incom.vary=NULL,   # If NULL; incom is not varied across the rundays; else, strength with which incom varies between rundays, i.e. variance of random intercepts (d2)
                                  strengths.beta1.vary=NULL,...){ # If NULL; beta1 is not varied across the rundays; else, strength with which incom varies between rundays, i.e. ratio of variance of random intercepts (d2) versus variance of x

  if(is.vector(datx)) {
    n <- length(datx)
    p <- 1 } else {
    n <- nrow(datx)
    p <- ncol(datx) }
  whichmiss <- which(incoms!=0)
  pmiss <- length(whichmiss)
    
  if(is.null(rundays)) {
    if(n.rundays==1) rundays = rep(1,n) else rundays <- sample(cut(seq(1,n), breaks=n.rundays,labels=F))  } else n.rundays=length(unique(rundays))

  newthress <- sapply(whichmiss,function(j) {
    if(p>1) x <- datx[,j]  else x <- datx
    so <- sort(x)
    thres <- incoms[j]*n
    ganz <- floor(thres)
    rest <- thres-ganz
    if(rest>0) newthres <- ganz + sample(c(0,1), 1, prob=c(1-rest,rest)) else newthres <- ganz
    newthres})

  if(!is.null(strengths.incom.vary)& any(!is.na(strengths.incom.vary))) {
    if(pmiss==1) {
      incommat <- data.frame(incoms[whichmiss]+rnorm(n.rundays,mean=0,sd=sqrt(strengths.incom.vary[1]))) } else {
      incommat <- rmvnorm(n.rundays,mean=incoms[whichmiss], sigma=diag(sqrt(strengths.incom.vary)) %*% (diag(1-cors.incom,pmiss)+cors.incom) %*% diag(sqrt(strengths.incom.vary))) }   
    incommat[incommat<0] <- 0
    incommat[incommat>1] <- 1} else incommat <-  matrix(rep(incoms[whichmiss],n.rundays),byrow=T,nrow=n.rundays)
    
  incommat1 <- sapply(1:ncol(incommat),function(j) incommat[,j]*incoms[j]/sum(incommat[,j])*nrow(incommat))
  incomabs <- floor(incommat1 *as.numeric(table(rundays)))
  incomabs <- sapply(1:ncol(incommat), function(j) {
    while(sum(incomabs[,j])>newthress[j]) incomabs[,j] <- incomabs[,j]-1
    incomabs[,j]})
  restmat <- as.matrix(incommat1*as.numeric(table(rundays))-incomabs )
  restmiss <- newthress - apply(incomabs,2,sum)
  for (j in 1:ncol(incommat)) {
    if(restmiss[j]>0) for (i in 1:restmiss[j]) {
      ids.to.miss <- which(rmultinom(1,1,prob=restmat[,j])[,1]==1) 
      incomabs[ids.to.miss,j] <- incomabs[ids.to.miss,j] +1 
      restmat[ids.to.miss,j] <- restmat[ids.to.miss,j] -1
      if(restmat[ids.to.miss,j] <0) restmat[ids.to.miss,j] <- 0
      if(sum(restmat[,j])==0) restmat[,j] <- incommat[,j] }
    }

  if(!is.null(strengths.beta1.vary) & any(!is.na(strengths.beta1.vary))) {
    if(pmiss==1) {
      beta1mat <- data.frame(betas1[whichmiss]+rnorm(n.rundays,mean=0,sd=sqrt(strengths.beta1.vary[1]))) } else { 
      beta1mat <- rmvnorm(n.rundays,mean=betas1[whichmiss], sigma=diag(sqrt(strengths.beta1.vary)) %*% (diag(1-cors.beta1,pmiss)+cors.beta1) %*% diag(sqrt(strengths.beta1.vary))) }   
    beta1mat[beta1mat>0] <- 0} else beta1mat <-  matrix(rep(betas1[whichmiss],n.rundays),byrow=T,nrow=n.rundays)

  for (j in 1:ncol(incommat)) {
    if(sum(incommat[,j])>0) {
      Probvec <- rep(NA,n)
      if(p>1) x <- datx[,whichmiss[j]] else x <- datx
      for (i in 1:n.rundays) {
        if(incommat[i,j]==0) { Prob_i <- rep(0,length(which(rundays==i))) } else {
          f <- function(beta0) mean(logistic(beta0+beta1mat[i,j]*x[rundays==i]))-incommat[i,j]
          beta0 <- tryCatch(uniroot(f,c(-500,500),tol=1E-7)$root , error=function(e) NULL)
          if(is.null(beta0)) beta0 <- mean(logit(incommat[i,j])-beta1mat[i,j]*x[rundays==i])
          lp <- beta0 + beta1mat[i,j]*x[rundays==i]
          Prob_i  <- as.numeric(lp>0)-sign(lp)*logistic(-abs(lp)) }   
        Probvec[rundays==i] <- Prob_i}  
        if(sum(Probvec)==0) Probvec <- Probvec+0.1 
        ids.to.miss <- which(rmultinom(1,1,prob=Probvec)[,1]==1)
        if(1 == incomabs[rundays[ids.to.miss],j]) ids.to.remove <- which(rundays==rundays[ids.to.miss]) else ids.to.remove <- ids.to.miss
        for (k in 2:newthress[j]) {
          which.in <- setdiff(1:n,c(ids.to.miss,ids.to.remove)) 
          if(sum(Probvec[which.in])==0) Probvec[which.in] <- Probvec[which.in]+0.1 
          new.id.to.miss <- which.in[which(rmultinom(1,1,prob=Probvec[which.in])[,1]==1)]
          ids.to.miss <- c(ids.to.miss, new.id.to.miss) 
          if(length(which(rundays[ids.to.miss]==rundays[new.id.to.miss])) == incomabs[rundays[new.id.to.miss],j]) ids.to.remove <- c(ids.to.remove,which(rundays==rundays[new.id.to.miss])) else ids.to.remove <- c(ids.to.remove,new.id.to.miss)}
        if(p>1) datx[ids.to.miss,whichmiss[j]] <- NA else datx[ids.to.miss] <- NA }  }
  datx}


### 7.  gen.miss.mcar: generate random missingness (only mechanisms implemented to be possible on top of the other mechanisms)
gen.miss.mcar <- function(datx,                        # Data vector/matrix
                          incoms,...) {                # Percentage of incomplete entries to be introduced into x (to be achieved overall, if x already contained missings)
                          
  if(is.vector(datx)) {
    n <- length(datx)
    p <- 1 } else {
    n <- nrow(datx)
    p <- ncol(datx) }

  # # # # # # # # # # # # # # # # # # # # # 
  one.fun <- function(j) {
    incom <- incoms[j]
    if(p>1) x <- datx[,j] else x <- datx
    if(incom>0) {      
      nmiss <- incom*n
      ganz <- floor(nmiss)
      rest <- nmiss-ganz 
      if(rest>0) nmiss <- ganz + sample(c(0,1), 1, prob=c(1-rest,rest)) else nmiss <- ganz
  
      nmiss0 <- length(which(is.na(x)))
      nmiss1 <- nmiss-nmiss0
    
      ids.to.miss <- sample(which(!is.na(x)), nmiss1)
      x[ids.to.miss] <- NA }
    x}
  # # # # # # # # # # # # # # # # # # # # # 

  if(p==1)  res <- one.fun(1) else     res <- sapply(1:p,function(i) one.fun(i)) 
  res}

  
### 8. Selection function
gen.miss <- function(mech,...) {
  switch(mech,
    lod=gen.miss.lod(...),
    trend=gen.miss.trend(...),
    lod.runday=gen.miss.trend.runday(...),
    trend.runday=gen.miss.trend.runday(...),
    mcar=gen.miss.mcar(...))   }
    
### 9. median correction
med.corr <- function(x,rundays,last.pheno=F) {
  meds <- apply(x[,1:(ncol(x)-as.numeric(last.pheno)),drop=F],2,function(y) tapply(y,rundays,median))
  meds <- data.frame(meds)
  meds <- meds[rundays,]
  if(last.pheno) meds <- data.frame(meds,0)
  x-meds }
  