# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# R functions to impute and then model "metabolomics" data with missing values using different imputation strategies
#
# Parts of this are experimental
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #




### 1. nie: conduct RC method as conducted by Nie et al.
nie <- function(x,                                    # normally distributed covariate vector
                tr=min(x,na.rm=T)-1E-16,              # truncation point, it is assumed that missings in x if values are below tr
                return.exp = F) {                     # logical, indicating whether only E(X|X<tr) (TRUE) or the imputed data vector (FALSE) should be returned
  xcc <- x[!is.na(x)]
  fit <- tryCatch(mle.tmvnorm(matrix(xcc),lower=tr,start=list(mu=mean(xcc),sigma=matrix(var(xcc)))),error=function(e) NULL)
  if(!is.null(fit)) {
    su <- fit@coef
    mu <- su[1]
    sigma <- su[2]
    xmis <- mu-sigma*dnorm((tr-mu)/sigma)/pnorm((tr-mu)/sigma)
    if(is.na(xmis)) xmis <- tr} else xmis <- NA       
  if(return.exp) ret <- xmis
  if(!return.exp) {
    x[is.na(x)] <- xmis
    ret <- x}
  return(ret)}
  
  
### 2. trunc_uni: conduct univariate TSI, i.e., impute missing values of a single variable x by drawing from the truncated normal distribution  
trunc_uni <- function(x,                              # normally distributed covariate vector
                      tr=min(x,na.rm=T)-1E-16,        # truncation point, it is assumed that missings in x if values are below tr
                      seed=NULL,                      # seed for random number generator
                      outlier.treatment =T,
                      outlier.back=T,                 # after imputation, replace outlier positions by their original values
                      algorithm="gibbs") {            # algorithm in rtmvnorm
  if(any(is.na(x))) {
    if(!is.null(seed)) set.seed(seed)
  
    # remove outliers before imputation (to increase stability)
    Lower <- mean(x,na.rm=T) - 4*sd(x,na.rm=T)
    Upper <- mean(x,na.rm=T) + 4*sd(x,na.rm=T)
    if(outlier.treatment) {
      outliers <- which(x>Upper | x<Lower)
      if(length(outliers)>=1) {
        save.x.outliers <- x[outliers]
        x[outliers] <- NA
        tr <- min(x,na.rm=T)-1E-16} } 

    xcc <- x[!is.na(x)]
    fit <- tryCatch(mle.tmvnorm(matrix(xcc),lower=tr,start=list(mu=mean(xcc),sigma=matrix(var(xcc))),lower.bounds=Lower, upper.bounds=Upper,method="L-BFGS-B"),error=function(e) NULL)
    if(!is.null(fit)) {
      su <- fit@coef
      mu <- su[1]
      sigma <- su[2]
      nmis <- length(which(is.na(x)))
      if(algorithm=="gibbs") xmis <- rtmvnorm(nmis, mean=mu, sigma=diag(sigma^2,1), upper=tr, algorithm="gibbs",burn.in.samples=100) 
      if(algorithm=="rejection") xmis <- rtmvnorm(nmis, mean=mu, sigma=diag(sigma^2,1), upper=tr, algorithm="rejection")
      x[is.na(x)] <- xmis }   }  # if error, returns still with missings
      
  # remove -infinite values
  x[x==-Inf] <- min(x[x>-Inf])
  if(outlier.treatment & outlier.back)  if(length(outliers)>=1) x[outliers] <- save.x.outliers 
  return(x)}
  
### 3.a) trunc_uni.stab
trunc_uni.stab <- function(...)   {
  tryCatch(trunc_uni(...), error=function(e) NULL)
}
 
### 4. Functions from R package miceadds 
.milist <- function (mi.res){
    mi.list <- NULL
    M <- mi.res$m
    for (ii in 1:M) {
        mi.list[[ii]] <- complete(mi.res, action = ii)
    }
    return(mi.list)
}


.sub.micombine.cor <- function (cor.list, N, conf.level){
    fisher.cor.list <- as.list(1/2 * log((1 + cor.list)/(1 -
        cor.list)))
    var.fisher <- as.list(rep(1/(N - 3), length(cor.list)))
    fisher.cor.combine <- mitools:::MIcombine(fisher.cor.list, var.fisher)
    zr <- coef(fisher.cor.combine)
    zr.se <- sqrt(fisher.cor.combine$variance)[1, 1]
    t.zr <- zr/zr.se
    fisher2cor <- function(z) {
        (exp(2 * z) - 1)/(exp(2 * z) + 1)
    }
    res <- c(r = fisher2cor(zr), fisher_r = zr, fisher_rse = zr.se,
        fmi = fisher.cor.combine$missinfo, t = t.zr, p = 2 *
            pnorm(abs(t.zr), lower.tail = FALSE), fisher2cor(zr +
            qnorm((1 - conf.level)/2) * zr.se), fisher2cor(zr -
            qnorm((1 - conf.level)/2) * zr.se))
    names(res)[7] <- paste("lower", round(100 * conf.level, 2),
        sep = "")
    names(res)[8] <- paste("upper", round(100 * conf.level, 2),
        sep = "")
    res <- c(res, -(res[8] - res[7])/(2 * qnorm((1 - conf.level)/2)))
    names(res)[9] <- "rse"
    res <- res[c(1, 9, 2:8)]
    return(res)
}


micombine.cor.mod <- function (mi.res=NULL,mi.list=NULL, variables = 1:(ncol(mi.list[[1]])), conf.level = 0.95) {
    if(is.null(mi.list)) mi.list <- .milist(mi.res)
    N <- nrow(mi.list[[1]])
    VV <- length(variables)
    if (is.character(variables)) {
        variables <- which(colnames(mi.list[[1]]) %in% variables)
    }
    dfr <- NULL
    for (i in 1:(VV - 1)) {
        for (j in (i + 1):VV) {
            if (i != j) {
                ii <- variables[i]
                jj <- variables[j]
                if (i != j) {
                  cor.ii.jj <- unlist(lapply(mi.list, FUN = function(dat) {
                    cor(dat[, ii], dat[, jj])
                  }))
                  res.ii.jj <- .sub.micombine.cor(cor.list = cor.ii.jj,
                    N = N, conf.level = conf.level)
                  dfr <- rbind(dfr, c(ii, jj, res.ii.jj))
                }
            }
        }
    }
    vars <- colnames(mi.list[[1]])
    dfr1 <- dfr
    dfr <- rbind(dfr, dfr1[, c(2, 1, seq(3, ncol(dfr)))])
    if (VV == 2) {
        dfr <- dfr[, -c(1:2)]
    }
    else {
        dfr <- data.frame(variable1 = vars[dfr[, 1]], variable2 = vars[dfr[,
            2]], dfr[, -c(1:2)])
    }
    dfr
}

  
### 5. do.imp
do.imp <- function(dat,                         # data frame, 1. col: met1, 2. col: met2, 3. col: continous pheno w/o missings
                   datCom=NULL,                 # potentially corresponding complete data set
                   methods = c("complete","cca", "min","mean","nie","u.tsi","mi.avg.norm", "mi.avg.pmm",           
                               "u.tsi.runday","si.norm","si.pmm","si.pan.incl.runday",
                               "si.pmm.in.rundays","knn.sample","knn.variable",
                               "u.tsmi", "u.tsmi.runday",  "mi.norm","mi.pmm","mi.pan.incl.runday"),
                   analysis=c("mse1","mse2","mse3","mse4","corxy","pcorxy","lmx","lmy","logx"),  
                   m=20,                        # number of imputations for multiple imputation
                   which.aux=3:12,              # which are aux. variables or NULL
                   which.y = max(which.aux)+1,  # which is the response variable   
                   do.quickpred = F,            # logical indicating whether quickpred function should be used to create sparser pred matrix for mice() 
                   K=5,                         # vector of values K for KNN
                   min.in.runday=10,            # minimum number of observed values in each runday used as condition to perform TSI
                   return.da1=F) {              # logical indicating whether da1 should be returned (imputed data using only x) instead of regression/correlation results (e.g. for distribution assessment)


### Function to get standard error for correlation estimate

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
cor.test.plus <- function(x) {
  list(x, 
       SE = unname(sqrt((1 - x$estimate^2)/x$parameter)))
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


### Functions to perform analysis for single and multiple imputations

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  do.analysis.one <- function(da1,da2,da3,meth) { 
    if("mse1" %in% analysis)   mse1 <- tryCatch(c(mean((da1[is.na(dat[,1]),1]-datCom[is.na(dat[,1]),1])^2),NA,NA),error=function(e) rep(NA,3))    else mse1 <- rep(NA,3)
    if("mse2" %in% analysis)   mse2 <- tryCatch(c(mean((da1[is.na(dat[,2]),2]-datCom[is.na(dat[,2]),2])^2),NA,NA),error=function(e) rep(NA,3))    else mse2 <- rep(NA,3)
    if("mse3" %in% analysis)   mse3 <- tryCatch(c(mean((da2[is.na(dat[,1]),1]-datCom[is.na(dat[,1]),1])^2),NA,NA),error=function(e) rep(NA,3))    else mse3 <- rep(NA,3)
    if("mse4" %in% analysis)   mse4 <- tryCatch(c(mean((da3[is.na(dat[,1]),1]-datCom[is.na(dat[,1]),1])^2),NA,NA),error=function(e) rep(NA,3))    else mse4 <- rep(NA,3)
    if("corxy" %in% analysis) corxy <- tryCatch(unlist(cor.test.plus(cor.test(da1[,1],da1[,2],use="p")))[c("estimate.cor","SE","p.value")],error=function(e) rep(NA,3))    else corxy <- rep(NA,3)
    if("pcorxy" %in% analysis) pcorxy <- tryCatch(sapply(pcor(da1[,c(1,2,which.aux)])[c("estimate","statistic","p.value")], function(part) part[2,1]),error=function(e) rep(NA,3))    else pcorxy <- rep(NA,3)
    if("lmx" %in% analysis)   lmx <- tryCatch(summary(lm(da2[,which.y] ~ da2[,1]))$coef["da2[, 1]",c(1,2,4)],error=function(e)rep(NA,3))                  else lmx <- rep(NA,3)
    if("lmy" %in% analysis)   lmy <- tryCatch(summary(lm(da2[,1] ~ da2[,which.y]))$coef["da2[, which.y]",c(1,2,4)],error=function(e) rep(NA,3))                  else lmy <- rep(NA,3)
    if("logx" %in% analysis)   logx <- tryCatch(summary(glm(I(as.numeric(da3[,which.y]>median(da3[,which.y]))) ~ da3[,1],family=binomial))$coef["da3[, 1]",c(1,2,4)],error=function(e) rep(NA,3)) else logx <-rep(NA,3) 
    ret <-rbind(mse1=c("mse_1cor",mse1), mse2=c("mse_2cor",mse2), mse3=c("mse_1lm",mse3), mse4=c("mse_1log",mse4),corxy=c("corxy",corxy),pcorxy=c("pcorxy",pcorxy),lmx=c("lmx",lmx),lmy=c("lmy",lmy),logx=c("logx",logx))
    colnames(ret) <- c("analysis","estimate","SE","pvalue")
    ret <- data.frame(method=meth,ret)
    for (i in 1:2) ret[,i] <- as.character(ret[,i])
    for (i in 3:5) ret[,i] <- as.numeric(as.character(ret[,i]))
    return(ret) }
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  do.analysis.mult1 <- function(da1list,da2list,da3list,meth) {     
    if("mse1" %in% analysis)   mse1 <- tryCatch(c(mean(sapply(da1list,function(da1) mean((da1[is.na(dat[,1]),1]-datCom[is.na(dat[,1]),1])^2))),NA,NA),error=function(e) rep(NA,3))  else mse1 <- rep(NA,3)
    if("mse2" %in% analysis)   mse2 <- tryCatch(c(mean(sapply(da1list,function(da1) mean((da1[is.na(dat[,2]),2]-datCom[is.na(dat[,2]),2])^2))),NA,NA),error=function(e) rep(NA,3))  else mse2 <- rep(NA,3)
    if("mse3" %in% analysis)   mse3 <- tryCatch(c(mean(sapply(da2list,function(da2) mean((da2[is.na(dat[,1]),1]-datCom[is.na(dat[,1]),1])^2))),NA,NA),error=function(e) rep(NA,3))  else mse3 <- rep(NA,3)
    if("mse4" %in% analysis)   mse4 <- tryCatch(c(mean(sapply(da3list,function(da3) mean((da3[is.na(dat[,1]),1]-datCom[is.na(dat[,1]),1])^2))),NA,NA),error=function(e) rep(NA,3))  else mse4 <- rep(NA,3)
    if("corxy" %in% analysis) corxy <- tryCatch(micombine.cor.mod(mi.list=da1list,variables=1:2) [1,c("r","rse","p")],error=function(e) rep(NA,3))    else corxy <- rep(NA,3) 
    if("pcorxy" %in% analysis) pcorxy <- tryCatch(.sub.micombine.cor(na.omit(sapply(da1list,function(da1) 
           pcor(da1[,c(1,2,which.aux)])$estimate[2,1])),N=nrow(dat),conf.level=0.95)[c("r","rse","p")],error=function(e) rep(NA,3))    else pcorxy <- rep(NA,3)      
    if("lmx" %in% analysis)   lmx <- tryCatch(summary(pool(as.mira(lapply(da2list,function(da2) lm(da2[,which.y] ~ da2[,1])))))["da2[, 1]",c(1,2,5)],error=function(e) rep(NA,3))            else lmx <- rep(NA,3)
    if("lmy" %in% analysis)   lmy <- tryCatch(summary(pool(as.mira(lapply(da2list,function(da2) lm(da2[,1] ~ da2[,which.y])))))["da2[, which.y]",c(1,2,5)],error=function(e) rep(NA,3))            else lmy <- rep(NA,3)
    if("logx" %in% analysis)   logx <- tryCatch(summary(pool(as.mira(lapply(da3list,function(da3) glm(I(as.numeric(da3[,which.y]>median(da3[,which.y]))) ~ da3[,1],family=binomial)))))["da3[, 1]",c(1,2,5)],error=function(e) rep(NA,3))            else logx <- rep(NA,3)
    ret <-rbind(mse1=c("mse_1cor",mse1), mse2=c("mse_2cor",mse2), mse3=c("mse_1lm",mse3), mse4=c("mse_1log",mse4),corxy=c("corxy",corxy),pcorxy=c("pcorxy",pcorxy),lmx=c("lmx",lmx),lmy=c("lmy",lmy),logx=c("logx",logx))
    colnames(ret) <- c("analysis","estimate","SE","pvalue")
    ret <- data.frame(method=meth,ret)
    for (i in 1:2) ret[,i] <- as.character(ret[,i])
    for (i in 3:5) ret[,i] <- as.numeric(as.character(ret[,i]))
    return(ret) }
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

### Function to perform mice with some error handling

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
mice.stab <- function(X, m, method, force.in=0) {
  if(!any(apply(X,2,function(x) length(which(!is.na(x)))<2))) {
    if(do.quickpred) {
      pred <- quickpred(X,mincor=0.1, minpuc=0.25, include=names(X)[force.in], method="pearson") 
      if(method=="2l.pan") pred[,"rundays"] <- -2
      ret <- tryCatch(mice(X,pred=pred,m=m,method=method,printFlag=F),error=function(e) NULL) } else {
      if(method=="2l.pan") {
        pred <- mice(X,maxit=0)$pred
        pred[,"rundays"] <- -2
        ret <- tryCatch(mice(X,m=m,method=method,pred=pred,printFlag=F),error=function(e) NULL) } else
      ret <- tryCatch(mice(X,m=m,method=method,printFlag=F),error=function(e) NULL) } }     else ret <- NULL 
    ret} 

complete.stab <- function(miceout, j, dat) {
  if(!is.null(miceout)) ret <- complete(miceout,j) else  ret <- dat
  ret}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


if(!return.da1) results <- NULL else results <- list()

  incom.vars <- which(apply(dat[,c(1,2,which.aux)],2,function(x) any(is.na(x))))

if(do.quickpred) {
 force.in1 <- 1:2
 force.in2 <- force.in3 <- c(1,2+length(which.aux))    }  else {
 force.in1 <- force.in2 <- force.in3 <- 0 }



  ## complete
  if(!is.null(datCom) & "complete" %in% methods) {
    dat_ <- datCom
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=dat_,da2=dat_,da3=dat_,meth="complete")) else results <- c(results,complete=list(dat_))
    rm(dat_)}

  ## cca
  if("cca" %in% methods) {
    dat_ <- dat
    if(!return.da1) {
      da1 <- dat_[complete.cases(dat_[,1:2]),]
      da2 <- dat_[complete.cases(dat_[,1]),]
      results <- rbind(results,do.analysis.one(da1=da1,da2=da2,da3=da2,meth="cca"))} else results <- c(results,cca=list(dat_))
    rm(dat_)}

  ## min
  if("min" %in% methods)  {
    dat_ <- dat
    for (j in incom.vars) dat_[,j][is.na(dat_[,j])] <- min(dat_[,j],na.rm=T)   
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=dat_,da2=dat_,da3=dat_,meth="min")) else results <- c(results,min=list(dat_))
    rm(dat_)}

  ## mean
  if("mean" %in% methods)  {
    dat_ <- dat
    for (j in incom.vars) dat_[,j][is.na(dat_[,j])] <- mean(dat_[,j],na.rm=T)
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=dat_,da2=dat_,da3=dat_,meth="mean")) else results <- c(results,mean=list(dat_))
    rm(dat_)}
  
  if("mi.avg.norm" %in% methods)  {
    dat_ <- dat[,c(1:2,which.aux)]
    imp <- mice(dat_, method="norm", predictorMatrix=quickpred(dat_))
    for(i in 1:length(imp$imp))
    {
      if(length(imp$imp[[i]]>0))
      {
        for(j in 1:nrow(imp$imp[[i]]))
        {
          dat_[row.names(dat_)==names(apply(imp$imp[[i]],1,mean))[j],i]<-apply(imp$imp[[i]],1,mean)[j]
        }}
    }
    dat_ <- cbind(dat_, dat[,which.y])
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=dat_,da2=dat_,da3=dat_,meth="mi.avg.norm")) else results <- c(results,mi.avg.norm=list(dat_))
    rm(dat_)}
  
  if("mi.avg.pmm" %in% methods)  {
    dat_ <- dat[,c(1:2,which.aux)]
    imp <- mice(dat_, method="pmm", predictorMatrix=quickpred(dat_))
    for(i in 1:length(imp$imp))
    {
      if(length(imp$imp[[i]]>0))
      {
        for(j in 1:nrow(imp$imp[[i]]))
        {
          dat_[row.names(dat_)==names(apply(imp$imp[[i]],1,mean))[j],i]<-apply(imp$imp[[i]],1,mean)[j]
        }}
    }
    dat_ <- cbind(dat_, dat[,which.y])
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=dat_,da2=dat_,da3=dat_,meth="mi.avg.pmm")) else results <- c(results,mi.avg.pmm=list(dat_))
    rm(dat_)}

  ## nie
  if("nie" %in% methods) {
    dat_ <- dat
    for (j in incom.vars) dat_[,j] <-  nie(dat_[,j])
    if(any(is.na(dat_[,1]))) dat_ <- NULL
    da1 <- dat_
    if(!is.null(dat_)) if(any(is.na(dat_[,2]))) da1 <- NULL
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=da1,da2=dat_,da3=dat_,meth="nie")) else results <- c(results,nie=list(da1))
    rm(da1,dat_)}

  ## nie.runday
  if("nie.runday" %in% methods) {
    dat_ <- dat
    for (j in incom.vars) {
      counts <- tapply(dat[,j], dat_$rundays, function(x) length(which(!is.na(x)))>min.in.runday) 
      expvalues <- NULL    
      for (r in names(counts)[counts]) {
        res <- nie(dat_[dat_$rundays==r,j])
        expvalues <- c(expvalues,nie(dat_[dat_$rundays==r,j], return.exp=T)) 
        if(!is.null(res)) dat_[dat_$rundays==r,j] <- res} 
      if(any(is.na(dat_[,j]))) {
        if(length(expvalues)>0) {if(any(!is.na(expvalues))) dat_[is.na(dat_[,j]),j] <- mean(expvalues,na.rm=T)} else 
          dat_[is.na(dat_[,j]),j] <- nie(dat_[,j],return.exp=T) }}
    if(any(is.na(dat_[,1]))) dat_ <- NULL
    da1 <- dat_
    if(!is.null(dat_)) if(any(is.na(dat_[,2]))) da1 <- NULL
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=da1,da2=dat_,da3=dat_,meth="nie.runday")) else results <- c(results,nie.runday=list(da1))
    rm(da1,dat_)}

  ## u.tsi
  if("u.tsi" %in% methods) {
    dat_ <- dat
    for (j in incom.vars) {
      res <- trunc_uni.stab(dat_[,j])   
      if(!is.null(res)) dat_[,j] <- res}
    if(any(is.na(dat_[,1]))) dat_ <- NULL
    da1 <- dat_
    if(!is.null(dat_)) if(any(is.na(dat_[,2]))) da1 <- NULL
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=da1,da2=dat_,da3=dat_,meth="u.tsi")) else results <- c(results,u.tsi=list(da1))
    rm(da1,dat_)}

  ## u.tsi.runday
  if("u.tsi.runday" %in% methods) {
    dat_ <- dat
    for (j in incom.vars) {
      counts <- tapply(dat[,j], dat_$rundays, function(x) length(which(!is.na(x)))>min.in.runday)      
      for (r in names(counts)[counts]) {
        res <- trunc_uni.stab(dat_[dat_$rundays==r,j])
        if(!is.null(res)) dat_[dat_$rundays==r,j] <- res} }

    da1 <- da2 <- da3 <- dat_
    if("corxy" %in% analysis) {
      if(any(is.na(da1[,1:2]))) {
        miceout <- mice.stab(da1[,c(1:2,which.aux)],m=1,method="norm",force.in=force.in1)     
        da1[,c(1:2,which.aux)] <- complete.stab(miceout,1,da1) } }
    if("lmx" %in% analysis | "lmy" %in% analysis) {
      if(any(is.na(da2[,1]))) {
        miceout <- mice.stab(dat_[,c(1,which.aux,which.y)],m=1,method="norm",force.in=force.in2)   
        da2[,c(1,which.aux,which.y)] <- complete.stab(miceout,1,da2) } }    
    if("logx" %in% analysis) {
      if(any(is.na(da3[,1]))) {       
        miceout <- mice.stab(data.frame(da3[,c(1,which.aux),drop=F], as.numeric(da3[,which.y]>median(da3[,which.y]))),m=1,method="norm",force.in=force.in3) 
        da3[,c(1,which.aux)] <- complete.stab(miceout,1,da3)[,1:length(c(1,which.aux))] } }   
    if(any(is.na(da1[,1:2]))) da1 <- NULL
    if(any(is.na(da2[,1]))) da2 <- NULL
    if(any(is.na(da3[,1]))) da3 <- NULL    
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=da1,da2=da2,da3=da3,meth="u.tsi.runday")) else results <- c(results, u.tsi.runday=list(da1))
    rm(da1,da2,da3,dat_)}

  ## si.norm und si.pmm
  for (method in c("si.norm","si.pmm")) { 
    if(method %in% methods) {
      dat_ <- da1 <- da2 <- da3 <- dat
      meth <- sub("si.","",sub("mi.","",sub(".incl.runday","",sub(".in.rundays","",method))))
      if("corxy" %in% analysis) {
        if(any(is.na(da1[,1:2]))) {
          miceout <- mice.stab(da1[,c(1:2,which.aux)],m=1,method=meth,force.in=force.in1) 
          da1[,c(1:2,which.aux)] <- complete.stab(miceout,1,da1) } }     
      if("lmx" %in% analysis | "lmy" %in% analysis) {
        if(any(is.na(da2[,1]))) {
          miceout <- mice.stab(da2[,c(1,which.aux,which.y)],m=1,method=meth,force.in=force.in2) 
          da2[,c(1,which.aux,which.y)] <- complete.stab(miceout,1,da2) } }     
      if("logx" %in% analysis) {
        if(any(is.na(da3[,1]))) {
          miceout <- mice.stab(data.frame(da3[,c(1,which.aux),drop=F], as.numeric(da3[,which.y]>median(da3[,which.y]))),m=1,method=meth,force.in=force.in3) 
          da3[,c(1,which.aux)] <- complete.stab(miceout,1,da3)[,1:length(c(1,which.aux))] } }     
    if(any(is.na(da1[,1:2]))) da1 <- NULL
    if(any(is.na(da2[,1]))) da2 <- NULL
    if(any(is.na(da3[,1]))) da3 <- NULL
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=da1,da2=da2,da3=da3,meth=method)) else {
      results <- c(results, list(da1))
      names(results)[length(results)] <- method}
    rm(da1,da2,da3,dat_)} }


  ## si.pan.incl.runday  
  if("si.pan.incl.runday" %in% methods) {
      dat_ <- dat
      meth <- "2l.pan"
      da1 <- da2 <- da3 <- dat_ 
      if("corxy" %in% analysis) {
        if(any(is.na(da1[,1:2]))) {
          miceout <- mice.stab(da1[,c(1:2,which.aux,ncol(da1))],m=1,method=meth,force.in=c(force.in1,3+length(which.aux))) 
          da1[,c(1:2,which.aux,ncol(da1))] <- complete.stab(miceout,1,da1) }   } 
      if("lmx" %in% analysis | "lmy" %in% analysis) {
        if(any(is.na(da2[,1]))) {
          miceout <- mice.stab(da2[,c(1,which.aux,which.y,ncol(da2))],m=1,method=meth,force.in=c(force.in2,3+length(which.aux))) 
          da2[,c(1,which.aux,which.y,ncol(da2))] <- complete.stab(miceout,1,da2) }}     
      if("logx" %in% analysis) {
        if(any(is.na(da3[,1]))) {
          miceout <- mice.stab(data.frame(da3[,c(1,which.aux,ncol(da3)),drop=F], as.numeric(da3[,which.y]>median(da3[,which.y]))),m=1,method=meth, force.in=c(force.in3,3+length(which.aux)))   
          da3[,c(1,which.aux,ncol(da3))] <- complete.stab(miceout,1,da3)[,1:length(c(1,which.aux,ncol(da3)))] }  }     
    if(any(is.na(da1[,1:2]))) da1 <- NULL
    if(any(is.na(da2[,1]))) da2 <- NULL
    if(any(is.na(da3[,1]))) da3 <- NULL
    if(!return.da1) results <- rbind(results,do.analysis.one(da1=da1,da2=da2,da3=da3,meth="si.pan.incl.runday")) else {
      results <- c(results, si.pan.incl.runday=list(da1))}
    rm(da1,da2,da3,dat_)} 

    
  ## knn.variable  
  if("knn.variable" %in% methods) {
    da1list <- da2list <- rep(list(dat),length(K)) 
  
    if("corxy" %in% analysis) { 
      da1 <- dat[,c(1:2,which.aux)]      
      D2 <- as.matrix(dist(t(scale(da1)),upper=T,diag=T))
      D2[D2==0] <- NA
      for (j in incom.vars) {
          comobs <- complete.cases(da1[,j])
          Mean <- mean(da1[,j],na.rm=T)
          SD <- sd(da1[,j],na.rm=T)
          if(any(!is.na(D2[,j]))) {
            KNNvars <- order(D2[,j],na.last=NA)
            KNNvars <- KNNvars[sapply(KNNvars, function(jj) any(!is.na(da1[!comobs,jj])))] }    else KNNvars <- NULL
          
          da1list <- lapply(1:length(da1list),function(ii) {
            k <- K[ii]
            da <- da1list[[ii]]
            KNNvars_sel <- KNNvars[1:min(k,length(KNNvars))]
            if(any(!is.na(D2[,j])) & length(KNNvars)>=1) da[!comobs,c(1:2,which.aux)[j]]  <- 
              sapply(1:length(which(!comobs)),function(co) {
                da_all <- scale(da)[co,c(1:2,which.aux)[KNNvars]]
                da_sel <- scale(da)[co,c(1:2,which.aux)[KNNvars_sel]]
                if(any(!is.na(da_sel))) ret <- sum((da_sel*SD+Mean)*exp(-D2[KNNvars_sel,j]),na.rm=T)/sum(exp(-D2[KNNvars_sel,j])[!is.na(da_sel)],na.rm=T) else ret <- na.omit(da_all)[1]*SD+Mean
                ret})
            if(any(is.na(da[,j]))) {      
               still.incom <-  !complete.cases(da[,j]) 
               da[still.incom,j] <- mean(da[!still.incom,j],na.rm=T) }
     da}  ) } }

        
    if("lmx" %in% analysis | "lmy" %in% analysis | "logx" %in% analysis) {
      da2 <- dat[,c(1,which.aux)]     
      incom.vars2 <- which(apply(dat[,c(1,which.aux)],2,function(x) any(is.na(x))))
      D2 <- as.matrix(dist(t(scale(da2)),upper=T,diag=T))
      D2[D2==0] <- NA
      for (j in incom.vars2) {
          comobs <- complete.cases(da2[,j])
          Mean <- mean(da2[,j],na.rm=T)
          SD <- sd(da2[,j],na.rm=T)
          if(any(!is.na(D2[,j]))) {
            KNNvars <- order(D2[,j],na.last=NA)
            KNNvars <- KNNvars[sapply(KNNvars, function(jj) any(!is.na(da2[!comobs,jj])))] }    else KNNvars <- NULL
          
          da2list <- lapply(1:length(da2list),function(ii) {
            k <- K[ii]
            da <- da2list[[ii]]
            KNNvars_sel <- KNNvars[1:min(k,length(KNNvars))]
            if(any(!is.na(D2[,j])) & length(KNNvars)>=1) da[!comobs,c(1,which.aux)[j]]  <- 
              sapply(1:length(which(!comobs)),function(co) {
                da_all <- scale(da)[co,c(1,which.aux)[KNNvars]]
                da_sel <- scale(da)[co,c(1,which.aux)[KNNvars_sel]]
                if(any(!is.na(da_sel))) ret <- sum((da_sel*SD+Mean)*exp(-D2[KNNvars_sel,j]),na.rm=T)/sum(exp(-D2[KNNvars_sel,j])[!is.na(da_sel)],na.rm=T) else ret <- na.omit(da_all)[1]*SD+Mean
                ret})            
            if(any(is.na(da[,c(1,which.aux)[j]]))) {      
               still.incom <-  !complete.cases(da[,c(1,which.aux)[j]]) 
               da[still.incom,c(1,which.aux)[j]] <- mean(da[!still.incom,c(1,which.aux)[j]],na.rm=T) }
            da}  ) } }
 
    if(!return.da1) for(ii in 1:length(K)) {
      k <- K[ii]
      results <- rbind(results,do.analysis.one(da1=da1list[[ii]],da2=da2list[[ii]],da3=da2list[[ii]],meth=paste("knn.variable_",k,sep="")))}    else {
      results <- c(results,list(da1list))
      names(results)[length(results)] <- "knn.variable" }
    rm(da1list,da2list) }          

  ## knn.sample.euc
  if("knn.sample.euc" %in% methods) {
  
    da1list <- da2list <- rep(list(dat),length(K))
  
    if("corxy" %in% analysis) { 
      da1 <- dat[,c(1:2,which.aux)] 
      incom.obs <- which(apply(da1,1,function(x) any(is.na(x))))
      D2 <- as.matrix(dist(scale(da1),upper=T,diag=T))
      D2[D2==0] <- NA
      for (i in incom.obs){                                   
          comvars <-  complete.cases(as.numeric(da1[i,]))
          if(any(!is.na(D2[i,]))) {
            KNNids <- order(D2[i,],na.last=NA) 
            KNNids <- KNNids[sapply(KNNids,function(j) any(!is.na(da1[j,!comvars])))]
            }  else KNNids <- KNNids_naomit <- NULL
          da1list <- lapply(1:length(da1list),function(ii) {
            k <- K[ii] 
            da <-  da1list[[ii]]
            if(!is.null(KNNids)) KNNids_sel <- KNNids[1:min(k,length(KNNids))]
            if(any(!is.na(D2[i,])) & length(KNNids)>=1) da[i,c(1:2,which.aux)[!comvars]] <- 
              sapply(1:length(which(!comvars)), function(co) {                
                da_all <- da[KNNids,c(1:2,which.aux)[co]]
                da_sel <- da[KNNids_sel,c(1:2,which.aux)[co]]
                if(any(!is.na(da_sel))) ret <- sum(da_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(da_sel)],na.rm=T) else ret <- na.omit(da_all)[1]
                ret})  
            
            if(any(is.na(da[i,]))) {       
              still.incom <-  !complete.cases(as.numeric(da[i,c(1:2,which.aux)])) 
              da[i,c(1:2,which.aux)[still.incom]] <- apply(da[,c(1:2,which.aux)[still.incom],drop=F],2,mean,na.rm=T)}  
          da}  ) } }
          
          
    if("lmx" %in% analysis | "lmy" %in% analysis | "logx" %in% analysis) {
      da2 <- dat[,c(1,which.aux)]     
      incom.obs <- which(apply(da2,1,function(x) any(is.na(x))))
      D2 <- as.matrix(dist(scale(da2),upper=T,diag=T)) # eucl
      D2[D2==0] <- NA
      for (i in incom.obs){                                   
          comvars <-  complete.cases(as.numeric(da2[i,]))
          if(any(!is.na(D2[i,]))) {
            KNNids <- order(D2[i,],na.last=NA) 
            KNNids <- KNNids[sapply(KNNids,function(j) any(!is.na(da2[j,!comvars])))]
            }   else KNNids <- KNNids_naomit <- NULL
          da2list <- lapply(1:length(da2list),function(ii) {
            k <- K[ii] 
            da <-  da2list[[ii]]
            if(!is.null(KNNids)) KNNids_sel <- KNNids[1:min(k,length(KNNids))]
            if(any(!is.na(D2[i,])) & length(KNNids)>=1) da[i,c(1,which.aux)[!comvars]] <-               
                sapply(1:length(which(!comvars)), function(co) {                
                da_all <- da[KNNids,c(1:2,which.aux)[co]]
                da_sel <- da[KNNids_sel,c(1:2,which.aux)[co]]
                if(any(!is.na(da_sel))) ret <- sum(da_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(da_sel)],na.rm=T) else ret <- na.omit(da_all)[1]
                ret}) 
                
            if(any(is.na(da[i,]))) {     
              still.incom <-  !complete.cases(as.numeric(da[i,c(1,which.aux)])) 
              da[i,c(1,which.aux)[still.incom]] <- apply(da[,c(1,which.aux)[still.incom],drop=F],2,mean,na.rm=T)}   
         da}  ) } }

    if(!return.da1) for(ii in 1:length(K)) {
      k <- K[ii]
      results <- rbind(results,do.analysis.one(da1=da1list[[ii]],da2=da2list[[ii]],da3=da2list[[ii]],meth=paste("knn.sample.euc_",k,sep="")))}    else {
      results <- c(results,list(da1list))
      names(results)[length(results)] <- "knn.sample.euc" }
    rm(da1list,da2list) }          

  ## knn.sample.euc.sel - use only variables with cor>0.1 for distance computation

  if("knn.sample.euc.sel" %in% methods) {
  
    da1list <- da2list <- rep(list(dat),length(K)) 
    cor.cutoff <- 0.2
  
    if("corxy" %in% analysis) { 
      da1 <- dat[,c(1:2,which.aux)]  
      incom.obs <- which(apply(da1,1,function(x) any(is.na(x))))
      Cor <- cor(da1,use="p")
      D2list <- lapply(incom.vars, function(j) {
        varsel <- which(abs(Cor[j,])>cor.cutoff)  
        if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
        if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
        D2 <- as.matrix(dist(scale(da1[,varsel])),upper=T,diag=T) 
        if(any(is.na(D2))) {
          D2a <- as.matrix(dist(scale(da1)),upper=T,diag=T)*sqrt(length(varsel)/ncol(da1)) 
          D2[is.na(D2)] <- D2a[is.na(D2)] }
        diag(D2) <- NA
        D2})
      names(D2list) <- incom.vars
      for (i in incom.obs){
        comvars <-  complete.cases(as.numeric(da1[i,]))
        for (j in which(!comvars)) {
          D2 <- D2list[[as.character(j)]]                                 
          if(any(!is.na(D2[i,]))) {
            KNNids <- order(D2[i,],na.last=NA)
            KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(da1[ii,j])))] 
            }  else KNNids  <- NULL
          da1list <- lapply(1:length(da1list),function(ii) {
            k <- K[ii] 
            da <-  da1list[[ii]]
            if(!is.null(KNNids)) KNNids_sel <- intersect(KNNids[1:min(k,length(KNNids))],KNNids_naomit)
            if(length(KNNids_sel)<1) KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))] else 
              if(length(which(sapply(KNNids_sel,function(ii) !is.na(da1[ii,j])))) < floor(k/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))] 
            if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
                da_sel <- da[KNNids_sel,c(1:2,which.aux)[j]]
                da[i,c(1:2,which.aux)[j]] <- sum(da_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(da_sel)],na.rm=T) }
            da}) }}
      da1list <- lapply(da1list, function(da) {
        da[,c(1:2,which.aux)] <- apply(da[,c(1:2,which.aux)],2, function(x) {
          if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
          x}) 
        da}) } 
          
          
    if("lmx" %in% analysis | "lmy" %in% analysis | "logx" %in% analysis) {
      da2 <- dat[,c(1,which.aux)]     
      incom.obs <- which(apply(da2,1,function(x) any(is.na(x))))
      Cor <- cor(da2,use="p")
      D2list <- lapply(incom.vars, function(j) {
        varsel <- which(abs(Cor[j,])>cor.cutoff)  
        if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
        if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
        D2 <- as.matrix(dist(scale(da2[,varsel])),upper=T,diag=T) 
        D2[D2==0] <- NA
        D2})
      names(D2list) <- incom.vars
      for (i in incom.obs){                                   
        comvars <-  complete.cases(as.numeric(da2[i,]))
        for (j in which(!comvars)) {
          D2 <- D2list[[as.character(j)]]                                 
          if(any(!is.na(D2[i,]))) {
            KNNids <- order(D2[i,],na.last=NA)
            KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(da2[ii,j])))] 
            }  else KNNids  <- NULL
          da2list <- lapply(1:length(da2list),function(ii) {
            k <- K[ii] 
            da <-  da2list[[ii]]
            if(!is.null(KNNids)) KNNids_sel <- intersect(KNNids[1:min(k,length(KNNids))],KNNids_naomit)
            if(length(KNNids_sel)<1) KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))] else 
              if(length(which(sapply(KNNids_sel,function(ii) !is.na(da2[ii,j])))) < floor(k/2) )  KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))] 
            if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
                da_sel <- da[KNNids_sel,c(1,which.aux)[j]]
                da[i,c(1,which.aux)[j]] <- sum(da_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(da_sel)],na.rm=T) }
            da}) }}
      da2list <- lapply(da2list, function(da) {
        da[,c(1,which.aux)] <- apply(da[,c(1,which.aux)],2, function(x) {
          if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
          x})
        da }) } 

    if(!return.da1) for(ii in 1:length(K)) {
      k <- K[ii]
      results <- rbind(results,do.analysis.one(da1=da1list[[ii]],da2=da2list[[ii]],da3=da2list[[ii]],meth=paste("knn.sample.euc.sel_",k,sep="")))}    else {
      results <- c(results,list(da1list))
      names(results)[length(results)] <- "knn.sample.euc.sel" }
    rm(da1list,da2list) }          


  ## u.tsmi                        
  if("u.tsmi" %in% methods) {
    dalist <- rep(list(dat),m)
    dafun <- function(da) {
      for (j in incom.vars) {
        res <- trunc_uni.stab(da[,j]) 
        if(!is.null(res)) da[,j] <- res}

      if(any(is.na(da[,1]))) da <- NULL                                  
      da}
    dalist <- lapply(dalist, function(da) tryCatch(dafun(da),error=function(e) NULL))
    da1list <- lapply(dalist, function(da) {
      if(!is.null(da)) if(any(is.na(da[,2]))) da <- NULL
      da })
    if(!return.da1) results <- rbind(results,do.analysis.mult1(da1list=da1list,da2list=dalist,da3list=dalist,meth="u.tsmi"))  else results <- c(results,u.tsmi=list(da1list))
    rm(dalist,da1list)}


  ## u.tsmi.runday     
  if("u.tsmi.runday" %in% methods) {
    dalist <- da1list <- da2list <- da3list <- rep(list(dat),m)
    dafun <- function(da) {
      for (j in incom.vars) {
        counts <- tapply(da[,j], da$rundays, function(x) length(which(!is.na(x)))>min.in.runday)      
        for (r in names(counts)[counts]) {
          res <- trunc_uni.stab(da[da$rundays==r,j])
          if(!is.null(res)) da[da$rundays==r,j] <- res}}
      da}
    dalist <- lapply(dalist, function(da) tryCatch(dafun(da),error=function(e) NULL))
    if("corxy" %in% analysis) {
      da1list <- lapply(dalist, function(da) {
        if(!is.null(da))    if(any(is.na(da[,1:2]))) {
          miceout.l <- mice.stab(da[,c(1:2,which.aux)], m=1,method="norm",force.in=force.in1)
          da[,c(1:2,which.aux)] <- complete.stab(miceout.l,1,da)}
        da}) }
    if("lmx" %in% analysis | "lmy" %in% analysis) {
      da2list <- lapply(dalist, function(da) {
        if(!is.null(da)) if(any(is.na(da[,1]))) {
          miceout.l <- mice.stab(da[,c(1,which.aux,which.y)], m=1,method="norm",force.in=force.in2)
          da[,c(1,which.aux,which.y)] <- complete.stab(miceout.l,1,da)}
        da})  } 
    if("logx" %in% analysis) {
      da3list <- lapply(dalist, function(da) {
        if(!is.null(da))  if(any(is.na(da[,1]))) {
          miceout.l <- mice.stab(data.frame(da[,c(1,which.aux),drop=F], as.numeric(da[,which.y]>median(da[,which.y]))),m=1,method="norm",force.in=force.in3) 
          da[,c(1,which.aux)] <- complete.stab(miceout.l,1,da)[,1:length(c(1,which.aux))]}
        da}) }
    if(any(sapply(da1list,function(da) is.na(da[,1:2])))) da1list <- NULL
    if(any(sapply(da2list,function(da) is.na(da[,1])))) da2list <- NULL
    if(any(sapply(da3list,function(da) is.na(da[,1])))) da3list <- NULL
    if(!return.da1) results <- rbind(results,do.analysis.mult1(da1list=da1list,da2list=da2list,da3list=da3list,meth="u.tsmi.runday")) else results <- c(results,u.tsmi.runday=list(da1list))
    rm(dalist,da1list,da2list,da3list)}

  ## mi.norm und mi.pmm
  for (method in c("mi.norm","mi.pmm")) { 
    if(method %in% methods) {
      da1list <- da2list <- da3list <- rep(list(dat),m)
      meth <- sub("si.","",sub("mi.","",sub(".incl.runday","",sub(".in.rundays","",method))))
      if("corxy" %in% analysis) {
        if(any(is.na(dat[,1:2]))) {
          miceout <- mice.stab(dat[,c(1:2,which.aux)],m=m,method=meth,force.in=force.in1) 
          da1list <- lapply(1:m, function(l) {
            da1list[[l]][,c(1:2,which.aux)] <- complete.stab(miceout,l,da1list[[l]])
            da1list[[l]]}) }} 
      if("lmx" %in% analysis | "lmy" %in% analysis) {
        if(any(is.na(dat[,1]))) {
          miceout <- mice.stab(dat[,c(1,which.aux,which.y)],m=m,method=meth,force.in=force.in2)
          da2list <- lapply(1:m, function(l) {
            da2list[[l]][,c(1,which.aux,which.y)] <- complete.stab(miceout,l,da2list[[l]])
            da2list[[l]]}) }} 
      if("logx" %in% analysis) {
        if(any(is.na(dat[,1]))) {
          miceout <- mice.stab(data.frame(dat[,c(1,which.aux),drop=F], as.numeric(dat[,which.y]>median(dat[,which.y]))),m=m,method=meth,force.in=force.in3)
          da3list <- lapply(1:m, function(l) {
            da3list[[l]][,c(1,which.aux)] <- complete.stab(miceout,l,da3list[[l]])[,1:length(c(1,which.aux))]
            da3list[[l]]}) }} 
    if(any(sapply(da1list,function(da) is.na(da[,1:2])))) da1list <- NULL
    if(any(sapply(da2list,function(da) is.na(da[,1])))) da2list <- NULL
    if(any(sapply(da3list,function(da) is.na(da[,1])))) da3list <- NULL
    if(!return.da1) results <- rbind(results,do.analysis.mult1(da1list=da1list,da2list=da2list,da3list=da3list,meth=method))  else {
      results <- c(results,list(da1list))
      names(results)[length(results)] <- method}
    rm(da1list,da2list,da3list)} }


  ## mi.pan.incl.runday
  if("mi.pan.incl.runday" %in% methods) {
      dat_ <- dat
      da1list <- da2list <- da3list <- rep(list(dat_),m)
      meth <- "2l.pan"
      if("corxy" %in% analysis) {
        if(any(is.na(dat_[,1:2]))) {
          miceout <- mice.stab(dat_[,c(1:2,which.aux,ncol(dat_))],m=m,method=meth, force.in=c(force.in1,3+length(which.aux)))
          da1list <- lapply(1:m, function(l) {
            da1list[[l]][,c(1:2,which.aux,ncol(dat_))] <- complete.stab(miceout,l,da1list[[l]])
            da1list[[l]]}) }} 
      if("lmx" %in% analysis | "lmy" %in% analysis) {
        if(any(is.na(dat_[,1]))) {
          miceout <- mice.stab(dat_[,c(1,which.aux,which.y,ncol(dat_))],m=m,method=meth, force.in=c(force.in2,3+length(which.aux)))
          da2list <- lapply(1:m, function(l) {
            da2list[[l]][,c(1,which.aux,which.y,ncol(dat_))] <- complete.stab(miceout,l,da2list[[l]])
            da2list[[l]]}) }} 
      if("logx" %in% analysis) {
        if(any(is.na(dat_[,1]))) {
          miceout <- mice.stab(data.frame(dat_[,c(1,which.aux,ncol(dat_)),drop=F], as.numeric(dat_[,which.y]>median(dat_[,which.y]))),m=m,method=meth, force.in=c(force.in3,3+length(which.aux)))
          da3list <- lapply(1:m, function(l) {
            da3list[[l]][,c(1,which.aux,ncol(dat_))] <- complete.stab(miceout,l,da3list[[l]])[,1:length(c(1,which.aux,ncol(dat_)))]
            da3list[[l]]}) }} 
    if(any(sapply(da1list,function(da) is.na(da[,1:2])))) da1list <- NULL
    if(any(sapply(da2list,function(da) is.na(da[,1])))) da2list <- NULL
    if(any(sapply(da3list,function(da) is.na(da[,1])))) da3list <- NULL
    if(!return.da1) results <- rbind(results,do.analysis.mult1(da1list=da1list,da2list=da2list,da3list=da3list,meth="mi.pan.include.runday")) else {
      results <- c(results,mi.pan.incl.runday=list(da1list))}
    rm(dat_,da1list,da2list,da3list)} 


results}
