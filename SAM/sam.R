rm(list=ls())
library(stockassessment)
fit = fitfromweb("NEACod-2023")
ssbplot(fit)
fbarplot(fit, partial = FALSE)

library(RTMB)
load("SAM/NEAcod2023.RData")

#Define SSB and F-bar functions
ssbFUN <- function(logN, logFF, M, SW, MO, PF, PM){
  nrow <- nrow(logN)
  ncol <- ncol(logN)
  ret <- numeric(nrow)
  for(y in 1:nrow){
    for(a in 1:ncol){
      ret[y] = ret[y]+SW[y,a]*MO[y,a]*exp(logN[y,a])*exp(-PF[y,a]*exp(logFF[y,a])-PM[y,a]*M[y,a])
    }
  }
  return(ret);
}

fbarFUN <- function(logFF,FbarRange,minAge){
  nrow <- nrow(logFF)
  ageI <- c(FbarRange[1]:FbarRange[2])-(minAge -1)
  ncol <- ncol(logFF)
  ret <- numeric(nrow)
  for(y in 1:nrow){
    tt <- 0
    for(a in ageI){
      tt <- tt + exp(logFF[y,a])
    }
    ret[y] = ret[y]+tt/length(ageI)
  }
  return(ret);
}

f<-function(par){
  getAll(par, dat)

  sdR <- exp(logSdLogN[1])
  sdS <- exp(logSdLogN[-1])
  sdF <- exp(logSdLogFsta[keyVarF[1,unique(keyF+1)]+1])
  sdO <- exp(logSdLogObs)

  logobs[is.na(logobs)] <- missing  ##patch missing
  logFF <- logF[,keyF+1] ## expand F
  
  ssb <- ssbFUN(logN,logFF,M,SW,MO,PF,PM)
  fbar <- fbarFUN(logFF, FbarRange = fbarRange, minAge = minAge)
  nll <- 0

  #Recruitment
  for(y in 2:nrow(M)){
    predN <- logN[y-1,1]
    nll <- nll - dnorm(logN[y,1],predN,sdR,TRUE)
  }

  #Survival process
  for(y in 2:nrow(M)){
    for(a in 2:ncol(M)){
      predN <- logN[y-1,a-1]-exp(logFF[y-1,a-1])-M[y-1,a-1]
      if(a==ncol(M)){
        predN <- log(exp(predN)+exp(logN[y-1,a]-exp(logFF[y-1,a])-M[y-1,a]))
      }
      nll <- nll - dnorm(logN[y,a],predN,sdS,TRUE)
    }
  }

  #Fishing mortality
  SigmaF <- matrix(0,ncol(logF),ncol(logF))
  if(fcormode==0){
    diag(SigmaF) <- sdF*sdF
  }
  for(y in 2:nrow(M)){
    nll <- nll - dmvnorm(logF[y,],logF[y-1,],SigmaF,log=TRUE)
  }
  
  #Observations
  logPred <- numeric(length(logobs))
  for(i in 1:length(logobs)){
    y <- aux[i,1]-minYear+1
    f <- aux[i,2]
    a <- aux[i,3]-minAge+1
    Z <- exp(logFF[y,a])+M[y,a]
    Z_P <- exp(logFF[y,a:dim(logN)[2]])+M[y,a:dim(logN)[2]]
    if(fleetTypes[f]==0){
      logPred[i] <- logN[y,a]-log(Z)+log(1-exp(-Z))+logFF[y,a]
    }
    if(fleetTypes[f]==2){
      if(a < maxAgePerFleet[f] - minAge + 1){
        logPred[i] <- logFpar[keyQ[f,a]+1]+logN[y,a]-Z*sampleTimes[f]  
      }else{
        logPred[i] <- logFpar[keyQ[f,a]+1]+
          log(sum(exp(logN[y,a:dim(logN)[2]]-Z_P*sampleTimes[f])))
      }
    }
  }
  Svec <- list()
  for(f in 1:nrow(idx1)){
    thisdim <- sum(!is.na(keySd)[f,])
    S <- matrix(0,thisdim,thisdim)
    if(covType[f]==0){
      diag(S) <- sdO[na.omit(keySd[f,])+1]^2
    }else if(covType[f]==1){
      dist <- numeric(thisdim)
      d=2;
      for(a in 1:ncol(keyIGAR)){
        if(!is.na(keyIGAR[f,a])){
          dist[d] <- dist[d-1]+exp(transfIRARdist[keyIGAR[f,a]+1])
          d <- d+1
        }
      }
      sdvec <- sdO[na.omit(keySd[f,])+1]
      for(i in 1:nrow(S)){
        for(j in 1:(i-1)){
          S[i,j] <- sdvec[i]*sdvec[j]*(0.5^(dist[i]-dist[j]));
          S[j,i] <- S[i,j];
        }
        S[i,i] <- sdvec[i]^2;
      }
    }
    Svec[[f]] <- S;
  }
  for(f in 1:nrow(idx1)){
    for(y in 1:ncol(idx1)){
      if(!is.na(idx1[f,y])){
        idx <- (idx1[f,y]+1):(idx2[f,y]+1)
        nll <- nll -  dmvnorm(logobs[idx],logPred[idx],Svec[[f]],log=TRUE)
      }
    }
  }
  logSSB = log(ssb)
  logFbar = log(fbar)
  ADREPORT(logSSB)
  ADREPORT(logFbar)
  nll
}


#run the model
obj <- RTMB::MakeADFun(f, par, random=c("logN", "logF", "missing"))
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(eval.max=1000, iter.max=1000,trace = 1))
sdr<-sdreport(obj)

fit$opt$objective-opt$objective

#Get SSB and Fbar
rl<-as.list(sdr, "Est",report = TRUE)
rlsd<-as.list(sdr, "Std",report = TRUE)

ssbplot(fit)
lines(dat$years, exp(rl$logSSB),col = 'red',lwd = 2)
lines(dat$years, exp(rl$logSSB + 2*rlsd$logSSB),col = 'red',lwd = 2,lty = 2)
lines(dat$years, exp(rl$logSSB - 2*rlsd$logSSB),col = 'red',lwd = 2,lty = 2)

fbarplot(fit,partial = FALSE)
lines(dat$years, exp(rl$logFbar),col = 'red',lwd = 2)
lines(dat$years, exp(rl$logFbar + 2*rlsd$logFbar),col = 'red',lwd = 2,lty = 2)
lines(dat$years, exp(rl$logFbar - 2*rlsd$logFbar),col = 'red',lwd = 2,lty = 2)
