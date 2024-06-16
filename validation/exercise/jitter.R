library(RTMB)
y = readRDS("validation/exercise/ar1Latent.rds")

data = list(y = y)
par = list(logsd = 0,
           logitRho = 0,
           gamma = rep(0,length(y)))


nll = function(par){
  getAll(par,data)
  sd = exp(logsd)
  rho = 2/(1 + exp(-logitRho))-1

  y = OBS(y)

  nll = 0
  nll = nll -dnorm(gamma[1],0,sqrt(sd*sd/(1-rho^2)),TRUE);
  for(i in 2:length(gamma)){
    nll = nll -dnorm(gamma[i],rho*gamma[i-1],sd,TRUE);
  }


  nll = nll -sum(dpois(y,exp(gamma),TRUE))


  sumGamma = sum(exp(gamma)/length(gamma))
  ADREPORT(sumGamma)
  return(nll)
}


obj = MakeADFun(nll,par,random = "gamma")
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
sdrep = sdreport(obj,getJointPrecision = TRUE)
pl = as.list(sdrep,what = "Est")
plsd = as.list(sdrep,what = "Std")
rl = as.list(sdrep,what = "Est", report = TRUE)
rlsd = as.list(sdrep,what = "Std", report = TRUE)


#Store everyting we need in a list
fit = list(obj = obj,
           opt = opt,
           sdrep = sdrep,
           pl = pl,
           plsd = plsd,
           rl = rl,
           rlsd = rlsd,
           data = data)


#Jitter
jit = function(fit,nojit = 10,sd = 0.2){
  parJit <- lapply(1:nojit, function(i)relist(unlist(par)+rnorm(length(unlist(par)),sd=sd), par))
  data = fit$data
  jj = lapply(parJit, function(p){
    obj = MakeADFun(nll,p,random = "gamma")
    opt <- nlminb(obj$par,obj$fn,obj$gr)
    sdrep = sdreport(obj)
    pl = as.list(sdrep,what = "Est")
    rl = as.list(sdrep,what = "Est",report = TRUE)
    ret = list(opt = opt,
               pl = pl,
               rl = rl)
    ret
  })
  maxabsdiff <- apply(abs(do.call(cbind, lapply(jj, function(f)unlist(f$pl)-unlist(fit$pl)))),1,max)
  maxlist <- relist(maxabsdiff, fit$pl)
  ret <- as.data.frame(unlist(lapply(maxlist,function(x)if(length(x)>0)max(x) else NULL)))
  logLik <- max(abs(unlist(lapply(jj, function(f) f$opt$objective -fit$opt$objective))))
  ret <- rbind(ret,  logLik=logLik)
  sumGamma <- max(abs(unlist(lapply(jj, function(f) f$rl$sumGamma -fit$rl$sumGamma))))
  ret <- rbind(ret,  sumGamma = sumGamma)
  names(ret) <- "max(|delta|)"
  attributes(ret)$runs = jj
  attributes(ret)$fit = fit
  return(ret)
}

jj = jit(fit,nojit = 10,sd = 0.2)
jj

