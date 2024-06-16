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


plot(y,xlim = c(1,100))

#Store everyting we need in a list
fit = list(obj = obj,
           opt = opt,
           sdrep = sdrep,
           pl = pl,
           plsd = plsd,
           rl = rl,
           rlsd = rlsd,
           data = data)



#Check.consisticy
test = checkConsistency(fit$obj,n = 200)
summary(test)
