library(RTMB)
logObs = readRDS("validation/exercise/logObs.RData")
data = list(logObs = logObs)
par = list(logsdRho = 0,
           logitRho = 0,
           logsdObs = 0,
           gamma = rep(0,length(logObs)))


nll = function(par){
  getAll(par,data)
  sd = exp(logsdRho)
  rho = 2/(1 + exp(-logitRho))-1
  sdObs = exp(logsdObs)
  beta = sdObs
    
  logObs = OBS(logObs)
  
  nll = 0
  nll = nll -dnorm(gamma[1],0,sqrt(sd*sd/(1-rho^2)),TRUE);
  for(i in 2:length(gamma)){
    nll = nll -dnorm(gamma[i],rho*gamma[i-1],sd,TRUE);
  }
  
  nll = nll -sum(dnorm(logObs,gamma,sdObs,TRUE))
  return(nll)
}

obj = MakeADFun(nll,par,random = "gamma")
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
sdrep = sdreport(obj,getJointPrecision = TRUE)
pl = as.list(sdrep,what = "Est")
plsd = as.list(sdrep,what = "Std")


#Calculate residuals

