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
rl = as.list(sdrep,what = "Est",report = TRUE)
rlsd = as.list(sdrep,what = "Std",report = TRUE)

#Store everyting we need in a list
fit = list(obj = obj,
           opt = opt,
           sdrep = sdrep,
           par = par,
           pl = pl,
           plsd = plsd,
           rl = rl,
           rlsd = rlsd,
           data = data)


#Simulation#########################
nsim = 20
simStudy = list()
for(i in 1:nsim){
  data = list(y = fit$obj$simulate()$y)
  objSim = MakeADFun(nll,par,random = "gamma")
  optSim = nlminb(objSim$par,objSim$fn,objSim$gr)
  sdrepSim = sdreport(objSim)
  plSim = as.list(sdrepSim, what = "Est")
  simStudy[[i]] = list(obj = objSim,
                       opt = optSim,
                       sdrep = sdrepSim,
                       pl = plSim,
                       data = data)
}


library(plotrix)
parEst = c(fit$pl$logsd,fit$pl$logitRho)
parEstU = c(fit$pl$logsd + 2*fit$plsd$logsd,fit$pl$logitRho+ 2*fit$plsd$logitRho)
parEstL = c(fit$pl$logsd - 2*fit$plsd$logsd,fit$pl$logitRho- 2*fit$plsd$logitRho)
plotCI(parEst,li = parEstL,ui = parEstU,lwd = 2,
       ylab = "Estimate", xlab = "Parameter",xaxt='n',
       ylim = c(min(parEstL)-1,max(parEstU) + 1)
)
axis(1, at = 1:2,labels = c(expression(log~sigma),expression(logit~rho)))

for(i in 1:length(simStudy)){
  parEstSim = c(simStudy[[i]]$pl$logsd,simStudy[[i]]$pl$logitRho)
  points(parEstSim,col = 'red')
}

plotCI(parEst,li = parEstL,ui = parEstU,add = TRUE,lwd = 2)

