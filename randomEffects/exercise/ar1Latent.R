library(RTMB)
y = readRDS("randomEffects/exercise/ar1Latent.rds")

#Set up and estimate model
dat = list(y = y)
par = list(logsd = 0,
           logitRho = 0,
           gamma = rep(0,length(y)))


nll = function(par){
  getAll(par,dat)
  sd = exp(logsd)
  rho = 2/(1 + exp(-logitRho))-1

  #TODO, implement the likelihood 
  
  return(nll)
}


