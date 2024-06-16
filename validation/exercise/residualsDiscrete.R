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
  return(nll)
}


obj = MakeADFun(nll,par,random = "gamma")
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
sdrep = sdreport(obj,getJointPrecision = TRUE)

#OSA
res = RTMB::oneStepPredict(obj,
                           method="oneStepGeneric",
                           discrete=TRUE,
                           range=c(0,Inf))
acf(res$residual)
plot(res$residual)
qqnorm(res$residual)
abline(0,1)

map = list(logitRho = as.factor(NA))
obj = MakeADFun(nll,par,random= "gamma",map = map)
opt = nlminb(obj$par,obj$fn,obj$gr)
res = RTMB::oneStepPredict(obj,
                           method="oneStepGeneric",
                           discrete=TRUE,
                           range=c(0,Inf))
acf(res$residual)
plot(res$residual)
qqnorm(res$residual)
abline(0,1)
