library(RTMB)
#Read data
y = readRDS("randomEffects/exercise/ar1Latent.rds")
plot(y)

missing =numeric(0)#60:90 
dat = list(y = y,procedure = 1,missing = missing)
par = list(logsd = 0,
           logitRho = 0,
           gamma = rep(0,length(y)))


nll = function(par){
  getAll(par,dat)
  sd = exp(logsd)
  rho = 2/(1 + exp(-logitRho))-1

  y = OBS(y)

  nll = 0
  if(procedure ==1){
    nll = nll -dnorm(gamma[1],0,sqrt(sd*sd/(1-rho^2)),TRUE);
    for(i in 2:length(gamma)){
      nll = nll -dnorm(gamma[i],rho*gamma[i-1],sd,TRUE);
    }
  }else if(procedure==2){
    nll = -dautoreg(gamma, mu=0, phi=rho, scale=sqrt(sd*sd/(1-rho^2)), log=TRUE)
  }

  if(length(missing)>0){
    observed = seq(1,length(y))[-missing]
  }else{
    observed = seq(1,length(y))
  }
  nll = nll -sum(dpois(y[observed],exp(gamma[observed]),TRUE))

  sumExpGamma = sum(exp(gamma))
  ADREPORT(sumExpGamma)

  return(nll)
}


obj = MakeADFun(nll,par,random = "gamma")
opt <- nlminb(obj$par,obj$fn,obj$gr, control = list(trace = 1))
sdrep = sdreport(obj,getJointPrecision = TRUE)
#sdrep = sdreport(obj,bias.correct = TRUE)
pl = as.list(sdrep,what = "Est")
plsd = as.list(sdrep,what = "Std")
rl = as.list(sdrep,what = "Est",report = TRUE)
rlsd = as.list(sdrep,what = "Std",report = TRUE)

rl$sumExpGamma[1] + 2*rlsd$sumExpGamma[1]*c(-1,1)


gammaTruth = readRDS("randomEffects/exercise/ar1LatentGamma.rds")
sum(exp(gammaTruth))

plot(y,xlim = c(0,100))
lines(exp(gammaTruth),lwd = 3,col = 'green')

lines(exp(pl$gamma),col = 'red',lwd = 3)
lines(exp(pl$gamma + 2*plsd$gamma),lty = 2,col = 'red')
lines(exp(pl$gamma - 2*plsd$gamma),lty = 2,col = 'red')







##########################################
#Not with markov. Very time consuming for high dimension of gamma.
data = list(y = y)
par = list(logsd = 0,
           logitRho = 0,
           epsilon = rep(0,length(y)))

nllBad = function(par){
  getAll(par,data)

  sd = exp(logsd)
  rho = 2/(1 + exp(-logitRho))-1
  gamma = rep(0,length(epsilon))

  nll = -dnorm(epsilon[1],0,sqrt(sd*sd/(1-rho^2)),TRUE);
  gamma[1] = epsilon[1]
  for(i in 2:length(gamma)){
    gamma[i] = rho*gamma[i-1] + epsilon[i]
    nll = nll -dnorm(epsilon[i],0,sd,TRUE);
  }
  nll = nll -sum(dpois(y,exp(gamma),TRUE))
  return(nll)
}
objBad = MakeADFun(nllBad,par,random = "epsilon")
optBad = nlminb(objBad$par,objBad$fn,objBad$gr)

#Illustrate that we get exactly the same results
opt$objective - optBad$objective
opt$par-optBad$par

#look at sparsness
library(SparseM)
sdrepBad = sdreport(objBad,getJointPrecision = TRUE)
nameIndex = which(colnames(sdrepBad$jointPrecision)=="epsilon")
Q = sdrepBad$jointPrecision[nameIndex,nameIndex]
png("figures/sparsityResiduals.png")
image(Q[1:100,1:100], main = "Sparsity increments")
dev.off()


