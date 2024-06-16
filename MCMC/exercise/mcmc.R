library(RTMB)
y = readRDS("MCMC/exercise/ar1Latent.rds")

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


#mcmc
obj = MakeADFun(nll,par,random = "gamma")
library(tmbstan)
mcmc = tmbstan(obj, chains=1,iter = 10000)


#Run mle
obj = MakeADFun(nll,par,random = "gamma")
system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
sdrep = sdreport(obj,getJointPrecision = TRUE)
sdr <- summary(sdrep)
pl = as.list(sdrep,what = "Est")
plsd = as.list(sdrep,what = "Std")


mc <- extract(mcmc, pars=names(obj$par),
              inc_warmup=FALSE, permuted=FALSE)

pdf("figures/mcmcAll.pdf")
npar<-dim(mc)[3]
layout(rbind(rep(1,npar),c(2:(npar+1))))
matplot(mc[,1,], type="l")

post = as.matrix(mcmc)
plot(density(post[,'logsd']), lwd = 2)
x <- seq(min(post[,'logsd']), max(post[,'logsd']), length=50)
y <- dnorm(x,mean = pl$logsd, plsd$logsd)
lines(x,y, type = "l", lwd = 2, col = 'red')

plot(density(post[,'logitRho']), lwd = 2)
x <- seq(min(post[,'logitRho']), max(post[,'logitRho']), length=50)
y <- dnorm(x,mean = pl$logitRho, plsd$logitRho)
lines(x,y, type = "l", lwd = 2, col = 'red')
dev.off()

#shinystan
library(shinystan)
launch_shinystan(mcmc)



png("figures/mcmc.png",width = 700,heigh = 300)
matplot(mc[,1,], type="l")
dev.off()


