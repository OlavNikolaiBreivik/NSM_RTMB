library(RTMB)
dat = list()
dat$Y = c(13, 5, 28, 28, 15, 4, 13, 4, 10, 17, 11, 13, 12, 17, 3)
par = list()
par$logsize = 0
par$logitp = 0

f = function(par){
  size= exp(par$logsize)
  phi= plogis(par$logitp)
  nll= -sum(dnbinom(dat$Y,size,phi,log=TRUE))
  ADREPORT(size)
  return(nll)
}

obj = ...
opt = ...
sdrep = ...