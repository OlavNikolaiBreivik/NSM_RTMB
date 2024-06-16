library(RTMB)
rm(list = ls())

par = list(rho = 0.9)
data  = list(n = 10)

f = function(par){
  "[<-" = ADoverload("[<-")
  getAll(par,data)
  
  eta = rep(100,n)
  x = rep(rho,n)
  for(i in 2:n){
    #    #Do nothing
  }
  #  eta[1] = eta[1] - x[1]
  eta[1] = eta[1] - x[1]
  
  return(sum(eta))
}

f(par)
obj = MakeADFun(f,par)
