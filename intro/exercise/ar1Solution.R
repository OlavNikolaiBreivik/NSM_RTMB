library(RTMB)
y = readRDS("intro/exercise/ar1.Rds")

#Set up and estimate model
data = list(y = y)
par = list(logsd = 0,
           logitRho = 0)

f<-function(par){
  getAll(data,par)
  rho = 2/(1 + exp(-logitRho))-1
  timeSteps <-length(y)
  sd <- exp(logsd)
  
  nll<-0
  
  if(code==0){
    nll <- nll -dnorm(y[1],0,sqrt(sd*sd/(1-rho*rho)),log=TRUE)
    for(i in 2:timeSteps){    
      nll <- nll -dnorm(y[i],rho*y[i-1],sd,log=TRUE)
    }
  }
  
  if(code==1){
    nll <- nll - dautoreg(y,phi=rho,scale=sqrt(sd*sd/(1-rho*rho)), log=TRUE)
  }
  
  if(code==2){
    t<-1:length(y)  
    S <- sd^2/(1-rho^2)*rho^abs(outer(t,t,"-"))  
    nll <- nll - dmvnorm(y,0,S,log=TRUE)
  }
  
  if(code==3){
    Q<-Matrix::spMatrix(timeSteps,timeSteps)
    diag(Q)<-c(1,rep(1+rho^2, timeSteps-2),1)
    Q[cbind(2:timeSteps,1:(timeSteps-1))] <- -rho
    Q[cbind(1:(timeSteps-1), 2:timeSteps)] <- -rho
    #Q <- Q/sd/sd or scale
    nll <- nll - dgmrf(y, Q=Q, scale=sd, log=TRUE)
  }
  ADREPORT(rho) 
  nll
}
data$code=0
obj <- MakeADFun(f,par, silent=TRUE)
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
c(opt$obj, opt$par)
data$code=1
obj <- MakeADFun(f,par, silent=TRUE)
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
c(opt$obj, opt$par)
data$code=2;par$logitRho = 0.01
obj <- MakeADFun(f,par, silent=TRUE)
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
c(opt$obj, opt$par)
data$code=3
obj <- MakeADFun(f,par, silent=TRUE)
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
c(opt$obj, opt$par)


sdrep = sdreport(obj)
pl = as.list(sdrep,what = "Est")
plsd = as.list(sdrep,what = "Std")

# 95%CI intervals
c(exp(pl$logsd),exp(pl$logsd + 2*plsd$logsd*c(-1,1)))
c(2/(1 + exp(-pl$logitRho))-1,
  2/(1 + exp(-pl$logitRho + 2*plsd$logitRho*c(1,-1)))-1)

par$logitRho = 10
map = list(logitRho = as.factor(NA))
obj2 = MakeADFun(f,par,map = map)
opt2 = nlminb(obj2$par,obj2$fn,obj2$gr)

LR = -2*(opt$objective -opt2$objective)
1-pchisq(LR,df = 1)


