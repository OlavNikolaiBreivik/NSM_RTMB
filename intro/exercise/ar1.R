library(RTMB)
y = readRDS("intro/exercise/ar1.Rds")
plot(y,type = 'l',lwd = 2,cex.lab = 1.6, cex.axis = 1.5)

#Set up and estimate model
data = list(y = y)
par = list(logsd = 0,
           logitRho = 0)

f<-function(par){
  getAll(dat,par)
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
    #use RTMB::dautoreg
  }
  
  ADREPORT(rho) 
  nll
}

dat$code=0
#...
dat$code=1
#...

