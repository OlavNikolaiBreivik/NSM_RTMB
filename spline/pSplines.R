library(RTMB)
library(mgcv) #Use gam
library(Matrix) #Use sparse matrices

#Load the data
Vegetation <- read.table(file = "Vegetation.txt", header = TRUE, dec = ".")
Vegetation = Vegetation[!is.na(Vegetation$Richness),]

#Set up spline structure by using mgcv
gam_setup = gam(Richness ~ s(ROCK, bs = "cs") +
      s(LITTER, bs = "cs") + s(BARESOIL, bs = "cs") +
      s(FallPrec, bs = "cs") + s(SprTmax, bs = "cs"),
    data = Vegetation,fit=FALSE)

#Extrtact penelization matrices
S_ROCK = gam_setup$smooth[[1]]$S[[1]]
S_LITTER = gam_setup$smooth[[2]]$S[[1]]
S_BARESOIL = gam_setup$smooth[[3]]$S[[1]]
S_FallPrec = gam_setup$smooth[[4]]$S[[1]]
S_SprTmax = gam_setup$smooth[[5]]$S[[1]]

S_list = list(S_ROCK,S_LITTER,S_BARESOIL,S_FallPrec,S_SprTmax)
S_combined = .bdiag(S_list)         # join S's in sparse matrix
Sdims = unlist(lapply(S_list,nrow)) # Find dimension of each S


#For report, used for constructing plots
ROCK=seq(min(Vegetation$ROCK),max(Vegetation$ROCK),by = 1)
LITTER=seq(min(Vegetation$LITTER),max(Vegetation$LITTER),by = 1)
BARESOIL=seq(min(Vegetation$BARESOIL),max(Vegetation$BARESOIL),by = 1)
FallPrec=seq(min(Vegetation$FallPrec),max(Vegetation$FallPrec),by = 0.2)
SprTmax=seq(min(Vegetation$SprTmax),max(Vegetation$SprTmax),by = 0.2)

rockReport = PredictMat(gam_setup$smooth[[1]],data = data.frame(ROCK))
litterReport = PredictMat(gam_setup$smooth[[2]],data = data.frame(LITTER))
soilReport = PredictMat(gam_setup$smooth[[3]],data = data.frame(BARESOIL))
fallReport = PredictMat(gam_setup$smooth[[4]],data = data.frame(FallPrec))
sprReport = PredictMat(gam_setup$smooth[[5]],data = data.frame(SprTmax))

designMatrixForReport = list(rockReport,litterReport,soilReport,fallReport,sprReport)

#Define data object
data = list(Y = Vegetation$Richness, # Response
            X = gam_setup$X[,-1],  # Design matrix, without intercept
            S = S_combined,      # Combined penalty matrix
            Sdims = Sdims,
            designMatrixForReport = .bdiag(designMatrixForReport))

#Define parameter object
par = list(
  beta0 = 0,  # Intercept
  beta = rep(0,sum(Sdims)),  # Spline coefficients
  log_lambda = rep(rep(1,length(Sdims))), #Log spline penalization coefficients
  log_sigma = 0
)


#Define objective function
f = function(par){
  getAll(par,data)
    
  sigma = exp(log_sigma);
  lambda = exp(log_lambda);
  nll=0;
  k=1;  
  for(i in 1:length(Sdims)){
    m_i = Sdims[i]
    beta_i = beta[k:(k+m_i-1)]; #Recover betai
    S_i = S[k:(k+m_i-1),k:(k+m_i-1)]#Recover Si
    nll = nll- (0.5*m_i*log_lambda[i] - 0.5*lambda[i]*t(beta_i)%*% S_i  %*%beta_i);
    k = k+ m_i;
  }
  mu = beta0 + X%*%beta;
  nll = nll- sum(dnorm(Y, mu, sigma, TRUE));
  
  splineForReport = designMatrixForReport%*%beta;
  ADREPORT(splineForReport);
  ADREPORT(beta);
  return(nll)
}

#Fit model
obj = RTMB::MakeADFun(f, par,random="beta", silent = TRUE)
opt = nlminb(obj$par,obj$fn,obj$gr)
sdrep = sdreport(obj)

#Fit model with mgcv and compare splines
gam_setup = gam(Richness ~ s(ROCK, bs = "cs") +
                  s(LITTER, bs = "cs") + s(BARESOIL, bs = "cs") +
                  s(FallPrec, bs = "cs") + s(SprTmax, bs = "cs"),
                data = Vegetation,method = "ML",fit=TRUE) #Note that we use maximum likelihood method in the gam function.
par(mfrow=c(2,3))

#Plot splines
rl = as.list(sdrep,what = "Est",report = TRUE)
rlsd = as.list(sdrep,what = "Std",report = TRUE)
par(mfrow=c(2,3))
start = 1
stop = start + length(ROCK) -1
plot(gam_setup,select = 1,main = "Spline for ROCK")
lines(ROCK, rl$splineForReport[start:stop],col = 'red')
lines(ROCK, rl$splineForReport[start:stop]- 2*rlsd$splineForReport[start:stop], lty=2,col = 'red')
lines(ROCK, rl$splineForReport[start:stop]+ 2*rlsd$splineForReport[start:stop], lty=2,col = 'red')

start = stop +1
stop = start+ length(LITTER)-1
plot(gam_setup,select = 2,main = "Spline for Litter")
lines(LITTER, rl$splineForReport[start:stop], col = 'red')
lines(LITTER, rl$splineForReport[start:stop]- 2*rlsd$splineForReport[start:stop], lty=2,col = 'red')
lines(LITTER, rl$splineForReport[start:stop]+ 2*rlsd$splineForReport[start:stop], lty=2,col = 'red')

start = stop +1
stop = start+ length(BARESOIL)-1
plot(gam_setup,select = 3,main = "Spline for BARESOIL")
lines(BARESOIL, rl$splineForReport[start:stop], col = 'red')
lines(BARESOIL, rl$splineForReport[start:stop]- 2*rlsd$splineForReport[start:stop], col = 'red')
lines(BARESOIL, rl$splineForReport[start:stop]+ 2*rlsd$splineForReport[start:stop], col = 'red')

start = stop +1
stop = start+ length(FallPrec)-1
plot(gam_setup,select = 4,main = "Spline for FallPrec")
lines(FallPrec, rl$splineForReport[start:stop], col = 'red')
lines(FallPrec, rl$splineForReport[start:stop]- 2*rlsd$splineForReport[start:stop], col = 'red')
lines(FallPrec, rl$splineForReport[start:stop]+ 2*rlsd$splineForReport[start:stop], col = 'red')

start = stop +1
stop = start+ length(SprTmax)-1
plot(gam_setup,select = 5,main = "Spline for SprTmax")
lines(SprTmax, rl$splineForReport[start:stop], col = 'red')
lines(SprTmax, rl$splineForReport[start:stop]- 2*rlsd$splineForReport[start:stop], col = 'red')
lines(SprTmax, rl$splineForReport[start:stop]+ 2*rlsd$splineForReport[start:stop], col = 'red')



