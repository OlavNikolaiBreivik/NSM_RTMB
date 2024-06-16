
#Simulate data#########
set.seed(123)
n = 500
sd =1
rho = 0.98
gamma = rep(0,n)

gamma[1] =0# rnorm(1,0,sqrt(sd*sd/(1-rho^2))) 
for(i in 2:length(gamma)){
  gamma[i] = rho*gamma[i-1] + rnorm(1,0,sd)
}

alpha = 20
beta = 1.15

y = exp(rnorm(n, gamma, sd = sqrt(log(alpha* exp(gamma)^(beta-2) + 1))))

plot(gamma)

plot(log(y))

logObs = log(y)
saveRDS(logObs,file = "logObs.RData")
