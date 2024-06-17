#Save latent AR1 data
rm(list = ls())
library(RTMB)

set.seed(123)
n = 500
sd =0.5
rho = 0.9
gamma = rep(0,n)

gamma[1] = rnorm(1,0,sqrt(sd*sd/(1-rho^2)))
for(i in 2:length(gamma)){
  gamma[i] = rho*gamma[i-1] + rnorm(1,0,sd)
}
y = rpois(n,exp(gamma))
png("randomEffects/figures/latentAR1.png")
plot(y,xlim = c(1,100), width = 600)
lines(exp(gamma),col = 'green',lwd = 3)
dev.off()

saveRDS(y,file = "randomEffects/exercise/ar1Latent")
saveRDS(gamma,file = "randomEffects/exercise/ar1LatentGamma")
