## Example
# Data
set.seed(12345)    
F1 <- gl(2, 100)
F2 <- gl(2, 1, 200)
x <- runif(200, 0, 20)

# Parameters
a <- c(0.1, -0.1)
b <- c(1, 2, 3, 4)
s <- 0.1

# Simulated observations
y <- a[as.integer(F1)]*x + b[as.integer(F1:F2)] + rnorm(200, sd=s)

# Model runs
logLik(lm(y ~ -1 + F1:F2 + F1:x)) 
logLik(lm(y ~ -1 + F1:x + F1:F2))

#Do it manually
nll <- function(theta){
  pred <- theta[1:2][as.integer(F1)]*x + theta[3:6][as.integer(F1:F2)]
  -sum(dnorm(y, pred, exp(theta[7]), log=TRUE))
}
opt <- nlminb(numeric(7), nll)
cat("\'Manual log Lik.\' ", -opt$objective, "\n")

##





















# Model runs###########################################
r1 <- logLik(lm(y ~ -1 + F1:F2 + F1:x))  # Expected
r2 <- logLik(lm(y ~ -1 + F1:x + F1:F2))  # Unexpected - not the same


# Slight variations returns what I originally expected  
r3 <- logLik(lm(y ~ -1 + I(F1):x + F1:F2))  
r4 <- logLik(lm(y ~ -1 + F1:x + I(F1:F2)))

r1
# 'log Lik.' 184.5824 (df=7)
r2
# 'log Lik.' -263.9513 (df=5)
r3
# 'log Lik.' 184.5824 (df=7)
r4
# 'log Lik.' 184.5824 (df=7)

# We can also do it manually just to check:

nll <- function(theta){
  pred <- theta[1:2][as.integer(F1)]*x + theta[3:6][as.integer(F1:F2)]
  -sum(dnorm(y, pred, exp(theta[7]), log=TRUE))
}
opt <- nlminb(numeric(7), nll)
cat("\'Manual log Lik.\' ", -opt$objective, "\n")
# 'Manual log Lik.'  184.5824
