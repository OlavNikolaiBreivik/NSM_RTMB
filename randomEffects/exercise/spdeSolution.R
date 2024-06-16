library(RTMB)
dat_l = readRDS("randomEffects/exercise/haddockTiny.rds")

locUTM = cbind(dat_l$UTMX,dat_l$UTMY)
#Define mesh
mesh <- fmesher::fm_mesh_2d(locUTM,
  max.edge = c(30 , 70)
)

spde <- fmesher::fm_fem(mesh)
spdeMatrices = list(c0 = spde$c0, g1 = spde$g1,g2 = spde$g2)
A <- fmesher::fm_basis(
  mesh,
  loc = locUTM)


predUTM = readRDS("randomEffects/exercise/predPoints.rds")

A_pred <- fmesher::fm_basis(
  mesh,
  loc = predUTM)

dat <- list(catch = dat_l$catch,
             distance = dat_l$distance,
             spdeMatrices = spdeMatrices,
             A = A,
             A_pred = A_pred)

par <- list(beta= 0,
              logTau= 0,
              logKappa = -4.5,
              gamma = rep(0, mesh$n))


f <- function(par) {
  getAll(par,dat)
  kappa <- exp(logKappa)
  tau <- exp(logTau)

  catch = OBS(catch)

  #Precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  Q <- tau^2*(kappa^4 * spdeMatrices$c0 + 2 * kappa^2 * spdeMatrices$g1 + spdeMatrices$g2)
  ## GMRF prior
  nll <- -dgmrf(gamma, 0, Q, log=TRUE)
  ## data likelihood
  delta = as.vector(A%*%gamma)
  nll = nll - sum(dpois(catch,exp(beta +delta)*distance, TRUE))

  delta_pred = as.vector(A_pred%*%gamma)
  mu_pred = delta_pred
#  logAbundanceIndex = log(sum(exp(mu_pred))/length(mu_pred))
  logAbundanceIndex = log(mean(exp(mu_pred)))
  ADREPORT(logAbundanceIndex)

  range <- sqrt(8)/kappa
  ADREPORT(range)
  sigma <- 1/sqrt(4*pi*exp(2*logTau + 2*logKappa))
  ADREPORT(sigma)
  return(nll)
}

obj <- RTMB::MakeADFun(f, par, random = "gamma")
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 1))
rep = RTMB::sdreport(obj)
pl = as.list(rep,what = "Est")#Parameter estimates
plsd = as.list(rep,what = "Std") #Parameter sd

rl = as.list(rep,what = "Est",report = TRUE)
rlsd = as.list(rep,what = "Std",report = TRUE)

sim = obj$simulate()
#Plot spatial effect, read map
# Interpolate z values onto a regular grid
png("randomEffects/figures/spatialCatchSolution.png")
interp_surface <- akima::interp(mesh$loc[,1], mesh$loc[,2], pl$gamma,nx = 200,ny = 200)
fields::image.plot(interp_surface,
           col=viridis::turbo(500),
           xlab = '', ylab = '',
           xlim = c(0,1400),ylim = c(7550,8600),
           main = "",useRaster=T, horizontal=T
)
title(xlab = "Easting (km)", line = 2.5, cex.lab=1.8, xpd=NA)
title(ylab = "Northing (km)", line = 2.5, cex.lab=1.8, xpd=NA)
title(main = "Spatial MAP", line = 1, cex.main=2.5, xpd=NA)
points(locUTM,cex = log(dat$catch)/2)
points(locUTM[which(dat$catch==0),],pch = 4)
newmap <- maps::map("world", c("Norway","Russia","Iceland", "Greenland","Sweden","Finland"),fill = TRUE,plot = FALSE, col = "transparent")
IDs <- sapply(strsplit(newmap$names, ":"), function(x) x[1])
map.sp <- maptools::map2SpatialPolygons(
  newmap, IDs = IDs,
  proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
map.sp = sp::spTransform(map.sp, sp::CRS("+proj=utm +zone=35 ellps=WGS84 +units=km"))
library(sp)
plot(map.sp, add=T, col='grey')
dev.off()

png("randomEffects/figures/spatialCatchSolution.png")
interp_surface <- akima::interp(mesh$loc[,1], mesh$loc[,2], plsd$gamma,nx = 200,ny = 200)
fields::image.plot(interp_surface,
                   col=viridis::turbo(500),
                   xlab = '', ylab = '',
                   xlim = c(0,1400),ylim = c(7550,8600),zlim = range(plsd$gamma),
                   main = "",useRaster=T, horizontal=T
)
title(xlab = "Easting (km)", line = 2.5, cex.lab=1.8, xpd=NA)
title(ylab = "Northing (km)", line = 2.5, cex.lab=1.8, xpd=NA)
title(main = "Spatial MAP", line = 1, cex.main=2.5, xpd=NA)
points(locUTM,cex = log(dat$catch)/2)
points(locUTM[which(dat$catch==0),],pch = 4)
newmap <- maps::map("world", c("Norway","Russia","Iceland", "Greenland","Sweden","Finland"),fill = TRUE,plot = FALSE, col = "transparent")
IDs <- sapply(strsplit(newmap$names, ":"), function(x) x[1])
map.sp <- maptools::map2SpatialPolygons(
  newmap, IDs = IDs,
  proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
map.sp = sp::spTransform(map.sp, sp::CRS("+proj=utm +zone=35 ellps=WGS84 +units=km"))
library(sp)
plot(map.sp, add=T, col='grey')
dev.off()




#Simulate#######
obj$simulate() #Simulate latent effects given estimated model parameters

