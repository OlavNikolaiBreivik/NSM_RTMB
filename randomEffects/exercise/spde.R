library(RTMB)
dat_l = readRDS("randomEffects/exercise/haddockTiny.rds")
locUTM = cbind(dat_l$UTMX,dat_l$UTMY)

mesh <- fmesher::fm_mesh_2d(locUTM,
  max.edge = c(20 , 50)
)
plot(mesh)

spde <- fmesher::fm_fem(mesh)
spdeMatrices = list(c0 = spde$c0, g1 = spde$g1,g2 = spde$g2)
A <- fmesher::fm_basis(
  mesh,
  loc = locUTM)

dat <- list(catch = dat_l$catch,
             distance = dat_l$distance,
             spdeMatrices = spdeMatrices,
             A = A)

par <- list(beta= 0,
              logTau= 0,
              logKappa = -4.5,
              gamma = rep(0, mesh$n))


f <- function(par) {
  getAll(par,dat)
  kappa <- exp(logKappa)
  tau <- exp(logTau)

  #Precision matrix Q, e.g. Lindgren + Rue 2015 JSS p4:
  Q <- tau^2*(kappa^4 * spdeMatrices$c0 + 2 * kappa^2 * spdeMatrices$g1 + spdeMatrices$g2)
  nll = ...
  #Linear inerpolation
  delta = as.vector(A%*%gamma)
  
  nll = nll - ...
  return(nll)
}



#Plot spatial effect, read map
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

