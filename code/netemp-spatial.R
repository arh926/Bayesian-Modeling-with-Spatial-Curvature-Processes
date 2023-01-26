require(sp)
require(spBayes)
data("NETemp.dat")
require(rgdal)
require(spWombling)
# usa_shape <- raster::getData("GADM",country="USA", level=2)
# usam <- subset(usa_shape, NAME_1 != "Alaska" & NAME_1 != "Hawaii")
# usam.crs <- spTransform(usam,CRS("+proj=utm +zone=10 +datum=WGS84"))

sputm <- SpatialPoints(NETemp.dat[,c("UTMX","UTMY")], proj4string=CRS("+proj=utm +zone=10 +datum=WGS84")) 
spgeo <- spTransform(sputm, CRS("+proj=longlat +datum=WGS84"))
NETemp.dat$LON <- coordinates(spgeo)[,1]
NETemp.dat$LAT <- coordinates(spgeo)[,2]

NETemp.bdry <- chull(NETemp.dat$LON,NETemp.dat$LAT)
NETemp.dat[NETemp.bdry,c("LON","LAT")]
NETemp.shp <- raster::spPolygons(as.matrix(NETemp.dat[NETemp.bdry,c("LON","LAT")],nc=2), crs=sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))

# pdf("/Users/aritrahader/Dropbox/PhD Research/paper 2/plots/netemp-data.pdf", width=17, height=6)
mat <- matrix(c(1,2,3), nr=1,nc=3, byrow=T)
layout(mat,
       widths = c(5,5,2),
       heights = c(3,3))
hist(NETemp.dat$y.1, xlab=latex2exp::TeX("Temperature (in $\\degree$C)"), main="", col="lightblue",breaks=50)
sp_plot(col.seq.length = 11,
        col.text = "Spectral",
        data_frame = cbind(NETemp.dat$LON,
                           NETemp.dat$LAT,
                           NETemp.dat$y.1),
        xlab = latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        shape = NETemp.shp,
        contour.plot = T)
# dev.off()

coords <- cbind(NETemp.dat$LON,
                NETemp.dat$LAT)
N <- nrow(coords)
y <- NETemp.dat$y.1
X <- matrix(rep(1,N),nc=1)
D = 1:N
# cannot be randomly started
z_init= rep(0,N)
lower_phis=3/max(dist(coords))
upper_phis=300
lower_nu=1
upper_nu=3
shape_sigma=2
scale_sigma=1
shape_tau=2
scale_tau=1
mean_beta=rep(0,ncol(X))
prec_beta=1e-6*diag(ncol(X))
steps_init=1
stepnu_init=1
niter=1e4
nburn=niter/2
report=1e2
DtD <- diag(N)
cov.type <- "matern2"


# random starting points
phis_init = 0.2 # runif(1,lower_phis,upper_phis)
sigma2_init = 50 # 1/rgamma(1,shape=shape_sigma,rate=1/scale_sigma)
tau2_init = 1 # 1/rgamma(1,shape=shape_tau,rate=1/scale_tau)
beta_init = rep(0,ncol(X))
nu_init=1.1
mc_sp <- hlmBayes_sp(coords = coords,
                     y=y,
                     X=X,
                     z_init=z_init,
                     D = D,
                     phis_init=phis_init,
                     lower_phis=lower_phis,
                     upper_phis=upper_phis,
                     sigma2_init=sigma2_init,
                     shape_sigma=shape_sigma,
                     scale_sigma=scale_sigma,
                     tau2_init=tau2_init,
                     shape_tau=shape_tau,
                     scale_tau=scale_tau,
                     beta_init=beta_init,
                     mean_beta=mean_beta,
                     prec_beta=prec_beta,
                     steps_init=steps_init,
                     niter=niter,
                     report=report,
                     cov.type = cov.type,
                     nu_est = F,
                     nu_init = nu_init,
                     lower_nu = lower_nu,
                     upper_nu = upper_nu,
                     stepnu_init = stepnu_init)

model_summary <- hlm_summary(chain = mc_sp,nburn = nburn,thin = 1)
coef <- model_summary$summary.pars; round(coef,4)

z <- model_summary$summary.latent
z$sig <- apply(z,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})

# Checks: containment Probability
1-sum(apply(cbind(y,z[,"lower.hpd"]+coef["post_beta","lower.hpd"], z[,"upper.hpd"]+coef["post_beta","upper.hpd"]),
            1,
            function(x){
              if(x[1]>=x[2] & x[1]<=x[3]) return(0)
              else return(1)
            }))/N
y.hat <- z[,"median"]+coef["post_beta","median"]
y.upper.hpd <- z[order(y.hat),"upper.hpd"]+coef["post_beta","upper.hpd"]
y.lower.hpd <- z[order(y.hat),"lower.hpd"]+coef["post_beta","lower.hpd"]


# MSE
round(sqrt(mean((y-y.hat)^2)),2)

gp0 <- lm(y.1 ~ 1, data = NETemp.dat)
summary(gp0)
round(coef,3);round(coefficients(gp0),3)

mat <- matrix(c(1:8), nr=2,nc=4, byrow=T)
layout(mat,
       widths = c(rep(c(5,3),2)),
       heights = c(rep(3,4)))

sp_plot(11,"Spectral",cbind(coords,z[,1]),points.plot = T,contour.plot = T,sig = z$sig)
sp_plot(11,"Spectral",cbind(coords,y.hat),points.plot = T,contour.plot = T,sig = z$sig) 
sp_plot(11,"Spectral",cbind(coords,y),points.plot = T,contour.plot = T)
sp_plot(11,"Spectral",cbind(coords,y-y.hat),points.plot = T,contour.plot = T)

# legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
# plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
# text(x=2, y = seq(0.01,0.99,l=6), labels = round(seq(min(boston.c$CMEDV),max(boston.c$CMEDV),l=6),2),cex=2)
# rasterImage(legend_image, 0, 0, 1,1)

# RMSE
par(mfcol=c(1,1))
round(sqrt(mean((y-y.hat)^2)),4)
plot(y=y,x=y.hat,ylab="Observed",xlab="Estimated");abline(c(0,1))
plot(y-y.hat,ylab="Error");abline(h=0)

grid.points <- as.matrix(expand.grid(seq(min(NETemp.dat$LON),max(NETemp.dat$LON),by=0.5),
                                     seq(min(NETemp.dat$LAT),max(NETemp.dat$LAT),by=0.5)),nc=2)
tmp <- over(SpatialPoints(grid.points,proj4string = CRS(proj4string(NETemp.shp))),NETemp.shp)
grid.points <- grid.points[!is.na(tmp),]

# sp_plot(11,"Spectral",cbind(coords,y))
# points(grid.points, cex=0.2)

samples <- (nburn+1):niter
gradient_est <- spatial_gradient(coords=coords,
                                 cahin = mc_sp,
                                 grid.points = grid.points,
                                 samples = samples,
                                 cov.type = cov.type,
                                 ncores=5,
                                 nbatch = 100,
                                 nburn = nburn,
                                 niter = niter)
grad.s1.hpd <- data.frame(gradient_est$grad1.est)
grad.s1.hpd$sig <- apply(grad.s1.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

grad.s2.hpd <- data.frame(gradient_est$grad2.est)
grad.s2.hpd$sig <- apply(grad.s2.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

grad.s11.hpd <- data.frame(gradient_est$grad11.est)
grad.s11.hpd$sig <- apply(grad.s11.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  
grad.s12.hpd <- data.frame(gradient_est$grad12.est)
grad.s12.hpd$sig <- apply(grad.s12.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  
grad.s22.hpd <- data.frame(gradient_est$grad22.est)
grad.s22.hpd$sig <- apply(grad.s22.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

cex.pt <- 1
cex.txt <- 1.5
pdf("netemp-spatial.pdf",width=19.5,height=7)
mat <- matrix(c(1:12), nr=2,nc=6, byrow=T)
layout(mat,
       widths = c(rep(c(5,1.5),3)),
       heights = c(rep(3,6)))
sp_plot(11,"Spectral",shape=NETemp.shp,cbind(coords,z[,1]), points.plot = T,contour.plot = T,legend=T,sig = z$sig)
sp_plot(11,"Spectral",shape=NETemp.shp,cbind(grid.points,grad.s1.hpd[,1]),  points.plot = T,contour.plot = T,legend=T,sig = grad.s1.hpd$sig)
sp_plot(11,"Spectral",shape=NETemp.shp,cbind(grid.points,grad.s2.hpd[,1]),  points.plot = T,contour.plot = T,legend=T,sig = grad.s2.hpd$sig)
sp_plot(11,"Spectral",shape=NETemp.shp,cbind(grid.points,grad.s11.hpd[,1]),  points.plot = T,contour.plot = T,legend=T,sig = grad.s11.hpd$sig)
sp_plot(11,"Spectral",shape=NETemp.shp,cbind(grid.points,grad.s12.hpd[,1]),  points.plot = T,contour.plot = T,legend=T,sig = grad.s12.hpd$sig)
sp_plot(11,"Spectral",shape=NETemp.shp,cbind(grid.points,grad.s22.hpd[,1]),  points.plot = T,contour.plot = T,legend=T,sig = grad.s22.hpd$sig)
dev.off()
# Posterior Surface Analysis
det_est <- eigen_est_1 <- eigen_est_2 <- laplace_est <- div_est <- theta1 <- theta2 <- matrix(NA,nr=5000,nc=nrow(grid.points))
for(i in 1:5000){
  for(j in 1:nrow(grid.points)){
    hess.mat <- matrix(c(gradient_est$grad.s11.mcmc[i,j],
                         gradient_est$grad.s12.mcmc[i,j],
                         gradient_est$grad.s12.mcmc[i,j],
                         gradient_est$grad.s22.mcmc[i,j]),nr=2,nc=2,byrow = T)
    det_est[i,j] <- det(hess.mat)
    eigen.hess.mat <- eigen(hess.mat)$values
    eigen_est_1[i,j] <- eigen.hess.mat[1]; eigen_est_2[i,j] <- eigen.hess.mat[2]
    laplace_est[i,j] <- sum(diag(hess.mat))
    div_est[i,j] <- gradient_est$grad.s1.mcmc[i,j]+gradient_est$grad.s2.mcmc[i,j]
    # principal curvature
    t1 <- gradient_est$grad.s1.mcmc[i,j]
    t2 <- gradient_est$grad.s2.mcmc[i,j]
    theta1[i,j] <- atan(t2/t1)
    #if(t1<0 & t2>0) theta1[i,j] <- pi/2+theta1[i,j]
    #else if(t1<0 & t2<0) theta1[i,j] <- pi+theta1[i,j]
    #else if(t1>0 & t2<0) theta1[i,j] <- 3*pi/2+theta1[i,j]
    t1 <- gradient_est$grad.s12.mcmc[i,j]
    t2 <- gradient_est$grad.s22.mcmc[i,j]-gradient_est$grad.s11.mcmc[i,j]
    theta2[i,j] <- 0.5*atan(2*t1/t2)
    #if(t1<0 & t2>0) theta2[i,j] <- pi/2+theta2[i,j]
    #else if(t1<0 & t2<0) theta2[i,j] <- pi+theta2[i,j]
    #else if(t1>0 & t2<0) theta2[i,j] <- 3*pi/2+theta2[i,j]
  }
}
slist <- split(1:5000, ceiling(seq_along(1:5000)/(length(1:5000)/100)))
det_est_val <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(det_est[x,],2,median))),2,median),
                                HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(det_est[x,],2,median))))))
det_est_val$sig <- apply(det_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
eigen_est_val1 <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(eigen_est_1[x,],2,median))),2,median),
                                   HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(eigen_est_1[x,],2,median))))))
eigen_est_val1$sig <- apply(eigen_est_val1,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
eigen_est_val2 <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(eigen_est_2[x,],2,median))),2,median),
                                   HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(eigen_est_2[x,],2,median))))))
eigen_est_val2$sig <- apply(eigen_est_val2,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
laplace_est_val <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(laplace_est[x,],2,median))),2,median),
                                    HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(laplace_est[x,],2,median))))))
laplace_est_val$sig <- apply(laplace_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
div_est_val <- cbind.data.frame(apply(do.call(rbind,lapply(slist, function(x) apply(div_est[x,],2,median))),2,median),
                                HPDinterval(as.mcmc(do.call(rbind,lapply(slist, function(x) apply(div_est[x,],2,median))))))
div_est_val$sig <- apply(div_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
# theta1_est <- apply(REdaS::rad2deg(theta1),2,median)
# theta2_est <- apply(REdaS::rad2deg(theta2),2,median)
# par(mfcol=c(1,1))
# sp_plot(11,"Spectral",cbind(coords,z[,1]),
#         shape=NETemp.shp,
#         contour.plot = T,
#         xlab=latex2exp::TeX("Longitude$\\degree$"),
#         ylab = latex2exp::TeX("Latitude$\\degree$"),
#         points.plot=F,
#         raster.surf = F,
#         sig = z$sig,
#         legend = F)
# # points(grid.points, col="red", pch=16)
# for(i in 1:nrow(grid.points)){
#   arrows(x0 = grid.points[i,1],
#          y0 = grid.points[i,2],
#          x1 = grid.points[i,1]+0.2*cos(theta1_est[i]),
#          y1 =grid.points[i,2]+0.2*sin(theta1_est[i]),
#          length = 0.05)
# }
# sp_plot(11,"Spectral",cbind(coords,z[,1]),
#         shape=NETemp.shp,
#         contour.plot = T,
#         xlab=latex2exp::TeX("Longitude$\\degree$"),
#         ylab = latex2exp::TeX("Latitude$\\degree$"),
#         points.plot=F,
#         raster.surf = F,
#         sig = z$sig,
#         legend = F)
# abline(v=grid.points[,1])
# abline(h=grid.points[,2])
# # points(grid.points, col="red", pch=16)
# for(i in 1:nrow(grid.points)){
#   arrows(x0 = grid.points[i,1],
#          y0 = grid.points[i,2],
#          x1 = grid.points[i,1]+0.2*cos(theta2_est[i]),
#          y1 =grid.points[i,2]+0.2*sin(theta2_est[i]),
#          length = 0.05)
# }
mat <- matrix(c(1:6), nr=1,nc=6, byrow=T)
layout(mat,
       widths = c(rep(c(5,3),3)),
       heights = c(rep(3,6)))
# rster.obj <- sp_plot(11,"Spectral",cbind(coords,z[,1]),
#                      shape=NETemp.shp,
#                      contour.plot = T,
#                      xlab=latex2exp::TeX("Longitude$\\degree$"),
#                      ylab = latex2exp::TeX("Latitude$\\degree$"),
#                      points.plot=T,
#                      raster.surf = T,
#                      sig = z$sig)

# sp_plot(11,"Spectral",cbind(grid.points/100,grad.s1.hpd[,1]),
#         shape = boston.shp,contour.plot = T,
#         xlab=latex2exp::TeX("Longitude$\\degree$"),
#         ylab = latex2exp::TeX("Latitude$\\degree$"),
#         points.plot = T, sig=grad.s1.hpd$signif)
# 
# 
# sp_plot(11,"Spectral",cbind(grid.points/100,grad.s2.hpd[,1]), shape=boston.shp,contour.plot = T,
#         xlab=latex2exp::TeX("Longitude$\\degree$"),
#         ylab = latex2exp::TeX("Latitude$\\degree$"),
#         points.plot = T, sig=grad.s2.hpd$signif)
# 
# sp_plot(11,"Spectral",cbind(grid.points/100,grad.s11.hpd[,1]),contour.plot = T, shape=boston.shp,
#         xlab=latex2exp::TeX("Longitude$\\degree$"),
#         ylab = latex2exp::TeX("Latitude$\\degree$"),
#         points.plot = T, sig=grad.s11.hpd$signif)
# 
# sp_plot(11,"Spectral",cbind(grid.points/100,grad.s12.hpd[,1]),contour.plot = T, shape=boston.shp,
#         xlab=latex2exp::TeX("Longitude$\\degree$"),
#         ylab = latex2exp::TeX("Latitude$\\degree$"),
#         points.plot = T, sig=grad.s12.hpd$signif)
# 
# 
# sp_plot(11,"Spectral",cbind(grid.points/100,grad.s22.hpd[,1]),contour.plot = T, shape=boston.shp,
#         xlab=latex2exp::TeX("Longitude$\\degree$"),
#         ylab = latex2exp::TeX("Latitude$\\degree$"),
#         points.plot = T, sig=grad.s22.hpd$signif)

# mat <- matrix(c(1:4), nr=1,nc=4, byrow=T)
# layout(mat,
#        widths = c(rep(c(5,3),2)),
#        heights = c(rep(3,4)))

# sp_plot(11,"Spectral",cbind(grid.points/100,eigen_est_val1[,1]),
#         contour.plot = T, shape=boston.shp,
#         xlab=latex2exp::TeX("(a) $\\widehat{\\lambda_1}$"),
#         points.plot=T,sig=eigen_est_val1$sig)
# 
# 
# sp_plot(11,"Spectral",cbind(grid.points/100,eigen_est_val2[,1]),
#         contour.plot = T, shape=boston.shp,
#         xlab=latex2exp::TeX("(b) $\\widehat{\\lambda_2}$"),
#         points.plot=T,sig=eigen_est_val2$sig)
# 
sp_plot(11,"Spectral",cbind(grid.points,det_est_val[,1]),
        contour.plot = T, shape=NETemp.shp,
        xlab=latex2exp::TeX("(c) $\\widehat{|\\nabla^2Y|}$"),
        points.plot=T,sig=det_est_val$sig)

sp_plot(11,"Spectral",cbind(grid.points,div_est_val[,1]),
        contour.plot = T, shape=NETemp.shp,
        xlab=latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        points.plot=T,sig=div_est_val$sig,
        grid=F)

sp_plot(11,"Spectral",cbind(grid.points,laplace_est_val[,1]),
        contour.plot = T, shape=NETemp.shp,
        xlab=latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        points.plot=T,sig=laplace_est_val$sig,
        grid=F)

# sp_plot(11,"Spectral",cbind(coords/100,y.hat),
#         contour.plot = T, shape=boston.shp,
#         xlab=latex2exp::TeX("(f) $\\widehat{y}$"),
#         points.plot=T,sig=z$sig)

par(mfcol=c(1,1))
rster.obj <- sp_plot(11,"Spectral",cbind(coords,z[,1]),
                     shape=NETemp.shp,
                     contour.plot = T,
                     xlab=latex2exp::TeX("Longitude$\\degree$"),
                     ylab = latex2exp::TeX("Latitude$\\degree$"),
                     points.plot=T,
                     raster.surf = T,
                     sig = z$sig,
                     legend = F)
# Wombling
x <- raster::rasterToContour(rster.obj)#,nlevels=20
# Test for wombling
x.levels <- as.numeric(as.character(x$level))
level.choice <- c(-6,6) #c(-2,2) #c(-10,10)
#c(x.levels[2],x.levels[length(x.levels)-1])
subset.points.1 <- subset(x,level==level.choice[1])
subset.points.2 <- subset(x,level==level.choice[2])
#1,2, 5
lines(subset.points.1@lines[[1]]@Lines[[1]]@coords,lwd=3.5)
lines(subset.points.1@lines[[1]]@Lines[[2]]@coords,lwd=3.5)
lines(subset.points.2@lines[[1]]@Lines[[1]]@coords,lwd=3.5)
lines(subset.points.2@lines[[1]]@Lines[[2]]@coords,lwd=3.5)
subset.points <- list(subset.points.1@lines[[1]]@Lines[[1]]@coords,
                      subset.points.1@lines[[1]]@Lines[[2]]@coords,
                      subset.points.2@lines[[1]]@Lines[[1]]@coords,
                      subset.points.2@lines[[1]]@Lines[[2]]@coords)

# Use cwomb-riemann or bayes_cwomb to compute wombling measures on each curve