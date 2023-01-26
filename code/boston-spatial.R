# setwd("~/Dropbox/Desktop-Files/spatiotemporal_gradient")
source('sp_plot-1.R')
require(spData)
require(spWombling)
data(boston)

boston.bdry <- chull(boston.c$LON,boston.c$LAT)
boston.c[boston.bdry,c("LON","LAT")]
boston.shp <- raster::spPolygons(as.matrix(boston.c[boston.bdry,c("LON","LAT")],nc=2), crs=sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))

# pdf("/Users/aritrahader/Dropbox/PhD Research/paper 2/plots/boston-data.pdf", width=17, height=6)
mat <- matrix(c(1,2,3), nr=1,nc=3, byrow=T)
layout(mat,
       widths = c(5,5,2),
       heights = c(3,3))
hist(boston.c$CMEDV, xlab="Median House Prices (in USD 1000)", main="", col="lightblue",breaks=50)
sp_plot(11,"Spectral",data_frame = cbind(boston.c$LON,
                                         boston.c$LAT,
                                         boston.c$CMEDV),
        xlab=latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        shape=boston.shp,
        contour.plot = T,
        points.plot = F,legend = F)
text(x=boston.c$LON, y=boston.c$LAT,
     labels = boston.c$TOWN,
     cex=0.3)
# legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
# plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
# text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(boston.c$CMEDV),max(boston.c$CMEDV),l=6),2)),cex=1.5)
# rasterImage(legend_image, 0, 0, 1,1)
# dev.off()

coords <- cbind(boston.c$LON,
                boston.c$LAT)*100 # scaling the coordinates
#cnames <- c("(Intercept)","CRIM","ZN","INDUS","CHAS","NOX.2","RM.2","AGE","lDIS","lRAD","TAX","PTRATIO","B","lLSTAT")
N <- nrow(coords)
y <- boston.c$CMEDV
X <- matrix(1,nr=N)
# X <- as.matrix(cbind(1,
#                      boston.c$CRIM,
#                      boston.c$ZN,
#                      boston.c$INDUS,
#                      boston.c$CHAS,
#                      boston.c$NOX^2,
#                      boston.c$RM^2,
#                      boston.c$AGE,
#                      log(boston.c$DIS),
#                      log(boston.c$RAD),
#                      boston.c$TAX,
#                      boston.c$PTRATIO,
#                      boston.c$B,
#                      log(boston.c$LSTAT)))
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
cov.type <- "matern2"

#n.samples <- 1e4
#starting <- list("phi"=3/max(dist(coords))+0.01, "sigma.sq"=50, "tau.sq"=1,"nu"=5/2)
#tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1,"nu"=0)
#priors.1 <- list("beta.Norm"=list(rep(0,1), diag(1000,1)),
#                 "phi.Unif"=c(3/max(dist(coords)), 3/0.01), "sigma.sq.IG"=c(2, 2),
#                 "tau.sq.IG"=c(2, 0.1),"nu.Unif"=c(2,3))
#cov.model <- "matern"
#n.report <- 500
#m.1 <- spLM(y~1, coords=coords, starting=starting,
#            tuning=tuning, priors=priors.1, cov.model=cov.model,
#            n.samples=n.samples, verbose=T, n.report=n.report)

#burn.in <- 0.5*n.samples
#m.1 <- spRecover(m.1, start=burn.in, verbose=FALSE)
#round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
#round(summary(m.1$p.beta.recover.samples)$quantiles[c(3,1,5)],2)
#m.1.w.summary <- summary(mcmc(t(m.1$p.w.recover.samples)))$quantiles[,c(3,1,5)]

# random starting points
phis_init = .07 # runif(1,lower_phis,upper_phis)
sigma2_init = 50 # 1/rgamma(1,shape=shape_sigma,rate=1/scale_sigma)
tau2_init = 1 # 1/rgamma(1,shape=shape_tau,rate=1/scale_tau)
beta_init = rep(0,ncol(X))
nu_init=1.1
mc_sp <- hlmBayes_sp(coords=coords,
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
model_summary <- hlm_summary(chain = mc_sp, nburn = nburn, thin = 1)
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

gp0 <- lm(CMEDV ~ 1, data = boston.c)
summary(gp0)
coef;coefficients(gp0)

mat <- matrix(c(1,2,3), nr=1,nc=3, byrow=T)
layout(mat,
       widths = c(5,5,2),
       heights = c(3,3))
#sp_plot(11,"Spectral",cbind(coords,m.1.w.summary[,1]))
#points(coords,
#       col=sapply(m.1.w.summary$significance, function(x){
#         if(x==0) return("grey")
#         if(x==-1) return("blue")
#         else return("red")
#       }), pch=16,cex=0.5)
sp_plot(11,"Spectral",cbind(coords/100,y.hat),contour.plot = T,shape=boston.shp, points.plot = T,sig = z$sig, legend = F) 
sp_plot(11,"Spectral",cbind(coords/100,y), contour.plot = T,shape=boston.shp, points.plot = T)

# RMSE
par(mfcol=c(1,1))
sqrt(mean((y-y.hat)^2))
plot(y,y.hat)
plot(y-y.hat)

# Spatial Gradient and curvature assessment
grid.points <- as.matrix(expand.grid(seq(min(boston.c$LON),max(boston.c$LON),by=0.01),
                                     seq(min(boston.c$LAT),max(boston.c$LAT),by=0.01)),nc=2)
tmp <- over(SpatialPoints(grid.points,proj4string = CRS(proj4string(boston.shp))),boston.shp)
grid.points <- grid.points[!is.na(tmp),]
# sp_plot(11,"Spectral",cbind(coords/100,y),shape=boston.shp,contour.plot=T)
# points(grid.points, cex=0.2)
grid.points <- grid.points*100
samples <- (nburn+1):niter
gradient_est <- spatial_gradient(coords=coords,
                                 grid.points = grid.points,
                                 samples = samples,
                                 post_mcmc = mc_sp,
                                 cov.type = cov.type,
                                 ncores=5,
                                 nbatch = 100,return.mcmc = T)
# matern 2: 10
# matern 1: 100
grad.s1.hpd <- data.frame(gradient_est$grad1.est)
grad.s1.hpd$signif <- apply(grad.s1.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

grad.s2.hpd <- data.frame(gradient_est$grad2.est)
grad.s2.hpd$signif <- apply(grad.s2.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

grad.s11.hpd <- data.frame(gradient_est$grad11.est)
grad.s11.hpd$signif <- apply(grad.s11.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  
grad.s12.hpd <- data.frame(gradient_est$grad12.est)
grad.s12.hpd$signif <- apply(grad.s12.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  
grad.s22.hpd <- data.frame(gradient_est$grad22.est)
grad.s22.hpd$signif <- apply(grad.s22.hpd,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})  

# Posterior Surface Analysis
det_est <- eigen_est_1 <- eigen_est_2 <- laplace_est <- div_est <- matrix(NA,nr=5000,nc=nrow(grid.points))
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

cex.pt <- 0.7
cex.txt <- 1.5
# pdf("/Users/aritrahader/Desktop/spatiotemporal_gradient/space/boston-spatial.pdf",width=24,height=10)
mat <- matrix(c(1:6), nr=1,nc=6, byrow=T)
layout(mat,
       widths = c(rep(c(5,3),3)),
       heights = c(rep(3,6)))
rster.obj <- sp_plot(11,"Spectral",cbind(coords/100,z[,1]),
                     shape=boston.shp,
                     contour.plot = T,
                     xlab=latex2exp::TeX("Longitude$\\degree$"),
                     ylab = latex2exp::TeX("Latitude$\\degree$"),
                     points.plot=T,
                     raster.surf = T,
                     sig = z$sig)

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
# sp_plot(11,"Spectral",cbind(grid.points/100,det_est_val[,1]),
#         contour.plot = T, shape=boston.shp,
#         xlab=latex2exp::TeX("(c) $\\widehat{|\\nabla^2Y|}$"),
#         points.plot=T,sig=det_est_val$sig)

sp_plot(11,"Spectral",cbind(grid.points/100,div_est_val[,1]),
        contour.plot = T, shape=boston.shp,
        xlab=latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        points.plot=T,sig=div_est_val$sig,
        grid=F)

sp_plot(11,"Spectral",cbind(grid.points/100,laplace_est_val[,1]),
        contour.plot = T, shape=boston.shp,
        xlab=latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        points.plot=T,sig=laplace_est_val$sig,
        grid=F)

# sp_plot(11,"Spectral",cbind(coords/100,y.hat),
#         contour.plot = T, shape=boston.shp,
#         xlab=latex2exp::TeX("(f) $\\widehat{y}$"),
#         points.plot=T,sig=z$sig)
# dev.off()

rster.obj <- sp_plot(11,"Spectral",cbind(coords/100,z[,1]),
                     shape=boston.shp,
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
level.choice <- c(-15,15) #c(-10,10)
#c(x.levels[2],x.levels[length(x.levels)-1])
subset.points.1 <- subset(x,level==level.choice[1])
subset.points.2 <- subset(x,level==level.choice[2])
#1,2, 5
lines(subset.points.1@lines[[1]]@Lines[[1]]@coords,lwd=3.5)
lines(subset.points.2@lines[[1]]@Lines[[1]]@coords,lwd=3.5)
lines(subset.points.2@lines[[1]]@Lines[[2]]@coords,lwd=3.5)
lines(subset.points.2@lines[[1]]@Lines[[5]]@coords,lwd=3.5)
subset.points <- list(subset.points.1@lines[[1]]@Lines[[1]]@coords,
                      subset.points.2@lines[[1]]@Lines[[1]]@coords,
                      subset.points.2@lines[[1]]@Lines[[2]]@coords,
                      subset.points.2@lines[[1]]@Lines[[5]]@coords)
