require(spWombling)
require(coda)
require(parallel)
###########################
# Synthetic Spatial Data  #
###########################
# set.seed(123456)
N <- 100
tau <- 1
# synthetic location
coords <- matrix(runif(2*N,0,10),nr=N,nc=2) 
# synthetic pattern 
sig2 = 5
phis = 3
Delta = as.matrix(dist(coords))

y = MASS::mvrnorm(1, rep(0,N), sig2*(1+phis*sqrt(5)*Delta+phis^2*5*Delta^2/3)*exp(-phis*sqrt(5)*Delta))
# y <- rnorm(N,10*(sin(3*pi*coords[,1])+cos(3*pi*coords[,2])),tau)
sp_plot(11,"Spectral",cbind(coords,y), contour.plot=T, points.plot = T)

set.seed(NULL)
y=y # response
X=matrix(1,nr=N) # intercept
XtX=crossprod(X,X) 
D = 1:N
# cannot be randomly started
z_init= rep(0,N)
lower_phis=3/max(dist(coords))
upper_phis=30
lower_nu=2
upper_nu=3
shape_sigma=2
scale_sigma=0.07
shape_tau=2
scale_tau=0.1
mean_beta=rep(0,ncol(X))
prec_beta=1e-6*diag(ncol(X))
Delta=as.matrix(dist(coords)) # instead of coordinates supply dist
steps_init=1
stepnu_init=0.001
niter=1e4
nburn = niter/2
report=1e2
DtD <- diag(N)
cov.type <- "matern2"

phis_init = runif(1,lower_phis,upper_phis)
sigma2_init = 1/rgamma(1,shape=shape_sigma,rate=1/scale_sigma)
tau2_init = 1/rgamma(1,shape=shape_tau,rate=1/scale_tau)
beta_init = MASS::mvrnorm(1,mean_beta,1/prec_beta)
nu_init=2.5
nu_est <- F
mc_sp <- hlmBayes_sp(coords = coords,
                     y=y,
                     X=X,
                     z_init=z_init,
                     D = D,DtD=DtD,
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
                     nu_est = nu_est,
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

# Model Selection
#beta.est=coef["beta","median"]
#sigma2.est=coef["sigma2","median"]
#phi.est=coef["phis","median"]
#tau2.est=coef["tau2","median"]
#z.est=z[,"median"]
#round(unlist(model_selection(y=y,
#                             xmat=X,
#                             beta.est=beta.est,
#                             sigma2.est=sigma2.est,
#                             phi.est=phi.est,
#                             tau2.est=tau2.est,
#                             z.est=z.est,
#                             Delta=Delta,
#                             cov.type = cov.type)),2)

#,"2.5%","97.5%")
# par(mfcol=c(1,2))

# MSE
y.hat <- z[,"median"]+coef["post_beta","median"]
round(sqrt(mean((y-y.hat)^2)),2)


# Gradient Analysis
grid.points <- as.matrix(expand.grid(seq(0,1,by=0.05)[-c(1,21)],
                                     seq(0,1,by=0.05)[-c(1,21)]),nc=2)
colnames(grid.points) <- NULL
# plot(grid.points, pch=16, xlab="x",ylab="y")
# abline(h=seq(0,1,by=0.05)[-c(1,21)], lty=2, col=2)
# abline(v=seq(0,1,by=0.05)[-c(1,21)], lty=2, col=2)
samples <- (nburn+1):niter#+(i-1)*nburn
gradient_est <- spatial_gradient(coords=coords,
                                 grid.points = grid.points,
                                 samples = samples,
                                 post_mcmc = mc_sp,
                                 cov.type = cov.type,
                                 ncores=6,
                                 nbatch = 200,
                                 return.mcmc = T)
pat1.res <- list(coords=coords,y=y,mcmc=mc_sp,gradient=gradient_est)
save(pat1.res, file="~/Desktop/spatiotemporal_gradient/space/pat1.RData")
load("~/Desktop/spatiotemporal_gradient/space/pat1.RData")
grad.tval <- cbind.data.frame(s1=30*pi*cos(3*pi*grid.points[,1]),
                              s2=-30*pi*sin(3*pi*grid.points[,2]),
                              s11=-90*pi^2*sin(3*pi*grid.points[,1]),
                              s12=0,
                              s22=-90*pi^2*cos(3*pi*grid.points[,2]))
# grad.tval <- cbind.data.frame(s1=30*pi*cos(3*pi*grid.points[,1])*cos(3*pi*grid.points[,2]),
#                               s2=-30*pi*sin(3*pi*grid.points[,1])*sin(3*pi*grid.points[,2]),
#                               s11=-90*pi^2*sin(3*pi*grid.points[,1])*cos(3*pi*grid.points[,2]),
#                               s12=-90*pi^2*cos(3*pi*grid.points[,1])*sin(3*pi*grid.points[,2]),
#                               s22=-90*pi^2*sin(3*pi*grid.points[,1])*cos(3*pi*grid.points[,2]))

grad.s1.hpd <- data.frame(cbind(grad.tval[,1],gradient_est$grad1.est))
grad.s1.hpd$signif <- apply(grad.s1.hpd,1,function(x){
  if(x[3]>0 & x[4]>0) return (1)
  if(x[3]<0 & x[4]<0) return (-1)
  else return(0)
})

grad.s2.hpd <- data.frame(cbind(grad.tval[,2],gradient_est$grad2.est))
grad.s2.hpd$signif <- apply(grad.s2.hpd,1,function(x){
  if(x[3]>0 & x[4]>0) return (1)
  if(x[3]<0 & x[4]<0) return (-1)
  else return(0)
})

grad.s11.hpd <- data.frame(cbind(grad.tval[,3], gradient_est$grad11.est))
grad.s11.hpd$signif <- apply(grad.s11.hpd,1,function(x){
  if(x[3]>0 & x[4]>0) return (1)
  if(x[3]<0 & x[4]<0) return (-1)
  else return(0)
})

grad.s12.hpd <- data.frame(cbind(grad.tval[,4],gradient_est$grad12.est))
grad.s12.hpd$signif <- apply(grad.s12.hpd,1,function(x){
  if(x[3]>0 & x[4]>0) return (1)
  if(x[3]<0 & x[4]<0) return (-1)
  else return(0)
})

grad.s22.hpd <- data.frame(cbind(grad.tval[,5],gradient_est$grad22.est))
grad.s22.hpd$signif <- apply(grad.s22.hpd,1,function(x){
  if(x[3]>0 & x[4]>0) return (1)
  if(x[3]<0 & x[4]<0) return (-1)
  else return(0)
})
# Coverage Probability
cat(1-round(sum(apply(grad.s1.hpd,1, function(x){
  ifelse(x[1]>=x[3] & x[1]<=x[4],F,T) 
}))/361,3),
1-round(sum(apply(grad.s2.hpd,1, function(x){
  ifelse(x[1]>=x[3] & x[1]<=x[4],F,T) 
}))/361,3),
1-round(sum(apply(grad.s11.hpd,1, function(x){
  if(x[1]>=x[3] & x[1]<=x[4]) return(0)
  else return(1)
}))/361,3),
1-round(sum(apply(grad.s12.hpd,1, function(x){
  if(x[1]>=x[3] & x[1]<=x[4]) return(0)
  else return(1)
}))/361,3),
1-round(sum(apply(grad.s22.hpd,1, function(x){
  if(x[1]>=x[3] & x[1]<=x[4]) return(0)
  else return(1)
}))/361,3),"\n")

mat <- matrix(c(1:12), nr=2,nc=6, byrow=T)
layout(mat,
       widths = rep(c(5,3),3),
       heights = c(3,3))
# Process
sp_plot(11,"Spectral",cbind(coords,y), contour.plot=T, points.plot = T)
# True Gradient along s1: 
sp_plot(11,"Spectral",cbind(grid.points,grad.tval[,1]), contour.plot = T, points.plot = T)
# True Gradient along s2: 
sp_plot(11,"Spectral",cbind(grid.points,grad.tval[,2]), contour.plot = T, points.plot = T)
# Fitted Process
sp_plot(11,"Spectral",cbind(coords,z[,"median"]+coef["post_beta","median"]), contour.plot=T, points.plot = T, sig = z$sig)
# Estimated Gradient along s1
sp_plot(11,"Spectral",cbind(grid.points,grad.s1.hpd[,2]), contour.plot = T, points.plot = T,sig = grad.s1.hpd$signif)
# Estimated Gradient along s2
sp_plot(11,"Spectral",cbind(grid.points,grad.s2.hpd[,2]), contour.plot = T, points.plot = T, sig = grad.s2.hpd$signif)

mat <- matrix(c(1:12), nr=2,nc=6, byrow=T)
layout(mat,
       widths = rep(c(5,3),3),
       heights = rep(3,6))
# True Curvature
sp_plot(11,"Spectral",cbind(grid.points,grad.tval[,3]), contour.plot = T, points.plot = T)
sp_plot(11,"Spectral",cbind(grid.points,grad.tval[,4]), contour.plot = T, points.plot=T)
sp_plot(11,"Spectral",cbind(grid.points,grad.tval[,5]), contour.plot = T, points.plot = T)
# Estimated Curvature
sp_plot(11,"Spectral",cbind(grid.points,grad.s11.hpd[,2]), contour.plot = T, points.plot = T,sig = grad.s11.hpd$signif)
sp_plot(11,"Spectral",cbind(grid.points,grad.s12.hpd[,2]), contour.plot = T, points.plot = T, sig=grad.s12.hpd$signif)
sp_plot(11,"Spectral",cbind(grid.points,grad.s22.hpd[,2]), contour.plot = T, points.plot = T, sig = grad.s22.hpd$signif)

##############################
# Posterior Surface Analysis # integrate this into the spatial gradient function
##############################
det_est <- eigen_est_1 <- eigen_est_2 <- laplace_est <- div_est <- matrix(NA,nr=100,nc=nrow(grid.points))
for(i in 1:100){
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
det_est_val <- cbind.data.frame(est=apply(det_est,2,median),HPDinterval(as.mcmc(det_est)))
det_est_val$sig <- apply(det_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
eigen_est_val1 <- cbind.data.frame(est=apply(eigen_est_1,2,median),HPDinterval(as.mcmc(eigen_est_1)))
eigen_est_val1$sig <- apply(eigen_est_val1,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
eigen_est_val2 <- cbind.data.frame(est=apply(eigen_est_2,2,median),HPDinterval(as.mcmc(eigen_est_2)))
eigen_est_val2$sig <- apply(eigen_est_val2,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
laplace_est_val <- cbind.data.frame(est=apply(laplace_est,2,median),HPDinterval(as.mcmc(laplace_est)))
laplace_est_val$sig <- apply(laplace_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
div_est_val <- cbind.data.frame(est=apply(div_est,2,median),HPDinterval(as.mcmc(div_est)))
div_est_val$sig <- apply(div_est_val,1,function(x){
  if(x[2]>0 & x[3]>0) return (1)
  if(x[2]<0 & x[3]<0) return (-1)
  else return(0)
})
det_true <- grad.s11.hpd[,1]*grad.s22.hpd[,1]-(grad.s12.hpd[,1])^2
eigen_true <- t(apply(cbind.data.frame(grad.s11.hpd[,1],(grad.s12.hpd[,1]),grad.s22.hpd[,1]),1, function(x){
  eigen(matrix(c(x[1],x[2],x[2],x[3]), nc=2,nr=2, byrow=T))$values
}))
div_true <- grad.s1.hpd[,1]+grad.s2.hpd[,1]
laplace_true <- grad.s11.hpd[,1]+grad.s22.hpd[,1]

pdf("/Users/aritrahader/Dropbox/PhD Research/paper 2/plots/post-surf-2.pdf",
    width=9.5)
mat <- matrix(c(1:12), nr=2,nc=6, byrow=T)
layout(mat,
       widths = rep(c(5,3),3),
       heights = rep(3,6))
sp_plot(11,"Spectral",cbind(grid.points,eigen_true[,1]), contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(a) $\\lambda_1$"))
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(eigen_true[,1]),max(eigen_true[,1]),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)
sp_plot(11,"Spectral",cbind(grid.points,eigen_true[,2]),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(b) $\\lambda_2$"))
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(eigen_true[,2]),max(eigen_true[,2]),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)
sp_plot(11,"Spectral",cbind(grid.points,det_true),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(c) $|\\nabla^2Y|$"))
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(det_true)/1e4,max(det_true)/1e4,l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)
sp_plot(11,"Spectral",cbind(grid.points,div_true),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(d) $div Y=\\nabla_1 Y+\\nabla_2 Y$"))
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(div_true),max(div_true),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)
sp_plot(11,"Spectral",cbind(grid.points,laplace_true),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(e) $\\Delta Y=\\nabla^2_{11}Y+\\nabla^2_{22}Y$"))
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(laplace_true)/10,max(laplace_true)/10,l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)
sp_plot(11,"Spectral",cbind(coords,y),
        xlab=latex2exp::TeX("(f) $y$"),contour.plot = T,ylab="")
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(y),max(y),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)


sp_plot(11,"Spectral",cbind(grid.points,eigen_est_val1[,1]),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(a) $\\widehat{\\lambda_1}$"))
points(grid.points,
       col=sapply(eigen_est_val1$sig, function(x){
         if(x==0) return("white")
         if(x==-1) return("cyan")
         else return("green")
       }), pch=16,cex=0.7)
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(eigen_est_val1[,1]),max(eigen_est_val1[,1]),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)

sp_plot(11,"Spectral",cbind(grid.points,eigen_est_val2[,1]),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(b) $\\widehat{\\lambda_2}$"))
points(grid.points,
       col=sapply(eigen_est_val2$sig, function(x){
         if(x==0) return("white")
         if(x==-1) return("cyan")
         else return("green")
       }), pch=16,cex=0.7)
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(eigen_est_val2[,1]),max(eigen_est_val2[,1]),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)

sp_plot(11,"Spectral",cbind(grid.points,det_est_val[,1]),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(c) $\\widehat{|\\nabla^2Y|}$"))
points(grid.points,
       col=sapply(det_est_val$sig, function(x){
         if(x==0) return("white")
         if(x==-1) return("cyan")
         else return("green")
       }), pch=16,cex=0.7)
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(det_est_val[,1])/1e4,max(det_est_val[,1])/1e4,l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)

sp_plot(11,"Spectral",cbind(grid.points,div_est_val[,1]),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(d) $\\widehat{div Y}=\\widehat{\\nabla_1 Y}+\\widehat{\\nabla_2 Y}$"))
points(grid.points,
       col=sapply(div_est_val$sig, function(x){
         if(x==0) return("white")
         if(x==-1) return("cyan")
         else return("green")
       }), pch=16,cex=0.7)
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(div_est_val[,1]),max(div_est_val[,1]),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)

sp_plot(11,"Spectral",cbind(grid.points,laplace_est_val[,1]),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(e) $\\widehat{\\Delta Y}=\\widehat{\\nabla^2_{11}Y}+\\widehat{\\nabla^2_{22}Y}$"))
points(grid.points,
       col=sapply(laplace_est_val$sig, function(x){
         if(x==0) return("white")
         if(x==-1) return("cyan")
         else return("green")
       }), pch=16,cex=0.7)
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(laplace_est_val[,1])/10,max(laplace_est_val[,1])/10,l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)

sp_plot(11,"Spectral",cbind(coords,y.hat),contour.plot = T,ylab="",
        xlab=latex2exp::TeX("(f) $\\widehat{y}$"))
points(coords,
       col=sapply(z$sig, function(x){
         if(x==0) return("white")
         if(x==-1) return("cyan")
         else return("green")
       }), pch=16,cex=0.7)
legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(y.hat),max(y.hat),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

# Check
##s1
par(mfrow=c(2,3))
y.hat <- z[,"median"]+coef["post_beta","median"]
y.upper.hpd <- z[order(y.hat),"upper.hpd"]+coef["post_beta","upper.hpd"]
y.lower.hpd <- z[order(y.hat),"lower.hpd"]+coef["post_beta","lower.hpd"]

plot(y.hat,y, type="n",
     ylim=c(min(y.lower.hpd),
            max(y.upper.hpd)),
     xlab=latex2exp::TeX("$\\widehat{y}(s)$"),
     ylab=latex2exp::TeX("$y(s)$")) 
polygon(c(y.hat[order(y.hat)],rev(y.hat[order(y.hat)])),
        c(y.upper.hpd, rev(y.lower.hpd)),
        col = 'grey80', border = NA)
points(y.hat,y)
abline(c(0,1),col="blue")
lines(y.hat[order(y.hat)],y.lower.hpd, lty="dashed", col="red")
lines(y.hat[order(y.hat)],y.upper.hpd, lty="dashed", col="red")
grid()
#legend("bottomright", inset=0.05, c("hpd","0-1 line"), lty=c("dashed","solid"),col=c("red","blue"))

grad.s1.upper.hpd <- grad.s1.hpd[,4]
grad.s1.hat <- grad.s1.hpd[,2]
grad.s1.lower.hpd <- grad.s1.hpd[,3]

plot(grad.s1.hat,grad.s1.hpd[,1], type="n",
     ylim=c(min(grad.s1.lower.hpd),max(grad.s1.upper.hpd)),
     xlab=latex2exp::TeX("$\\widehat{\\nabla_1 y}(s)$"),
     ylab=latex2exp::TeX("$\\nabla_1 y(s)$")) 
polygon(c(grad.s1.hat[order(grad.s1.hat)],
          rev(grad.s1.hat[order(grad.s1.hat)])),
        c(grad.s1.upper.hpd[order(grad.s1.hat)],
          rev(grad.s1.lower.hpd[order(grad.s1.hat)])),
        col = 'grey80', border = NA)
points(grad.s1.hat,grad.s1.hpd[,1])
abline(c(0,1),col="blue")
lines(grad.s1.hat[order(grad.s1.hat)],
      grad.s1.lower.hpd[order(grad.s1.hat)],
      lty="dashed", col="red")
lines(grad.s1.hat[order(grad.s1.hat)],
      grad.s1.upper.hpd[order(grad.s1.hat)],
      lty="dashed", col="red")
grid()

##s2
grad.s2.upper.hpd <- grad.s2.hpd[,4]
grad.s2.hat <- grad.s2.hpd[,2]
grad.s2.lower.hpd <- grad.s2.hpd[,3]

plot(grad.s2.hat,grad.s2.hpd[,1], type="n",
     ylim=c(min(grad.s2.lower.hpd),max(grad.s2.upper.hpd)),
     xlab=latex2exp::TeX("$\\widehat{\\nabla_2 y}(s)$"),
     ylab=latex2exp::TeX("$\\nabla_2 y(s)$")) 
polygon(c(grad.s2.hat[order(grad.s2.hat)],
          rev(grad.s2.hat[order(grad.s2.hat)])),
        c(grad.s2.upper.hpd[order(grad.s2.hat)],
          rev(grad.s2.lower.hpd[order(grad.s2.hat)])),
        col = 'grey80', border = NA)
points(grad.s2.hat,grad.s2.hpd[,1])
abline(c(0,1),col="blue")
lines(grad.s2.hat[order(grad.s2.hat)],
      grad.s2.lower.hpd[order(grad.s2.hat)],
      lty="dashed", col="red")
lines(grad.s2.hat[order(grad.s2.hat)],
      grad.s2.upper.hpd[order(grad.s2.hat)],
      lty="dashed", col="red")
grid()

# plot(c(-3,3),c(-3,3),type = 'n', axes = F,xlab = '', ylab = '', main = '')
# text(x = 0,y=2,paste("RMSE(s1)",round(sqrt(mean((grad.s1.hpd[,1]-grad.s1.hpd[,2])^2)),2),sep=": "))
# text(x = 0,y=1,paste("RMSE(s2)",round(sqrt(mean((grad.s2.hpd[,1]-grad.s2.hpd[,2])^2)),2),sep=": "))
# text(x = 0,y=0,paste("RMSE(s11)",round(sqrt(mean((grad.s11.hpd[,1]-grad.s11.hpd[,2])^2)),2),sep=": "))
# text(x = 0,y=-1,paste("RMSE(s22)",round(sqrt(mean((grad.s12.hpd[,1]-grad.s12.hpd[,2])^2)),2),sep=": "))
# text(x = 0,y=-2,paste("RMSE(s12)",round(sqrt(mean((grad.s22.hpd[,1]-grad.s22.hpd[,2])^2)),2),sep=": "))

grad.s11.upper.hpd <- grad.s11.hpd[,4]
grad.s11.hat <- grad.s11.hpd[,2]
grad.s11.lower.hpd <- grad.s11.hpd[,3]

plot(grad.s11.hat,grad.s11.hpd[,1],type="n",
     ylim=c(min(grad.s11.lower.hpd),max(grad.s11.upper.hpd)),
     xlab=latex2exp::TeX("$\\widehat{\\nabla_{11} y}(s)$"),
     ylab=latex2exp::TeX("$\\nabla_{11} y(s)$"))
polygon(c(grad.s11.hat[order(grad.s11.hat)],
          rev(grad.s11.hat[order(grad.s11.hat)])),
        c(grad.s11.upper.hpd[order(grad.s11.hat)],
          rev(grad.s11.lower.hpd[order(grad.s11.hat)])),
        col = 'grey80', border = NA)
points(grad.s11.hat,grad.s11.hpd[,1])
abline(c(0,1),col="blue")
lines(grad.s11.hat[order(grad.s11.hat)],
      grad.s11.lower.hpd[order(grad.s11.hat)],
      lty="dashed", col="red")
lines(grad.s11.hat[order(grad.s11.hat)],
      grad.s11.upper.hpd[order(grad.s11.hat)],
      lty="dashed", col="red")
grid()
grad.s12.upper.hpd <- grad.s12.hpd[,4]
grad.s12.hat <- grad.s12.hpd[,2]
grad.s12.lower.hpd <- grad.s12.hpd[,3]

plot(grad.s12.hat,grad.s12.hpd[,1],type="n",
     ylim=c(min(grad.s12.lower.hpd),max(grad.s12.upper.hpd)),
     xlab=latex2exp::TeX("$\\widehat{\\nabla_{12} y}(s)$"),
     ylab=latex2exp::TeX("$\\nabla_{12} y(s)$"))
polygon(c(grad.s12.hat[order(grad.s12.hat)],
          rev(grad.s12.hat[order(grad.s12.hat)])),
        c(grad.s12.upper.hpd[order(grad.s12.hat)],
          rev(grad.s12.lower.hpd[order(grad.s12.hat)])),
        col = 'grey80', border = NA)
points(grad.s12.hat,grad.s12.hpd[,1])
abline(c(0,1),col="blue")
lines(grad.s12.hat[order(grad.s12.hat)],
      grad.s12.lower.hpd[order(grad.s12.hat)],
      lty="dashed", col="red")
lines(grad.s12.hat[order(grad.s12.hat)],
      grad.s12.upper.hpd[order(grad.s12.hat)],
      lty="dashed", col="red")

grid()

grad.s22.upper.hpd <- grad.s22.hpd[,4]
grad.s22.hat <- grad.s22.hpd[,2]
grad.s22.lower.hpd <- grad.s22.hpd[,3]

plot(grad.s22.hat,grad.s22.hpd[,1],type="n",
     ylim=c(min(grad.s22.lower.hpd),max(grad.s22.upper.hpd)),
     xlab=latex2exp::TeX("$\\widehat{\\nabla_{22} y}(s)$"),
     ylab=latex2exp::TeX("$\\nabla_{22} y(s)$"))
polygon(c(grad.s22.hat[order(grad.s22.hat)],
          rev(grad.s22.hat[order(grad.s22.hat)])),
        c(grad.s22.upper.hpd[order(grad.s22.hat)],
          rev(grad.s22.lower.hpd[order(grad.s22.hat)])),
        col = 'grey80', border = NA)
points(grad.s22.hat,grad.s22.hpd[,1])
abline(c(0,1),col="blue")
lines(grad.s22.hat[order(grad.s22.hat)],
      grad.s22.lower.hpd[order(grad.s22.hat)],
      lty="dashed", col="red")
lines(grad.s22.hat[order(grad.s22.hat)],
      grad.s22.upper.hpd[order(grad.s22.hat)],
      lty="dashed", col="red")
grid()


#################################
# Wombling: Gradient Assessment #
#################################
# pdf("/Users/aritrah/OneDrive - University of Connecticut/PhD Research/paper 2/plots/womble.pdf",width=5,height = 5.5)
# sp_plot(11,"Spectral",cbind(coords,y))
####################### 
# creating boundaries #   
#######################
# select a trough boundary:: closed
# boundary_a <- locator(10)
# boundary_a <- t(do.call(rbind,boundary_a))
# points(boundary_a, pch=16, cex=0.5)
# lines(rbind(boundary_a,boundary_a[1,]))

# select a peak boundary:: closed
# boundary_b <- locator(10)
# boundary_b <- t(do.call(rbind,boundary_b))
# points(boundary_b, pch=16, cex=0.5)
# lines(rbind(boundary_b,boundary_b[1,]))

# select a contour:: closed
# boundary_c <- locator(10)
# boundary_c <- t(do.call(rbind,boundary_c))
# points(boundary_c, pch=16, cex=0.5)
# lines(rbind(boundary_c,boundary_c[1,]))

# select a contour:: open
# boundary_d <- locator(10)
# boundary_d <- t(do.call(rbind,boundary_d))
# points(boundary_d, pch=16, cex=0.5)
# lines(rbind(boundary_d))
# dev.off()
#pattern.1 <- list(coords=coords,
#                  y=y,
#                  ba.c=boundary_a,
#                  bb.c=boundary_b,
#                  bc.c=boundary_c,
#                  bd.o=boundary_d)
# save(pattern.1, file="/Users/aritrah/OneDrive - University of Connecticut/PhD Research/paper 2/plots/pattern-1.RData")
#load("~/Dropbox/PhD Research/paper 2/plots/pattern-1.RData")
#post_mcmc <- recover_mcmc(mcmc.obj = results,
#                          m.ind=c(0,0,0,0,1),
#                          nburn=1,
#                          combine.chains = T)

mat <- matrix(c(1:2), nr=1,nc=2, byrow=T)
layout(mat,
       widths = c(5,2),
       heights = c(3,3))

rastr.obj <- sp_plot(11,"Spectral",cbind(coords,y), contour.plot = T,raster.surf = T, legend = F)
x <- raster::rasterToContour(rastr.obj,nlevels=20)
# Test for wombling
x.levels <- as.numeric(as.character(x$level))
level.choice <- c(x.levels[2],x.levels[length(x.levels)-1])
subset.points.1 <- subset(x,level==level.choice[1])
subset.points.2 <- subset(x,level==level.choice[2])
x <- raster::rasterToContour(rastr.obj,nlevels=10)
subset.points.3 <- subset(x,level==15)
# curve-1
lines(subset.points.1@lines[[1]]@Lines[[1]]@coords,lwd=2.5)
# curve-2
lines(subset.points.2@lines[[1]]@Lines[[3]]@coords,lwd=2.5)
# curve 3
lines(subset.points.3@lines[[1]]@Lines[[3]]@coords,lwd=2.5)
pts.hull <- locator(20)
subset.points.4 <- bezierCurve(pts.hull$x,pts.hull$y,100)
# curve 4
lines(subset.points.4, type="l",lwd=2.5)
points(pts.hull, cex=0.5, col="white")
points(pts.hull, pch="+", cex=0.5)

legend_image <- as.raster(matrix(colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100), ncol=1))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=2, y = seq(0.01,0.99,l=6), labels = sprintf("%.2f",round(seq(min(y),max(y),l=6),2)))
rasterImage(legend_image, 0, 0, 1,1)

subset.points <- list(subset.points.1@lines[[1]]@Lines[[1]]@coords,
                   subset.points.2@lines[[1]]@Lines[[3]]@coords,
                   subset.points.3@lines[[1]]@Lines[[3]]@coords,
                   t(do.call(rbind,subset.points.4)))





gradient_est <- list()
for(i in 1:length(subset.points)){
  grid.points <- subset.points[[i]]
  gradient_est[[i]] <- spatial_gradient(coords=coords,
                                   grid.points = grid.points,
                                   post_mcmc = mc_sp,
                                   cov.type = cov.type,
                                   ncores=5,
                                   nbatch = 100,
                                   return.mcmc = T)
  cat("Boundary:",i,"\n")
}


cwomble_results <- c() 
for(i.1 in 1:length(subset.points)){
  x <- gradient_est[[i.1]]
  grid.points <- subset.points[[i.1]]
  norm_dir <- matrix(NA, nc=3,nr=nrow(grid.points)-1)
  for(i in 1:(nrow(grid.points)-1)){
    # compute normal direction
    tmp.n <- grid.points[(i+1),]-grid.points[i,]
    tmp.n <- rev(tmp.n/sqrt(sum(tmp.n^2)))
    tmp.n[1] <- -tmp.n[1]
    # compute line-segement length
    norm_dir[i,] <- c(tmp.n,dist(grid.points[c(i:(i+1)),]))
  }
  grad_true <- cbind.data.frame(grad.x.true=30*pi*cos(3*pi*grid.points[,1]),
                                grad.y.true=-30*pi*sin(3*pi*grid.points[,2]))
  grad_bndry <- matrix(NA,nc=nrow(grid.points)-1,nr=nrow(x$grad.s1.mcmc))
  for(i in 1:nrow(x$grad.s1.mcmc)){
    grad.bndry <- cbind(x$grad.s1.mcmc[i,-nrow(grid.points)],x$grad.s2.mcmc[i,-nrow(grid.points)])
    grad_bndry[i,] <- apply(cbind.data.frame(norm_dir[,c(1,2)],grad.bndry), 1,function(x){
      c(x[1],x[2])%*%c(x[3],x[4])
    })
  }
  grad_bndry <- as.mcmc(grad_bndry)
  grad_bndry_est <- cbind.data.frame(est=apply(grad_bndry,2,median),HPDinterval(grad_bndry))
  
  grad_bndry_est$sig <- apply(grad_bndry_est,1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })
  id.sig <- grad_bndry_est$sig%in%c(1,-1)
  col.segment <- sapply(grad_bndry_est$sig, function(x){
    if(x==0) "white"
    else if(x==1) "green"
    else if(x==-1) "cyan"
  })
  for(i in 1:(nrow(grid.points)-1)) lines(rbind(grid.points[i,],grid.points[(i+1),]), col=col.segment[i],lwd=3)
  
  # Average Gradient: computed only on the significant gradients
  true_val <- sum(apply(cbind.data.frame(norm_dir[,c(1,2)],grad_true[-nrow(grid.points),]), 1,function(x){
    c(x[1],x[2])%*%c(x[3],x[4])
  }))/sum(norm_dir[,3])
  grad.eval <- grad_true[-nrow(grid.points),]
  true_val_sig <- sum(apply(cbind.data.frame(norm_dir[id.sig,c(1,2)],grad.eval[id.sig,]), 1,function(x){
    c(x[1],x[2])%*%c(x[3],x[4])
  }))/sum(norm_dir[id.sig,3])
  avg_grad_bdry <- as.mcmc(apply(grad_bndry,1,sum)/sum(norm_dir[,3]))
  avg_grad_bdry_est <- c(median(avg_grad_bdry),apply(HPDinterval(grad_bndry),2,sum)/sum(norm_dir[,3]))
  avg_grad_bdry_sig <- as.mcmc(apply(grad_bndry[,id.sig],1,sum)/sum(norm_dir[,3]))
  avg_grad_bdry_est_sig <- c(median(avg_grad_bdry_sig),apply(HPDinterval(grad_bndry[,id.sig]),2,sum)/sum(norm_dir[,3]))
  cwomble_results <- rbind(cwomble_results,
                           c(true_val,avg_grad_bdry_est,avg_grad_bdry_est_sig,true_val_sig))
  
}

round(cwomble_results/1e3,2)

cwomble1_results <- c()
for(i.1 in 1:length(subset.points)){
    x <- gradient_est[[i.1]]
    # computing the unit Normal direction 
    # and the line segment lengths
    grid.points <- subset.points[[i.1]]
    norm_dir <- matrix(NA, nc=3,nr=nrow(grid.points)-1)
    for(i in 1:(nrow(grid.points)-1)){
      # compute normal direction
      tmp.n <- grid.points[(i+1),]-grid.points[i,]
      tmp.n <- rev(tmp.n/sqrt(sum(tmp.n^2)))
      tmp.n[2] <- -tmp.n[2]
      # compute line-segement length
      norm_dir[i,] <- c(tmp.n,dist(grid.points[c(i:(i+1)),]))
    }
    
    curve_bndry_true <- matrix(NA,nc=nrow(grid.points)-1,nr=1)
    for(i in 1:nrow(norm_dir)){
      hess.true <- matrix(c(-90*pi^2*sin(3*pi*grid.points[i,1]),0,0,-90*pi^2*cos(3*pi*grid.points[i,2])),nr=2,nc=2,byrow=T)
      curve_bndry_true[1,i] <- t(norm_dir[i,1:2])%*%hess.true%*%norm_dir[i,1:2]
    }
    
    curv_bndry <- matrix(NA,nc=nrow(grid.points)-1,nr=nrow(x$grad.s1.mcmc))
    for(i in 1:nrow(x$grad.s11.mcmc)){
      for(j in 1:nrow(norm_dir)){
        hess.bndry <- matrix(c(x$grad.s11.mcmc[i,j],
                               x$grad.s12.mcmc[i,j],
                               x$grad.s12.mcmc[i,j],
                               x$grad.s22.mcmc[i,j]),nc=2,nr=2,byrow=T)
        curv_bndry[i,j] <- t(norm_dir[j,1:2])%*%hess.bndry%*%norm_dir[j,1:2]
      }
    }
    curv_bndry <- as.mcmc(curv_bndry)
    curv_bndry_est <- cbind.data.frame(est=apply(curv_bndry,2,median),HPDinterval(curv_bndry))
    curv_bndry_est$sig <- apply(curv_bndry_est,1,function(x){
      if(x[2]>0 & x[3]>0) return (1)
      if(x[2]<0 & x[3]<0) return (-1)
      else return(0)
    })
    id.sig <- curv_bndry_est$sig%in%c(-1,1)
    col.segment <- sapply(curv_bndry_est$sig, function(x){
      if(x==0) "white"
      else if(x==1) "green"
      else if(x==-1) "cyan"
    })
    for(i in 1:(nrow(grid.points)-1)) lines(rbind(grid.points[i,],grid.points[(i+1),]), col=col.segment[i],lwd=3)
    # Average Curvature: computed only on the significant curvature
    true_val <- apply(curve_bndry_true,1,sum)/sum(norm_dir[,3])
    true_val_sig <- sum(curve_bndry_true[,id.sig])/sum(norm_dir[id.sig,3])
    avg_curv_bdry <- as.mcmc(apply(curv_bndry,1,sum)/sum(norm_dir[,3]))
    avg_curv_bdry_est <- c(median(avg_curv_bdry),apply(HPDinterval(curv_bndry),2,sum)/sum(norm_dir[,3]))
    avg_curv_bdry_sig <- as.mcmc(apply(curv_bndry[,id.sig],1,sum)/sum(norm_dir[id.sig,3]))
    if(sum(id.sig)!=0) avg_curv_bdry_est_sig <- c(median(avg_curv_bdry_sig),apply(HPDinterval(curv_bndry[,id.sig]),2,sum)/sum(norm_dir[id.sig,3]))
    else avg_curv_bdry_est_sig <- c(0,0,0)
    cwomble1_results <- rbind(cwomble1_results,
                              c(true_val,avg_curv_bdry_est,avg_curv_bdry_est_sig,true_val_sig))
}
round(cwomble1_results/1e4,2)

round(apply(do.call(rbind,lapply(gradient_est, function(x, grid.points){
  grid.points <- rbind(grid.points,grid.points[1,])
  norm_dir <- matrix(NA, nc=3,nr=nrow(grid.points)-1)
  for(i in 1:(nrow(grid.points)-1)){
    tmp.n <- grid.points[(i+1),]-grid.points[i,]
    tmp.n <- rev(tmp.n/sqrt(sum(tmp.n^2)))
    tmp.n[1] <- -tmp.n[1]
    # compute line-segement length
    norm_dir[i,] <- c(tmp.n,dist(grid.points[c(i:(i+1)),]))
  }
  
  grad.true <- cbind.data.frame(30*pi*cos(3*pi*grid.points[-11,1]),
                                -30*pi*sin(3*pi*grid.points[-11,2]))
  grad.ba <- cbind(x$grad1.est[,1], x$grad2.est[,1])
  grad.ba.lhpd <- cbind(x$grad1.est[,2], x$grad2.est[,2])
  grad.ba.uhpd <- cbind(x$grad1.est[,3], x$grad2.est[,3])
  true_val <- sum(apply(cbind.data.frame(norm_dir[,c(1,2)],grad.true), 1,function(x){
    c(x[1],x[2])%*%c(x[3],x[4])
  }))/sum(norm_dir[,3])
  est_val <- sum(apply(cbind.data.frame(norm_dir[,c(1,2)],grad.ba), 1,function(x){
    c(x[1],x[2])%*%c(x[3],x[4])
  }))/sum(norm_dir[,3])
  hpd <- c(sum(apply(cbind.data.frame(norm_dir[,c(1,2)],grad.ba.uhpd), 1,function(x){
    c(x[1],x[2])%*%c(x[3],x[4])
  }))/sum(norm_dir[,3]),sum(apply(cbind.data.frame(norm_dir[,c(1,2)],grad.ba.lhpd), 1,function(x){
    c(x[1],x[2])%*%c(x[3],x[4])
  }))/sum(norm_dir[,3]))
  
  post_surf <- matrix(NA,nc=4,nr=nrow(norm_dir))
  for(i in 1:nrow(norm_dir)){
    hess.true <- matrix(c(-90*pi^2*sin(3*pi*grid.points[i,1]),0,0,-90*pi^2*cos(3*pi*grid.points[i,2])),nc=2,nr=2,byrow=T)
    hess.ba <- matrix(c(x$grad11.est[i,1],
                        x$grad12.est[i,1],
                        x$grad12.est[i,1],
                        x$grad22.est[i,1]),nc=2,nr=2,byrow=T)
    hess.lhpd <- matrix(c(x$grad11.est[i,2],
                          x$grad12.est[i,2],
                          x$grad12.est[i,2],
                          x$grad22.est[i,2]),nc=2,nr=2,byrow=T)
    hess.uhpd <- matrix(c(x$grad11.est[i,3],
                          x$grad12.est[i,3],
                          x$grad12.est[i,3],
                          x$grad22.est[i,3]),nc=2,nr=2,byrow=T)
    true_val <- t(norm_dir[i,1:2])%*%hess.true%*%norm_dir[i,1:2]
    est_val <- t(norm_dir[i,1:2])%*%hess.ba%*%norm_dir[i,1:2]
    hpd <- c(t(norm_dir[i,1:2])%*%hess.lhpd%*%norm_dir[i,1:2],t(norm_dir[i,1:2])%*%hess.uhpd%*%norm_dir[i,1:2])
    post_surf[i,] <- c(true_val,est_val,hpd)
  }
  c(c(true_val,est_val,hpd[order(hpd)]),
        apply(post_surf,2,median)/sum(norm_dir[,3]))
}, grid.points=grid.points)),2,median),2)
