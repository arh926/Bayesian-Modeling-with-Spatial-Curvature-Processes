require(parallel)
require(sp)
require(spData)
require(spWombling)
data(meuse)
data(meuse.riv)
data(meuse.area)
# set coordinate system
coordinates(meuse) <- c("x", "y")
summary(meuse)
meuse <- spTransform(meuse, CRS("+proj=utm +zone=32 +datum=WGS84"))

meuse.shp <- raster::spPolygons(meuse.area,crs = CRS("+proj=utm +zone=32 +datum=WGS84"))

# pdf("~/Dropbox/PhD Research/paper 2/plots/meuse-spatial.pdf",width = 15, height=13)
mat <- matrix(c(1:8), nr=1,nc=8, byrow=T)
layout(mat,
       widths = rep(c(5,2),4),
       heights = rep(3,4))
for(metal in c("cadmium","copper","lead","zinc")){
  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"),)
  polygon(c(meuse.riv[1:88,1],meuse.riv[89:176,1]),
          c(meuse.riv[1:88,2],meuse.riv[89:176,2]),
          col = 'lightblue', border = "red")
  sp_plot(col.seq.length = 11,
          col.text = "Spectral",
          data_frame =  cbind(meuse@coords,meuse@data[,metal]),
          points.plot=T, legend = T,contour.plot = T,shape = meuse.shp,add = T)
}
# dev.off()

meuse <- meuse[-which(apply(meuse@data,1,function(x)sum(is.na(x)))==1),]
coords <- meuse@coords/10
Delta <- as.matrix(dist(coords))
y <- cbind(meuse@data[,"zinc"]/100,
           meuse@data[,"cadmium"],
           meuse@data[,"copper"]/10,
           meuse@data[,"lead"]/10)
X <- cbind(`(Intercept)`=1,
           elev = meuse@data$elev,
           om = meuse@data$om,
           dist.m = meuse@data$dist.m,
           flood.2 = as.numeric(meuse@data$ffreq==2),
           flood.3 = as.numeric(meuse@data$ffreq==3),
           soil.2  = as.numeric(meuse@data$soil==2),
           soil.3  = as.numeric(meuse@data$soil==3),
           lime.1 = as.numeric(meuse@data$lime==1))
N <- nrow(meuse@data)
D = 1:N
# cannot be randomly started
z_init= rep(0,N)
lower_phis=3/max(dist(coords))
upper_phis=30
lower_nu=1
upper_nu=3
shape_sigma=2
scale_sigma=1
shape_tau=2
scale_tau=100
mean_beta=rep(0,ncol(X))
prec_beta=1e-6*diag(ncol(X))
Delta=as.matrix(dist(coords)) # instead of coordinates supply dist
steps_init=1
stepnu_init=1
niter=1e4
nburn=niter/2
report=1e2
DtD <- diag(N)
cov.type <- "matern2"

phis_init = .007#runif(1,lower_phis,upper_phis)
sigma2_init = 50#1/rgamma(1,shape=shape_sigma,rate=1/scale_sigma)
tau2_init = 1#1/rgamma(1,shape=shape_tau,rate=1/scale_tau)
beta_init = rep(0,ncol(X))
nu_init=1.1
res_mc_sp <- mclapply(1:ncol(y), function(x){
  metal.conc <- y[,x]
  mc_sp <- hlmBayes_sp(coords = coords,
                       y=metal.conc,
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
},mc.cores = 4)
model_summary <- lapply(res_mc_sp, function(x) hlm_summary(chain = x,nburn = nburn,thin = 1))
coef <- lapply(model_summary, function(x) x$summary.pars)
lapply(coef,function(x) print(round(x,4)))
round(do.call(cbind,coef),4)

z <- lapply(model_summary, function(x) x$summary.latent)
for(i in 1:length(z)){
  z[[i]]$sig <- apply(z[[i]],1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })
}

for(i in 1:length(z)){
  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"),)
  polygon(c(meuse.riv[1:88,1],meuse.riv[89:176,1]),
          c(meuse.riv[1:88,2],meuse.riv[89:176,2]),
          col = 'lightblue', border = "red")
  sp_plot(col.seq.length = 11,
          col.text = "Spectral",
          data_frame =  cbind(meuse@coords,z[[i]][,"median"]),
          points.plot=T, legend = T,contour.plot = T,shape = meuse.shp,add = T,sig = z[[i]]$sig)
}
meuse.buffer <- rgeos::gBuffer(meuse.shp,width = 200,byid = T)
grid.meuse <- expand.grid(seq(178000,182000, by=200), seq(329000,334000, by=200))
tmp <- over(SpatialPoints(grid.meuse,proj4string = CRS(proj4string(meuse.buffer))),meuse.buffer)
grid.meuse <- matrix(unlist(grid.meuse[!is.na(tmp),]), nc=2, byrow = F)

# points(grid.meuse)
# plot.new()

samples <- (nburn+1):niter
gradient_est <- lapply(res_mc_sp, function(x){
  spatial_gradient(coords=coords,
                   grid.points = grid.meuse/10,
                   samples = samples,
                   post_mcmc = x,
                   cov.type = cov.type,
                   ncores=5,
                   nbatch = 100,
                   return.mcmc = T)
})
mat <- matrix(c(1:12), nr=2,nc=6, byrow=T)
layout(mat,
       widths = rep(c(5,3),3),
       heights = rep(3,6))
hvy_metals <- c("Cadmium","Copper","Lead","Zinc")
for(i in 1:length(gradient_est)){
  grad.s1.hpd <- data.frame(gradient_est[[i]]$grad1.est)
  grad.s1.hpd$signif <- apply(grad.s1.hpd,1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })  
  
  grad.s2.hpd <- data.frame(gradient_est[[i]]$grad2.est)
  grad.s2.hpd$signif <- apply(grad.s2.hpd,1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })  
  
  grad.s11.hpd <- data.frame(gradient_est[[i]]$grad11.est)
  grad.s11.hpd$signif <- apply(grad.s11.hpd,1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })  
  grad.s12.hpd <- data.frame(gradient_est[[i]]$grad12.est)
  grad.s12.hpd$signif <- apply(grad.s12.hpd,1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })  
  grad.s22.hpd <- data.frame(gradient_est[[i]]$grad22.est)
  grad.s22.hpd$signif <- apply(grad.s22.hpd,1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })  
  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"),)
  polygon(c(meuse.riv[1:88,1],meuse.riv[89:176,1]),
          c(meuse.riv[1:88,2],meuse.riv[89:176,2]),
          col = 'lightblue', border = "red")
  sp_plot(11,"Spectral",cbind(grid.meuse,grad.s1.hpd[,1]),
          shape = meuse.shp,contour.plot = T,
          xlab=latex2exp::TeX("Longitude$\\degree$"),
          ylab = latex2exp::TeX("Latitude$\\degree$"),
          points.plot = T, sig=grad.s1.hpd$signif,add=T)

  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"),)
  polygon(c(meuse.riv[1:88,1],meuse.riv[89:176,1]),
          c(meuse.riv[1:88,2],meuse.riv[89:176,2]),
          col = 'lightblue', border = "red")
  sp_plot(11,"Spectral",cbind(grid.meuse,grad.s2.hpd[,1]), 
          shape=meuse.shp,contour.plot = T,
          xlab=latex2exp::TeX("Longitude$\\degree$"),
          ylab = latex2exp::TeX("Latitude$\\degree$"),
          points.plot = T, sig=grad.s2.hpd$signif,add=T)
  plot(x=1:3,y=1:3,type="n",axes=F, xlab="",ylab="")
  text(x=2,y=2, labels=hvy_metals[i]);plot.new()
  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"),)
  polygon(c(meuse.riv[1:88,1],meuse.riv[89:176,1]),
          c(meuse.riv[1:88,2],meuse.riv[89:176,2]),
          col = 'lightblue', border = "red")
  sp_plot(11,"Spectral",cbind(grid.meuse,grad.s11.hpd[,1]),
          contour.plot = T, shape=meuse.shp,
          xlab=latex2exp::TeX("Longitude$\\degree$"),
          ylab = latex2exp::TeX("Latitude$\\degree$"),
          points.plot = T, sig=grad.s11.hpd$signif,add=T)
  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"),)
  polygon(c(meuse.riv[1:88,1],meuse.riv[89:176,1]),
          c(meuse.riv[1:88,2],meuse.riv[89:176,2]),
          col = 'lightblue', border = "red")
  sp_plot(11,"Spectral",cbind(grid.meuse,grad.s12.hpd[,1]),
          contour.plot = T, shape=meuse.shp,
          xlab=latex2exp::TeX("Longitude$\\degree$"),
          ylab = latex2exp::TeX("Latitude$\\degree$"),
          points.plot = T, sig=grad.s12.hpd$signif,add=T)

  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"),)
  polygon(c(meuse.riv[1:88,1],meuse.riv[89:176,1]),
          c(meuse.riv[1:88,2],meuse.riv[89:176,2]),
          col = 'lightblue', border = "red")
  sp_plot(11,"Spectral",cbind(grid.meuse,grad.s22.hpd[,1]),
          contour.plot = T, shape=meuse.shp,
          xlab=latex2exp::TeX("Longitude$\\degree$"),
          ylab = latex2exp::TeX("Latitude$\\degree$"),
          points.plot = T, sig=grad.s22.hpd$signif,add=T)
}
par(mfcol=c(1,1))
meuse.riv.window <- which(meuse.riv[1:88,1]>178000 & meuse.riv[1:88,1]<182000&meuse.riv[1:88,2]>329000 & meuse.riv[1:88,2]<334000)
meuse.riv.win <- meuse.riv[meuse.riv.window,]

samples <- (nburn+1):niter
gradient_est1 <- lapply(res_mc_sp, function(x){
  spatial_gradient(coords=coords,
                   grid.points = meuse.riv.win/10,
                   post_mcmc = x,
                   cov.type = cov.type,
                   ncores=5,
                   nbatch = 100,
                   return.mcmc = T)
})
# save(gradient_est1, file="~/Desktop/spatiotemporal_gradient/space/NETemp-level4-gradients.RData")
# save(gradient_est, file="~/Desktop/spatiotemporal_gradient/space/NETemp-level4-gradients.RData")
curve <- meuse.riv.win/10
tval <- sapply(1:(nrow(curve)-1), function(x) sqrt(sum((curve[(x+1),]-curve[x,])^2)))
uval <- t(sapply(1:(nrow(curve)-1), function(x) (curve[(x+1),]-curve[x,])))/tval

tu <- list(tval=tval,uval=uval)

cwomb_inf <- list()
for(i in 1:length(gradient_est1)){
  #curve <- subset.points[[i]]
  grad.s1.curve <- gradient_est1[[i]]$grad.s1.mcmc
  grad.s2.curve <- gradient_est1[[i]]$grad.s2.mcmc
  grad.s11.curve <- gradient_est1[[i]]$grad.s11.mcmc
  grad.s12.curve <- gradient_est1[[i]]$grad.s12.mcmc
  grad.s22.curve <- gradient_est1[[i]]$grad.s22.mcmc
  u.perp.mat <- t(apply(tu$uval,1,function(x){
    u.perp <- rev(x)
    u.perp[2] <- -u.perp[2]
    u.perp
  }))
  estval <- estval1 <- c()
  for(j in 1:nrow(grad.s1.curve)){
    estval <- rbind(estval,grad.s1.curve[j,-nrow(curve)]*u.perp.mat[,1]+grad.s2.curve[j,-nrow(curve)]*u.perp.mat[,2])
    estval1 <- rbind(estval1,grad.s11.curve[j,-nrow(curve)]*u.perp.mat[,1]^2+grad.s22.curve[j,-nrow(curve)]*u.perp.mat[,2]^2+2*grad.s12.curve[j,-nrow(curve)]*u.perp.mat[,1]*u.perp.mat[,2])
  }
  slist <- split(1:nrow(grad.s1.curve), ceiling(seq_along(1:nrow(grad.s1.curve))/(length(1:nrow(grad.s1.curve))/100)))
  estval.1 <- do.call(rbind,lapply(slist, function(x) apply(estval[x,],2,median)))
  estval1.1 <- do.call(rbind,lapply(slist, function(x) apply(estval1[x,],2,median)))
  
  # line-segment level
  estval.inf <- data.frame(unlist(round(t(apply(estval.1,2,function(x){
    x <- as.mcmc(x)
    c(median(x),HPDinterval(x))
  })),3)))
  estval.inf$sig <- apply(estval.inf,1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })
  
  estval1.inf <- data.frame(unlist(round(t(apply(estval1.1,2,function(x){
    x <- as.mcmc(x)
    c(median(x),HPDinterval(x))
  })),3)))
  estval1.inf$sig <- apply(estval1.inf,1,function(x){
    if(x[2]>0 & x[3]>0) return (1)
    if(x[2]<0 & x[3]<0) return (-1)
    else return(0)
  })
  
  tval <- tu$tval
  #trval <- 30*pi*cos(3*pi*(curve[-nrow(curve),1]))*u.perp.mat[,1]-30*pi*sin(3*pi*(curve[-nrow(curve),2]))*u.perp.mat[,2]
  #trval1 <- -90*pi^2*sin(3*pi*(curve[-nrow(curve),1]))*u.perp.mat[,1]^2-90*pi^2*cos(3*pi*(curve[-nrow(curve),2]))*u.perp.mat[,2]^2
  w1.c <- c(median(apply(estval,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval,1,function(x) sum(x*tval)/sum(tval)))))
  if(sum(estval.inf$sig==1)>1){
    w1.c.sig.p <- c(median(apply(estval.1[,estval.inf$sig==1],1,function(x) sum(x*tval[estval.inf$sig==1])/sum(tval[estval.inf$sig==1]))), HPDinterval(as.mcmc(apply(estval.1[,estval.inf$sig==1],1,function(x) sum(x*tval[estval.inf$sig==1])/sum(tval[estval.inf$sig==1])))))
  }else{
    w1.c.sig.p <- c(median(sapply(estval.1[,estval.inf$sig==1],function(x) sum(x*tval[estval.inf$sig==1])/sum(tval[estval.inf$sig==1]))), HPDinterval(as.mcmc(sapply(estval.1[,estval.inf$sig==1],function(x) sum(x*tval[estval.inf$sig==1])/sum(tval[estval.inf$sig==1])))))
  } 
  if(sum(estval.inf$sig==-1)>1){
    w1.c.sig.n <- c(median(apply(estval.1[,estval.inf$sig==-1],1,function(x) sum(x*tval[estval.inf$sig==-1])/sum(tval[estval.inf$sig==-1]))), HPDinterval(as.mcmc(apply(estval.1[,estval.inf$sig==-1],1,function(x) sum(x*tval[estval.inf$sig==-1])/sum(tval[estval.inf$sig==-1])))))
  }else{
    w1.c.sig.n <- c(median(sapply(estval.1[,estval.inf$sig==-1],function(x) sum(x*tval[estval.inf$sig==-1])/sum(tval[estval.inf$sig==-1]))), HPDinterval(as.mcmc(sapply(estval.1[,estval.inf$sig==-1],function(x) sum(x*tval[estval.inf$sig==-1])/sum(tval[estval.inf$sig==-1])))))
  } 
  #tw1.c <- sum(trval*tval)/sum(tval)
  w2.c <- c(median(apply(estval1,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval1,1,function(x) sum(x*tval)/sum(tval)))))
  if(sum(estval1.inf$sig==1)>1){
    w2.c.sig.p <- c(median(apply(estval1.1[,estval1.inf$sig==1],1,function(x) sum(x*tval[estval1.inf$sig==1])/sum(tval[estval1.inf$sig==1]))), HPDinterval(as.mcmc(apply(estval1.1[,estval1.inf$sig==1],1,function(x) sum(x*tval[estval1.inf$sig==1])/sum(tval[estval1.inf$sig==1])))))
  }else{
    w2.c.sig.p <- c(median(sapply(estval1.1[,estval1.inf$sig==1],function(x) sum(x*tval[estval1.inf$sig==1])/sum(tval[estval1.inf$sig==1]))), HPDinterval(as.mcmc(sapply(estval1.1[,estval1.inf$sig==1],function(x) sum(x*tval[estval1.inf$sig==1])/sum(tval[estval1.inf$sig==1])))))
  } 
  if(sum(estval1.inf$sig==-1)>1){
    w2.c.sig.n <- c(median(apply(estval1.1[,estval1.inf$sig==-1],1,function(x) sum(x*tval[estval1.inf$sig==-1])/sum(tval[estval1.inf$sig==-1]))), HPDinterval(as.mcmc(apply(estval1.1[,estval1.inf$sig==-1],1,function(x) sum(x*tval[estval1.inf$sig==-1])/sum(tval[estval1.inf$sig==-1])))))
  }else{
    w2.c.sig.n <- c(median(sapply(estval1.1[,estval1.inf$sig==-1],function(x) sum(x*tval[estval1.inf$sig==-1])/sum(tval[estval1.inf$sig==-1]))), HPDinterval(as.mcmc(sapply(estval1.1[,estval1.inf$sig==-1],function(x) sum(x*tval[estval1.inf$sig==-1])/sum(tval[estval1.inf$sig==-1])))))
  }
  #tw2.c <- sum(trval1*tval)/sum(tval)
  cwomb_inf[[i]] <- list(estval.inf,w1.c,estval1.inf,w2.c,
                         w1.c.sig.p,w1.c.sig.n,
                         w2.c.sig.p,w2.c.sig.n)
}
xtable::xtable(cbind(do.call(rbind,cwomb_inf[[1]][5:8]),
      do.call(rbind,cwomb_inf[[2]][5:8]),
      do.call(rbind,cwomb_inf[[3]][5:8]),
      do.call(rbind,cwomb_inf[[4]][5:8])),auto=T)

inf <- c()
for(i in 1:length(gradient_est1)){
  inf <- rbind(inf,c(cwomb_inf[[i]][[2]],
                     cwomb_inf[[i]][[4]],
                     cwomb_inf[[i]][[5]],
                     cwomb_inf[[i]][[6]],
                     cwomb_inf[[i]][[7]],
                     cwomb_inf[[i]][[8]]))
}
round(inf[,-c(1:6)],3)
for(i in 1:length(gradient_est1)){
  mat <- matrix(c(1:3), nr=1,nc=3, byrow=T)
  layout(mat,
         widths = c(rep(5,3)),
         heights = c(rep(3,3)))
  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"))
  
  sp_plot(col.seq.length = 11,
          col.text = "Spectral",
          data_frame = cbind(meuse@coords,z[[i]][,"median"]),
          shape=meuse.shp,
          contour.plot = T,
          xlab=latex2exp::TeX("Longitude$\\degree$"),
          ylab = latex2exp::TeX("Latitude$\\degree$"),
          points.plot=T,
          raster.surf = F,
          sig = z[[i]]$sig,
          legend = F,
          add=T)
  
  lines(curve*10,lwd=3.5,type="l")
  
  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"))
  
  sp_plot(col.seq.length = 11,
          col.text = "Spectral",
          data_frame = cbind(meuse@coords,z[[i]][,"median"]),
          shape=meuse.shp,
          contour.plot = T,
          xlab=latex2exp::TeX("Longitude$\\degree$"),
          ylab = latex2exp::TeX("Latitude$\\degree$"),
          points.plot=T,
          raster.surf = F,
          sig = z[[i]]$sig,
          legend = F,
          add=T)
  
  
  
  estval.inf <- cwomb_inf[[i]][[1]]
  col.segment <- sapply(estval.inf$sig, function(x){
    if(x==0) "white"
    else if(x==1) "green"
    else if(x==-1) "cyan"
  })
  
  for(j in 1:(nrow(curve)-1)){
    if(estval.inf$sig[j] %in%c(-1,1)){
      lines(rbind(curve[j,]*10,curve[(j+1),]*10), lwd=4.5)
    } 
    # else  lines(rbind(curve[j,],curve[(j+1),]),col="white", lwd=0.5)
  }
  for(j in 1:(nrow(curve)-1)){
    if(estval.inf$sig[j] %in%c(-1,1)){
      lines(rbind(curve[j,]*10,curve[(j+1),]*10),col=col.segment[j], lwd=0.5)
    }
  }
  
  plot(meuse@coords,type="n", xlim=c(178000,182000), ylim=c(329000,334000),
       xlab=latex2exp::TeX("UTM$_x$"),
       ylab=latex2exp::TeX("UTM$_y$"))
  
  sp_plot(col.seq.length = 11,
          col.text = "Spectral",
          data_frame = cbind(meuse@coords,z[[i]][,"median"]),
          shape=meuse.shp,
          contour.plot = T,
          xlab=latex2exp::TeX("Longitude$\\degree$"),
          ylab = latex2exp::TeX("Latitude$\\degree$"),
          points.plot=T,
          raster.surf = F,
          sig = z[[i]]$sig,
          legend = F,
          add=T)
  
  
  estval.inf <- cwomb_inf[[i]][[3]]
  col.segment <- sapply(estval.inf$sig, function(x){
    if(x==0) "white"
    else if(x==1) "green"
    else if(x==-1) "cyan"
  })
  
  for(j in 1:(nrow(curve)-1)){
    if(estval.inf$sig[j] %in%c(-1,1)){
      lines(rbind(curve[j,]*10,curve[(j+1),]*10), lwd=4.5)
    } 
    # else  lines(rbind(curve[j,],curve[(j+1),]),col="white", lwd=0.5)
  }
  for(j in 1:(nrow(curve)-1)){
    if(estval.inf$sig[j] %in%c(-1,1)){
      lines(rbind(curve[j,]*10,curve[(j+1),]*10),col=col.segment[j], lwd=0.5)
    }
  }
}


