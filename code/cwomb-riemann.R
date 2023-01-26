###################################
# Wombling : Riemann-Sum Approach #
###################################
samples <- (nburn+1):niter
gradient_est1 <- list()
for(i in 1:length(subset.points)){
  grid.points <- subset.points[[i]]
  gradient_est1[[i]] <- spatial_gradient(coords=coords,
                                         chain = mc_sp,
                                         grid.points = grid.points,#*100, for boston.shp
                                         cov.type = cov.type,
                                         ncores=5,
                                         nbatch = 100,
                                         nburn = nburn,
                                         niter= niter)
  cat("Boundary:",i,"\n")
}
# save(gradient_est, file="~/Desktop/spatiotemporal_gradient/space/NETemp-level4-gradients.RData")
# load("~/Desktop/spatiotemporal_gradient/space/NETemp-level1-gradients.RData")
tu <- lapply(subset.points, function(curve){
  tval <- sapply(1:(nrow(curve)-1), function(x) sqrt(sum((curve[(x+1),]-curve[x,])^2)))
  uval <- t(sapply(1:(nrow(curve)-1), function(x) (curve[(x+1),]-curve[x,])))/tval
  list(tval=tval,uval=uval)
})

cwomb_inf <- list()
for(i in 1:length(gradient_est1)){
  curve <- subset.points[[i]]
  grad.s1.curve <- gradient_est1[[i]]$grad.s1.mcmc
  grad.s2.curve <- gradient_est1[[i]]$grad.s2.mcmc
  grad.s11.curve <- gradient_est1[[i]]$grad.s11.mcmc
  grad.s12.curve <- gradient_est1[[i]]$grad.s12.mcmc
  grad.s22.curve <- gradient_est1[[i]]$grad.s22.mcmc
  u.perp.mat <- t(apply(tu[[i]]$uval,1,function(x){
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
  
  tval <- tu[[i]]$tval
  #trval <- 30*pi*cos(3*pi*(curve[-nrow(curve),1]))*u.perp.mat[,1]-30*pi*sin(3*pi*(curve[-nrow(curve),2]))*u.perp.mat[,2]
  #trval1 <- -90*pi^2*sin(3*pi*(curve[-nrow(curve),1]))*u.perp.mat[,1]^2-90*pi^2*cos(3*pi*(curve[-nrow(curve),2]))*u.perp.mat[,2]^2
  w1.c <- c(median(apply(estval,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval,1,function(x) sum(x*tval)/sum(tval)))))
  #tw1.c <- sum(trval*tval)/sum(tval)
  w2.c <- c(median(apply(estval1,1,function(x) sum(x*tval)/sum(tval))), HPDinterval(as.mcmc(apply(estval1,1,function(x) sum(x*tval)/sum(tval)))))
  #tw2.c <- sum(trval1*tval)/sum(tval)
  cwomb_inf[[i]] <- list(estval.inf,w1.c,estval1.inf,w2.c)
}

inf <- c()
for(i in 1:length(gradient_est1)){
  inf <- rbind(inf,c(cwomb_inf[[i]][[2]],
                     cwomb_inf[[i]][[4]]))
}
round(inf,3)

mat <- matrix(c(1:3), nr=1,nc=3, byrow=T)
layout(mat,
       widths = c(rep(5,3)),
       heights = c(rep(3,3)))
sp_plot(11,"Spectral",cbind(coords,z[,1]),#/100
        shape=NETemp.shp,#boston.shp
        contour.plot = T,
        xlab=latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        points.plot=T,
        raster.surf = F,
        sig = z$sig,
        legend = F)
lapply(subset.points, function(x){
  lines(x,lwd=3.5)
  return(0)})
# lines(subset.points.1@lines[[1]]@Lines[[1]]@coords,lwd=3.5)
# lines(subset.points.2@lines[[1]]@Lines[[1]]@coords,lwd=3.5)
#lines(subset.points.2@lines[[1]]@Lines[[2]]@coords,lwd=3.5)
#lines(subset.points.2@lines[[1]]@Lines[[5]]@coords,lwd=3.5)

sp_plot(11,"Spectral",cbind(coords,z[,1]),
        shape=NETemp.shp,
        contour.plot = T,
        xlab=latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        legend = F)
lapply(subset.points, function(x){
  lines(x,lwd=1.5, lty="dashed")
  return(0)})
# lines(subset.points.1@lines[[1]]@Lines[[1]]@coords,lwd=1.5, lty="dashed")
# lines(subset.points.2@lines[[1]]@Lines[[1]]@coords,lwd=1.5, lty="dashed")
# lines(subset.points.2@lines[[1]]@Lines[[2]]@coords,lwd=1.5, lty="dashed")
# lines(subset.points.2@lines[[1]]@Lines[[5]]@coords,lwd=1.5, lty="dashed")

for(i in 1:length(gradient_est1)){
  curve <- subset.points[[i]]
  estval.inf <- cwomb_inf[[i]][[1]]
  col.segment <- sapply(estval.inf$sig, function(x){
    if(x==0) "white"
    else if(x==1) "green"
    else if(x==-1) "cyan"
  })
  
  for(j in 1:(nrow(curve)-1)){
    if(estval.inf$sig[j] %in%c(-1,1)){
      lines(rbind(curve[j,],curve[(j+1),]), lwd=4.5)
    } 
    # else  lines(rbind(curve[j,],curve[(j+1),]),col="white", lwd=0.5)
  }
  for(j in 1:(nrow(curve)-1)){
    if(estval.inf$sig[j] %in%c(-1,1)){
      lines(rbind(curve[j,],curve[(j+1),]),col=col.segment[j], lwd=0.5)
    }
  }
}
sp_plot(11,"Spectral",cbind(coords,z[,1]),
        shape=NETemp.shp,
        contour.plot = T,
        xlab=latex2exp::TeX("Longitude$\\degree$"),
        ylab = latex2exp::TeX("Latitude$\\degree$"),
        legend = F)
lapply(subset.points, function(x){
  lines(x,lwd=1.5, lty="dashed")
  return(0)})
# lines(subset.points.1@lines[[1]]@Lines[[1]]@coords,lwd=1.5, lty="dashed")
# lines(subset.points.2@lines[[1]]@Lines[[1]]@coords,lwd=1.5, lty="dashed")
# lines(subset.points.2@lines[[1]]@Lines[[2]]@coords,lwd=1.5, lty="dashed")
# lines(subset.points.2@lines[[1]]@Lines[[5]]@coords,lwd=1.5, lty="dashed")

for(i in 1:length(gradient_est1)){
  curve <- subset.points[[i]]
  estval.inf <- cwomb_inf[[i]][[3]]
  col.segment <- sapply(estval.inf$sig, function(x){
    if(x==0) "white"
    else if(x==1) "green"
    else if(x==-1) "cyan"
  })
  
  for(j in 1:(nrow(curve)-1)){
    if(estval.inf$sig[j] %in%c(-1,1)){
      lines(rbind(curve[j,],curve[(j+1),]), lwd=4.5)
    } 
    # else  lines(rbind(curve[j,],curve[(j+1),]),col="white", lwd=0.5)
  }
  for(j in 1:(nrow(curve)-1)){
    if(estval.inf$sig[j] %in%c(-1,1)){
      lines(rbind(curve[j,],curve[(j+1),]),col=col.segment[j], lwd=0.5)
    }
  }
}


