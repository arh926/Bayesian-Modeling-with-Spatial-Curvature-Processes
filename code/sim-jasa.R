################################
# Simulation for JASA Revision #
################################
require(spWombling)
require(coda)
require(parallel)
###########################
# Synthetic Spatial Data  #
###########################
# set.seed(123456)
N = 100
nrep = 20
tau = 1
# synthetic location
min_sp = 0
max_sp = 1

# synthetic pattern 
sig2 = c(10, 100, 200)
phis = 3

niter = 1e4
nburn = niter/2
report = 1e2
X = matrix(1,nr = N) # intercept
XtX = crossprod(X,X) 
D = 1:N
DtD = diag(N)


# Gradient Analysis
grid.points = as.matrix(expand.grid(seq(min_sp, max_sp, length.out = 21)[-c(1,21)],
                                     seq(min_sp, max_sp, length.out = 21)[-c(1,21)]),
                         nc = 2)

colnames(grid.points) = NULL
samples = (nburn+1):niter

sim.jasa.results = list()
for(i in 1:length(sig2)){
  results = list()
  for(j in 1:nrep){
    coords = matrix(runif(2*N, min = min_sp, max = max_sp), nr = N, nc = 2) 
    Delta = as.matrix(dist(coords))
    cov.true = sig2[i]*(1+phis*sqrt(5)*Delta+phis^2*5*Delta^2/3)*exp(-phis*sqrt(5)*Delta) + tau^2*diag(N)
    
    y = MASS::mvrnorm(1, rep(0,N), cov.true)
    
    y = y # response
    # cannot be randomly started
    z_init = rep(0,N)
    lower_phis = 3/max(dist(coords))
    upper_phis = 30
    shape_sigma = 2
    scale_sigma = 1
    shape_tau = 2
    scale_tau = 1
    mean_beta = rep(0,ncol(X))
    prec_beta = 1e-6*diag(ncol(X))
    Delta = as.matrix(dist(coords)) # instead of coordinates supply dist
    steps_init = 1
    cov.type = "matern2"
    
    phis_init = runif(1, lower_phis, upper_phis)
    sigma2_init = 1/rgamma(1,shape = shape_sigma, rate = 1/scale_sigma)
    tau2_init = 1/rgamma(1,shape = shape_tau, rate = 1/scale_tau)
    beta_init = MASS::mvrnorm(1, mean_beta, 1/prec_beta)
    nu_init = 2.5
    nu_est = F
    mc_sp = hlmBayes_sp(coords = coords,
                         y = y,
                         X = X,
                         z_init = z_init,
                         D = D,DtD = DtD,
                         phis_init = phis_init,
                         lower_phis = lower_phis,
                         upper_phis = upper_phis,
                         sigma2_init = sigma2_init,
                         shape_sigma = shape_sigma,
                         scale_sigma = scale_sigma,
                         tau2_init = tau2_init,
                         shape_tau = shape_tau,
                         scale_tau = scale_tau,
                         beta_init = beta_init,
                         mean_beta = mean_beta,
                         prec_beta = prec_beta,
                         steps_init = steps_init,
                         niter = niter,
                         report = report,
                         cov.type = cov.type,
                         nu_est = nu_est,
                         nu_init = nu_init,
                         lower_nu = lower_nu,
                         upper_nu = upper_nu,
                         stepnu_init = stepnu_init, verbose = F)
    
    model_summary = hlm_summary(chain = mc_sp, nburn = nburn, thin = 1)
    coef = model_summary$summary.pars; round(coef,4)
    
    z = model_summary$summary.latent
    z$sig = apply(z,1,function(x){
      if(x[2]>0 & x[3]>0) return (1)
      if(x[2]<0 & x[3]<0) return (-1)
      else return(0)
    })
    
    
    
    gradient_est = spatial_gradient(coords = coords,
                                     grid.points = grid.points,
                                     samples = samples,
                                     post_mcmc = mc_sp,
                                     cov.type = cov.type,
                                     ncores = 6,
                                     nbatch = 200,
                                     return.mcmc = T)
    results[[j]] = list(coef = coef,
                        z = z,
                        grad.est = gradient_est)
    cat("sig2:\t", sig2[i], "replication:\t", j,"\n")
  }
  sim.jasa.results[[i]] = results
}
# Compute lengths of HPDs
lengths_hpd = lapply(sim.jasa.results, function(x){
  lapply(x, function(x.1){
    lapply(x.1$grad.est[1:5], function(x.2){
      x.2[,3]-x.2[,2]
    })
  })
})
# Summarize
sim_results = lapply(lengths_hpd, function(x){
  res = do.call(rbind, lapply(x, function(x.1){
    c(median(x.1$grad1.est),
      median(x.1$grad2.est),
      median(x.1$grad11.est),
      median(x.1$grad12.est),
      median(x.1$grad22.est))
  }))
  colnames(res) = c("g.1","g.2", "g.11", "g.12","g.22")
  res
})

par(mfcol=c(1,3))
boxplot(sim_results[[1]], ylim=c(0,1100), ylab = "median width of HPD (in units)", xlab = "sig2 = 10")
boxplot(sim_results[[2]], ylim=c(0,1100), ylab = "", xlab = "sig2 = 100")
boxplot(sim_results[[3]], ylim=c(0,1100), ylab = "", xlab = "sig2 = 200")

# signif_hpd = lapply(gradient_est.vals, function(x){
#   apply(x, 1, function(x.r){
#     if(x.r[2]>0 & x.r[3]>0) return (1)
#     if(x.r[2]<0 & x.r[3]<0) return (-1)
#     else return(0)
#   })
# })
mat = matrix(c(1:12), nr=2,nc=6, byrow=T)
layout(mat,
       widths = rep(c(5,3),3),
       heights = c(3,3))

sp_plot(11,"Spectral",cbind(coords,y), contour.plot=T, points.plot = T, sig = z$sig)
sp_plot(11,"Spectral",cbind(grid.points, gradient_est$grad1.est[,1]), contour.plot=T, points.plot = T, sig = signif_hpd$grad1.est)
sp_plot(11,"Spectral",cbind(grid.points, gradient_est$grad2.est[,1]), contour.plot=T, points.plot = T, sig = signif_hpd$grad2.est)
sp_plot(11,"Spectral",cbind(grid.points, gradient_est$grad11.est[,1]), contour.plot=T, points.plot = T, sig = signif_hpd$grad11.est)
sp_plot(11,"Spectral",cbind(grid.points, gradient_est$grad12.est[,1]), contour.plot=T, points.plot = T, sig = signif_hpd$grad12.est)
sp_plot(11,"Spectral",cbind(grid.points, gradient_est$grad22.est[,1]), contour.plot=T, points.plot = T, sig = signif_hpd$grad22.est)
