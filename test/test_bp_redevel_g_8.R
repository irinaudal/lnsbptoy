## Small simulation example to demonstrate the bias induced due to pi approximation.


##### Very useful function: 
# gives insight about ability to estimate theta without any missing data. Larget theta will inherintely imply larger sd.
#
# R

#xmin <- 1.39e-13 ; k <- 2.5 ; nsamp <- 100 ; sqrt(nsamp)/sum(log(rpareto(nsamp, k=k, x_min=xmin)/xmin))

# Model: 
# Y_i ~ Poisson(lambda = S_i*E/gamma)
# S_i ~ Pareto(theta,Smin,bp)
# theta ~ gamma(a,b)
# N ~ negbinom(alpha,beta)
# Smin ~ gamma(am,bm)
# bp ~ fcn of eta|Smin
# eta | Smin ~ N(mu,C), where eta = log(bp-Smin)
# I_i ~ Bin(1,g(lambda_i))
# g(lambda_i = S_i*E/gamma) = P(observe i-th source) = const
# pi(theta,Smin,bp) = \int g(lambda_i) * p(S_i | theta,Smin,bp)  dS_i 
#    \approx  g+e,   where e ~ N(0, sigma)
# OR
#    \approx via MC integration <- TODO: implement this later after the code validates
# sigma = 0.001, say is error that is incurred while approximating pi


# Start from clean workspace:
rm(list=ls())

print(Sys.info()["sysname"])

# Setup correct paths, R history etc...
if (Sys.info()["sysname"]=="Darwin"){
  # Note: main.dir should NOT have a trailing /
  main.dir   <- "~/Dropbox/Log_N_Log_S/R_Code"
  setwd(main.dir)
  libsrc.dir <- "."
  library(lnsbptoy)
  
} else if (Sys.info()["sysname"]=="Linux"){
  libdir <- "/home/isudal/R/x86_64-pc-linux-gnu-library/3.0/"
  library(lnsbptoy,lib.loc=libdir)
  main.dir   <- "/home/isudal/bplns" #lns_data/"
  setwd(main.dir)
  libsrc.dir <- "."
}

## Handle batch job arguments:
args <- commandArgs(TRUE)  # 1-indexed version is used now.
#args <- 3
cat(paste("Command-line arguments:\n"))
print(args)

###### define another function
# findInterval written in R as a binary search
findInterval2 <- function(x,v,all.inside=TRUE) {
  n <- length(v)
  if (x<v[1]) {
    if (all.inside){
      return(1)
    } else {
      return(0)
    }
  }
  if (x>=v[n]) {
    return (n)
  }
  i <- 1
  k <- n
  while({j = (k-i) %/% 2 + i; !(v[j] <= x && x < v[j+1])}) {
    if (x < v[j]){
      k <- j
    } else {
      i <- j+1
    }
  }
  return(j)
}
library(compiler)
#findInterval2 compilated with cmpfun
findInterval2 <- cmpfun(findInterval2)



###################
sim_start <- 1000
if (length(args)==0){
  sinkit <- TRUE #FALSE
  sim_num <- sim_start + 1
} else {
  sinkit <- TRUE
  sim_num <- sim_start + as.numeric(args[1])
}
set.seed(762*sim_num + 1330931 + 10231*92) 
###################

# Set simulation parameters
niter  <- 110000
burnin <-  10000
tune.iter <- 100
stop.tune <- 5000
verbose <- FALSE

print.every <- 10000
save.every  <- 1
save.progress <- FALSE
store_logPost <- FALSE
do.profile <- FALSE
save_all_draws <- ifelse(sim_num<1002,TRUE,FALSE)
do.mcmc.plots  <- ifelse(sim_num<1011,TRUE,FALSE)
precompute.pi.theta <- TRUE
    
model <- "bp"

fixed.theta <- FALSE #c(.6,1.2)   # Vector of values / FALSE
fixed.N     <- FALSE #150 #FALSE  # Value / FALSE
fixed.Smin  <- FALSE #1*10^-13
fixed.bp    <- FALSE #4*10^-13 #NULL
fixed.S.obs <- FALSE
fixed.S.mis <- FALSE

outer.dir.extra     <- paste("/bp_lns_toy_example_1bp_g_8_prec", sep="")
output.dir.extra    <- paste("/dataset_",sim_num,sep="")
output.dir.outer    <- paste(main.dir,outer.dir.extra,sep="")
output.dir          <- paste(output.dir.outer,output.dir.extra,sep="")
output.dir.wo.main  <- paste(outer.dir.extra,output.dir.extra,sep="")

sink.file           <- file.path(output.dir,"foo.txt")
profile_file        <- paste(output.dir,"/basicProfile.Rprof",sep="")
save.res.file       <- paste(output.dir,"/res_image.RData",sep="")
save.image.file     <- paste(output.dir,"/store_image.RData",sep="")

output.lnsplot     <- paste(output.dir,"/logN_logS_",model,".jpg",sep="")
output.mcmcplot    <- paste(output.dir,"/mcmc_draws_",model,"_%02d.jpg",sep="")
output.logPost     <- paste(output.dir,"/log_posterior.jpg",sep="")
output.priorsplot  <- paste(output.dir,"/model_priors.pdf",sep="")

if (!file.exists(output.dir.outer)){   #set directory - done after simulate()
  dir.create(output.dir.outer, recursive=TRUE)
  #   setwd(output.dir)
}
if (!file.exists(output.dir)){   #set directory - done after simulate()
  dir.create(output.dir, recursive=TRUE)
  #   setwd(output.dir)
}

# Use Normal additive error in approximation of the term pi(S)
sigma.vec <- c(0, 0.00001, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2)
subset <- 1 #c(1:9)
sigma.vec <- sigma.vec[subset]
length.jobs <- length(sigma.vec)

# Incompleteness
"g.func" <- function(lambda) {
  return(rep(0.8,length(lambda)))   # constant incompleteness
  #return(pnorm(lambda,mean=0,sd=40))   #mean=-25,sd=50))  # 1 minus exponential decay
}
nsamples  <- 10000
E <- 400000 #190000
gamma <- 1.6*(10^(-9))

# Precomputing pi parameters
length.theta <- 10 #80
length.Smin  <- 10 #300
length.bp    <- 10 #80

# Tuning parameters
v.th  <- 0.2 #proposal sd for theta
v.sm  <- 1.0 #proposal sd for Smin
v.bp  <- 0.5 #proposal sd for eta (bp)
v.so  <- 2.0*(10^(-13)) # proposal sd for S, flux

# Parameters for Normal priors for break-point(s):
mu <- -29.0 # -29.8 # c(-30) # dim: eta(bp)=m-1
C  <- 0.1  #0.2 # c(0.5) 
##############################
# Parameters for Gamma prior(s) for theta(s):
a <- c(68,185)  #c(72,160) #100  #c(10,30)   #c(12,10) #c(25,30)
b <- c(120,190)  #c(120,160)  #80   #c(20,30)   #c(18,10) #c(40,30)
##############################
# Parameters for Gamma prior(s) for minimum threshold parameter Smin(s):
Smin.prior.mean <- 1.0*10^(-13)    # Mean = am * bm
Smin.prior.var  <- (2.0*10^(-14))^2  #(2*10^(-14))^2 #(4.5*10^(-14))^2  # Var  = am * bm^2
am  <- Smin.prior.mean^2 / Smin.prior.var   # shape
bm  <- am / Smin.prior.mean                 # rate
##############################
# parameters for neg binom prior of N
#    Mean = alpha / beta
#    Var  = alpha * (1+beta) / beta^2
#    => beta   = prior.mean / (prior.var - prior.mean)
#    => alpha = beta * prior.mean  
N.prior.mean <- 400 #80 ##300
N.prior.var  <- 50^2 #50^2 #(1*100)^2
beta  <- N.prior.mean / (N.prior.var - N.prior.mean)
alpha <- beta * N.prior.mean
###################
# Which posterior quantiles should be stored:
q <- seq(0.01,0.99,by=0.01) #c(0,.025,.25,.5,.75,.975,1)
# Which posterior quantiles of Smis should be validated?
Smis.q.check <- q
max.lines           <- 1000
max.mcmc.plots      <- ifelse(sim_num<1005,"all",10)      # "all" or value
###################


########################################################################
####################### RUN ############################################
if (sinkit) {
  cat("Sinking to...\n")
  print(sink.file)
  sink(sink.file)
}

res <- simulate.mix(a=a,b=b,alpha=alpha,beta=beta,C=C,mu=mu,am=am,bm=bm,gamma=gamma,E=E,
                    fixed.theta=fixed.theta, fixed.N=fixed.N, fixed.Smin=fixed.Smin, fixed.bp=fixed.bp, 
                    g=g.func, verbose=FALSE) 

if (fixed.N)
  fixed.N <- res$par$N
if (fixed.Smin)
  fixed.Smin <- res$par$tau.[1]
if (any(fixed.bp != FALSE)){
  idx <- (fixed.bp != FALSE)
  if (any(fixed.bp == TRUE)){
    idx2 <- which(fixed.bp==TRUE)
    fixed.bp[idx2] <- res$par$tau.[idx2]
  } 
}  
if (any(fixed.theta != FALSE)){
  idx <- (fixed.theta != FALSE)
  if (any(fixed.theta == TRUE)){
    idx2 <- which(fixed.theta==TRUE)
    fixed.theta[idx2] <- res$par$theta[idx2]
  }
}  

cat("Data has been generated:\n")
print(res$obs)
print(res$par)
cat("All done. :)\n\n")

# State menbership proportions
if (length(a)>1){
  ct1 <- sum(res$par$S.obs < res$par$tau.[2])
  ct2 <- res$obs$n - ct1
  pct1 <- ct1/(ct1+ct2)
  pct2 <- ct2/(ct1+ct2)
  cat("The membership (observed-source) source counts are:\n")
  print(c(ct1,ct2))
  cat("The membership (observed-source) proportions are:\n")
  print(c(pct1,pct2))
  cat("\n")
}

# Store true parameters as vector
true.par <- unlist(res$par)
id.N  <- grep("N",names(true.par))
id.th <- grep("theta",names(true.par))
id.sm <- grep("tau.1",names(true.par))
id.bp <- grep("tau",names(true.par))[-1]
true.S.mis <- res$mis$S.mis


for (job_num in 1:length.jobs ) {
  
  pi <- NULL
  if (precompute.pi.theta){
    
    pi.file.location    <- paste(output.dir.outer,"/pi_theta_",job_num,"_8NormErr.RData",sep="")
    
    cat("pi.theta file location specified as:\n")
    print(pi.file.location)
    cat("Does the file already exist?\n")
    if (file.exists(pi.file.location)){
      # Case 1/4: precompute==TRUE, previous file exists, Smin==fixed/random.
      cat("Yes. Loading from file...\n")
      precompute.pi.theta <- FALSE
      # Note: the following default objects will be re-loaded with pi.theta: 
      #       pi,E,fixed.Smin,fixed.bp,gamma,g.func
      load(pi.file.location) 
      
    } else {
      cat("No. Need to pre-compute pi(theta).\n")
      cat("Note: pi(theta) will be saved, and can be used for future analyses...\n")
    }
    
    if (precompute.pi.theta){
      cat("Pre-computing pi(theta)...\n")
      
      # Split computation into segments of Smin
      Smin.prior.mean <-  1.0*10^(-13)    # Mean = am * bm
      Smin.prior.var  <- (8.0*10^(-14))^2  # Var  = am * bm^2
      am2  <- Smin.prior.mean^2 / Smin.prior.var   # shape
      bm2  <- am2 / Smin.prior.mean                 # rate
      
      pi <- pi.theta.compute(theta.grid=NULL, Smin.grid=NULL, bp.grid=NULL, gamma=gamma, g=g.func, 
                             E=E, a=a,b=b,C=C,mu=mu,am=am2,bm=bm2, 
                             length.theta=length.theta, length.Smin=length.Smin, length.bp=length.bp,
                             fixed.Smin=fixed.Smin, fixed.bp=fixed.bp,
                             theta.lo=c(0.3,0.7), theta.hi=c(0.9,1.3), 
                             Smin.hi=4*10^-13,
                             grid.eps=10^-2,
                             mc.integral=TRUE, nsamples=nsamples,
                             verbose=TRUE,
                             just.return.grid=FALSE)
  
      save(pi,E,fixed.Smin,fixed.bp,gamma,g.func, file=pi.file.location)
    }
  } # END precompute.pi.theta
  
  
  sigma <- sigma.vec[job_num]
  print(sigma)
  
  if(do.profile){
    Rprof(profile_file)
  }
  
  startTime <- proc.time()  # Get time for runtime measurement
  
  resmc <- analyze.mix(Y=res$obs$Y.obs.tot, niter=niter, burnin=burnin, v.so=v.so, v.th=v.th, v.sm=v.sm, v.bp=v.bp,
                       a=a, b=b, alpha=alpha, beta=beta, C=C,mu=mu, am=am,bm=bm, gamma=gamma, E=E,
                       nsamples=nsamples, g=g.func, pi=pi,
                       fixed.N=fixed.N, fixed.theta=fixed.theta, 
                       fixed.S.obs=fixed.S.obs, fixed.S.mis=fixed.S.mis, 
                       fixed.Smin=fixed.Smin, fixed.bp=fixed.bp,
                       print.every=print.every, 
                       save.every=save.every, save.progress=save.progress, save.progress.dir=output.dir,
                       store_logPost=store_logPost, sigma=sigma,
                       tune.iter=tune.iter, stop.tune=stop.tune, verbose=verbose)
  
  endMCMCTime <- proc.time()
  
  cat("All done with MCMC. :)\n\n")
  cat(sprintf("MCMC computation time (%.1f seconds)\n", round((endMCMCTime-startTime)[3], digits=1))) 
  
  if(do.profile){
    Rprof(NULL)
    print(summaryRprof(profile_file))
  }
  
  # Store statistics
  mcmc.stats <- summary(resmc$draws, quantiles=q)
  # Get posterior quantiles (credible intervals) of parameters, according to percentiles 'q'
  post.q <- summary(resmc$draws, quantiles=q)$quantiles
  
  input.pars <- resmc$inputs
  loglik     <- resmc$loglik
  eff.size   <- resmc$eff.size
  
  cat("Computing posterior summary statistics for Smis...\n")
  # First, check if we have any Smis samples
  if(length(unlist(resmc$draws.S.mis))==0 || length(true.S.mis)==0) {
    warning("There were no missing sources sampled in this simulation.")
    post.q.S.mis <- NULL
    true.par.S.mis <- NULL
  } else {
    # Set up convenient quantities for computing posterior statistics for S.mis:
    stats.S.mis <- matrix(numeric(0),nrow=niter-burnin,ncol=2+length(q))
    for (j in 1:length(resmc$draws.S.mis)){
      S.mis.vec <- resmc$draws.S.mis[[j]]
      stats.S.mis[j,] <- c(mean(S.mis.vec),sd(S.mis.vec), quantile(S.mis.vec,probs=Smis.q.check))    
    }
    # Get posterior quantiles (credible intervals) of S.mis statistics, according to percentiles 'q'
    post.q.S.mis <- t(apply(stats.S.mis,2, quantile, probs=q, na.rm=TRUE))
    rownames(post.q.S.mis) <- c("mean","sd",paste("q_",Smis.q.check,sep=""))
    # Compute summary statistics for the true Smis values:
    true.par.S.mis <- c(mean(true.S.mis), sd(true.S.mis), quantile(true.S.mis, probs=Smis.q.check))
    names(true.par.S.mis) <- c("mean","sd",paste("q_",Smis.q.check,sep=""))
  }
  
  cat("Names of resmc:"); print(names(resmc))
  cat("Names of draws:"); print(colnames(resmc$draws))
  cat("prop.accept.S=\n");     print(resmc$prop.accept$S.prop.accept)
  cat("prop.accept.theta=\n"); print(resmc$prop.accept$theta.prop.accept)
  cat("prop.accept.Smin=\n");  print(resmc$prop.accept$Smin.prop.accept)
  cat("prop.accept.bp=\n");    print(resmc$prop.accept$bp.prop.accept)
  cat("used sd,  v=\n"); print(resmc$prop.accept$v)
  cat("Eff.sample.size =\n"); print(eff.size)
  
  if(store_logPost){
    cat("Evaluating log(posterior) plot...\n")
    
    Smin <- res$par$tau.[1]
    bp   <- res$par$tau.[-1]
    theta.t <- res$par$theta
    verbose2 <- FALSE
    N.t  <- res$par$N
    n    <- res$obs$n
    S.obs.t   <- res$par$S.obs
    Y.obs.tot <- res$obs$Y.obs.tot
    pi.value <- pi.theta.get(theta=theta.t, Smin=Smin, bp=bp, gamma=gamma, E=E,
                             g=g.func, nsamples=nsamples, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
    pi.const <- ifelse(N.t==n,0,(N.t-n)*log(1.0-pi.value))
    p.N   <- lgamma(N.t+alpha)-lgamma(N.t-n+1) + pi.const - N.t*log(1+beta)
    p.th  <- sum(dgamma(x=theta.t, shape=a, rate=b, log=TRUE)) 
    p.Sm  <- dgamma(x=Smin, shape=am, rate=bm, log=TRUE)
    p.bp  <- dnorm(x=log(bp-Smin), mean=mu, sd=C, log=TRUE)
    lambda <- S.obs.t*E/gamma
    gp    <- log(g.func(lambda=lambda))
    p.S   <- dbrokenpareto(S.obs.t, x_min=Smin, k=theta.t, bp=bp, log=TRUE)
    p.Yt  <- dpois(x=Y.obs.tot, lambda=lambda, log=TRUE)
    
    terms <- list(p.N=p.N,p.th=p.th,p.Sm=p.Sm,p.bp=p.bp,p.gp=gp,p.S=p.S,p.Yt=p.Yt)
    true.logPost <- sum(unlist(terms))
    logPost.t <- resmc$draws[,"log.Poster"]
    
    jpeg(output.logPost, quality=100,width=960,height=960)
    plot(as.numeric(logPost.t),type='l',ylab="log(Posterior)",ylim=range(c(true.logPost,logPost.t)))
    abline(h=true.logPost,col="red")
    dev.off()
    
    ### ADD to true.par vector!
    true.par <- c(true.par, "log.Poster"=true.logPost)
    
    if(verbose){
      cat("log.Poster.values"); print(terms)
    } 
  }
  
  ########################################################
  ############## STORE RESULT ############################
  ## Write the final quantities to file:
  final_product <- list("post.q"=post.q,"true.par"=true.par,"post.q.S.mis"=post.q.S.mis,"true.S.mis"=true.par.S.mis)
  
  if (Sys.info()["sysname"]=="Darwin"){
    save.image(save.image.file)     
  } else if (Sys.info()["sysname"]=="Linux"){
    # Save only selected results
    if (!save_all_draws) {
      save(mcmc.stats, input.pars, eff.size, model, loglik, final_product, file=save.image.file)  
    } else {
      save(mcmc.stats, input.pars, eff.size, model, loglik, final_product, resmc, file=save.image.file)  
    }
  }
  
  cat("Making LNS plot...\n")
  
  # LogN-logS of posterior draws
  
  logS.truth <- log( c(res$par$S.obs,res$mis$S.mis) ,base=10)
  if(niter-burnin<1000){
    max.lines.jpg <- 0.1*(niter-burnin)
  } else{
    max.lines.jpg <- 2000
  }
  
  beginLNSPlotTime <- proc.time()
  
  jpeg(output.lnsplot, quality=100,width=960,height=960)
  plot.lns(resmc=resmc, logS.truth=logS.truth, max.lines=max.lines.jpg)
  dev.off()
  
  endLNSPlotTime <- proc.time()
  cat(sprintf("LogN-LogS plot time (%.1f seconds)\n", round((endLNSPlotTime-beginLNSPlotTime)[3], digits=1))) 
  
  if(do.mcmc.plots){
    cat("Making trace plots...\n")
    
    # Trace plots
    jpeg(output.mcmcplot, quality=100,width=960,height=960)
    par(mfrow=c(3,2))
    trace.mcmc.plot(resmc=resmc, res=res, true.par=true.par, max.mcmc.plots=max.mcmc.plots, remote.copy=FALSE)
    par(mfrow=c(1,1))
    dev.off()    
    
    
    #     cat("Making prior distr. plots...\n")
    #     pdf(output.priorsplot)
    # #     if(!do.prior.plots.complex) {
    #        prior.plots(resmc=resmc, true.par=true.par)
    # #     } else {
    #       posterior.prior.plots(resmc=resmc, true.par=true.par)
    #    # }
    #     dev.off()  
    #     
  }
  
  endTime <- proc.time()
  
  cat(paste("Finished MCMC sensitivity iteration ",job_num,"\n",sep=""))
  cat(sprintf("Elapsed time (%.1f seconds)\n", round((endTime-startTime)[3], digits=1))) 
  
  cat("\nAll done :)\n\n")
  
} ## end analyzing job (i.e. choice of sigma)


if (sinkit) {
  sink()
}

###############################
###############################


# 
# 
# ### Testing pi access from table:
# N.temp <- 100
# theta.temp <- seq(1.1, 1.8, length.out=N.temp)
# Smin.temp <- seq(4*10^-14, 4*10^-13, length.out=N.temp)
# bp.temp <- NULL
# 
# startTime <- proc.time()
# 
# "f1" <- function(theta.temp=theta.temp, Smin.temp=Smin.temp, bp.temp=bp.temp){
#   for(pp in 1:N.temp) {
#     pi.value <- pi.theta.get(pi=pi, theta=theta.temp[pp], Smin=Smin.temp[pp], bp=bp.temp, gamma=gamma, E=E,
#                              g=g.func, nsamples=nsamples, sigma=sigma, verbose=FALSE, method=1)     #marginal prob. of observing sources  
#   }
# }
# "f2" <- function(theta.temp=theta.temp, Smin.temp=Smin.temp, bp.temp=bp.temp){
#   for(pp in 1:N.temp) {
#     pi.value <- pi.theta.get(pi=pi, theta=theta.temp[pp], Smin=Smin.temp[pp], bp=bp.temp, gamma=gamma, E=E,
#                              g=g.func, nsamples=nsamples, sigma=sigma, verbose=FALSE, method=2)     #marginal prob. of observing sources  
#   }
# }
# 
# endTime <- proc.time()
# cat(sprintf("Retrieval of pi using recursive findInterval (%.4f seconds)\n", round((endTime-startTime)[3], digits=1))) 
# 
# 
# startTime <- proc.time()
# 
# for(pp in 1:N.temp) {
#   var.grid  <- pi$grid
#   new.coord <- matrix(c(Smin.temp[pp], bp.temp, theta.temp[pp]),nrow=1)
#   pi.value  <- pi$pi[nn2(data=var.grid,  query=new.coord,  k=1,  treetype="kd",  searchtype="priority", eps=10^-1)$nn.idx]
#      #marginal prob. of observing sources  
# }
# 
# endTime <- proc.time()
# cat(sprintf("Retrieval of pi using k-d trees (%.4f seconds)\n", round((endTime-startTime)[3], digits=1))) 
# 
# 
# print(pi.value)
# 
# 
# microbenchmark( f1(theta.temp, Smin.temp, bp.temp),
#                 f2(theta.temp, Smin.temp, bp.temp))
# 
# 
# 


