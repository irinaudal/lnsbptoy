## Small simulation example to demonstrate the bias induced due to pi approximation.

# Model: 
# Y_i ~ normal(Tau_i,1)
# Tau_i ~ gamma(Alpha,1)
# Alpha ~ gamma(20,8)
# I_i ~ Bin(1,g(Tau_i))
# g(Tau_i) = P(observe i-th source) = const
# pi(Alpha) = \int g(Tau_i) * p(Tau_i | Alpha)  dTau_i 
#    \approx  g+e,   where e ~ N(0, sigma)
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
  library(coda)
} else if (Sys.info()["sysname"]=="Linux"){
  libdir <- "/home/isudal/R/x86_64-pc-linux-gnu-library/3.0/"
  library(coda,lib.loc=libdir)
  main.dir   <- "/home/isudal/lns" #lns_data/"
  setwd(main.dir)
  libsrc.dir <- "."
}

## Handle batch job arguments:
args <- commandArgs(TRUE)  # 1-indexed version is used now.
#args <- 1
cat(paste("Command-line arguments:\n"))
print(args)

###################
sim_start <- 1000
if (length(args)==0){
  sinkit <- TRUE #FALSE
  sim_num <- sim_start + 1
} else {
  sinkit <- TRUE
  sim_num <- sim_start + as.numeric(args[1])
}
#set.seed(762*sim_num + 1330931) 
###################

# Set simulation parameters
length.datasets <- 200
plot.to.file    <- TRUE #FALSE #TRUE
niter  <- 55000
burnin <- 5000
tune.iter <- 100
stop.tune <- 3000

N.vec <- c(100,200,500)

prior.type.choice <- "gamma"                   # prior assumption: functional form
fixed.alpha.vec <- 2.5 #c(1.5, 2.5, 3.0, 3.5, 4.0)  # definition for gamma range
fixed.pars1 <- expand.grid("alpha"=fixed.alpha.vec,"N"=N.vec,"prior"=prior.type.choice)

prior.type.choice <- "normal"             # prior assumption: functional form
fixed.alpha.vec <- numeric(0) #c(0.1, 1.0, 2, 6, 10)  # definition for gamma range
fixed.pars2 <- expand.grid("alpha"=fixed.alpha.vec,"N"=N.vec,"prior"=prior.type.choice)

# merge settings
fixed.pars <- rbind(fixed.pars1,fixed.pars2)
N   <- fixed.pars$N[sim_num-sim_start]   # 100 # will be different
alpha <- fixed.pars$alpha[sim_num-sim_start] # 2.27 # will be different
prior.type <- fixed.pars$prior[sim_num-sim_start]

outer.dir.extra     <- paste("/TAU_Toy_example_3level_final_",prior.type, sep="")
output.dir.extra    <- paste("/data_N_",N,"_alpha_",alpha,sep="")
output.dir.outer    <- paste(main.dir,outer.dir.extra,sep="")
output.dir          <- paste(output.dir.outer,output.dir.extra,sep="")
output.dir.wo.main  <- paste(outer.dir.extra,output.dir.extra,sep="")
sink.file         <- file.path(output.dir, paste("MSE_result_N_",N,"_alpha_",alpha,".txt",sep=""))
output.plot.trace <- file.path(output.dir, paste("Plot_trace_N_",N,"_alpha_",alpha,"_%02d.jpg",sep=""))
output.plot.mse   <- file.path(output.dir.outer, paste("Plot_MSE_alpha_",alpha,"_%02d.jpg",sep=""))
if (!file.exists(output.dir.outer)){   #set directory - done after simulate()
  dir.create(output.dir.outer, recursive=TRUE)
  #   setwd(output.dir)
}
if (!file.exists(output.dir)){   #set directory - done after simulate()
  dir.create(output.dir, recursive=TRUE)
  #   setwd(output.dir)
}

# Use Normal additive error in approximation of the term pi(Tau)
sigma.vec <- c(0, 0.00001, 0.0001, 0.001, 0.004, 0.01, 0.05, 0.1, 0.2)
subset <- c(1:9)
sigma.vec <- sigma.vec[subset]
use.const <- FALSE
const     <- 0.02
# Incompleteness
"g" <- function(tau) {
  return(0.7)   # constant incompleteness
  #return(1 - exp(-tau))   # 1 minue exponential decay
}
# hyperparameters
sd.y <- 1
sd.tau <- 1
a    <- 20
b    <- 8
v.a  <- 0.5 #proposal sd for alpha
v.t  <- 0.7 #proposal sd for tau

if(prior.type=="normal"){
  "sampl.prior" <- function(n, a, b){
    return( rnorm(n, mean=a, sd=b) )
  }
  "log.dens.prior" <- function(x, a, b){
    return( dnorm(x, mean=a, sd=b, log=TRUE) )
  }
} else if(prior.type=="gamma"){
  "sampl.prior" <- function(n, a, b){
    return( rgamma(n, shape=a, rate=b))
  }
  "log.dens.prior" <- function(x, a, b){
    return( dgamma(x, shape=a, rate=b, log=TRUE) )
  }
}

"tune.v" <- function(v.a, v.t, prop.ct.tau, prop.ct.alpha, prop.state, tune.iter){
  # tune.v   Reduce or increase tuning SD parameters according for acceptance of proposals to be within 20-60%
  
  #   if (verbose) {
  #     cat("_____Tune v parameters_____\n") 
  #   } 
  n <- length(v.t)
  prop.state.t <- prop.ct.tau   - prop.state$prop.ct.tau
  prop.state.a <- prop.ct.alpha - prop.state$prop.ct.alpha
  tpa <- round(prop.state.t/tune.iter,3)*100      
  apa <- round(prop.state.a/tune.iter,3)*100      
  v.t[tpa<20] <- v.t[tpa<20]/2        #reject too many => reduce SD
  v.t[tpa>60] <- v.t[tpa>60]*2        #accept too many => increase SD     
  v.a[apa<20] <- v.a[apa<20]/2
  v.a[apa>60] <- v.a[apa>60]*2
  prop.state <- list("prop.ct.alpha"=prop.ct.alpha, "prop.ct.tau"=prop.ct.tau) #update current proposal state
  
  return(list("v.a"=v.a,"v.t"=v.t,"prop.state"=prop.state))
}

"met.tau" <- function(tau.curr, n, v.t, Yobs, sd.y, sd.tau, alpha.curr) {

  # Use Metropolis to sample a parameter tau_i
  tau.prop <- rnorm(n, mean=tau.curr, sd=v.t)
  if (any(tau.prop<=0)) { # bad proposal, cannot be negative
    p.tau.prop <- -Inf
  } else {
    p.tau.prop <- sum( dnorm(x=Yobs, mean=tau.prop, sd=sd.y, log=TRUE)  + 
                         dgamma(x=tau.prop, shape=alpha.curr, rate=sd.tau, log=TRUE) )
  }
  p.tau.curr <- sum( dnorm(x=Yobs, mean=tau.curr, sd=sd.y, log=TRUE)  + 
                       dgamma(x=tau.curr, shape=alpha.curr, rate=sd.tau, log=TRUE) )
  log.tau.alpha <- p.tau.prop - p.tau.curr
  if(is.nan(log.tau.alpha)) {
    log.tau.alpha <- -Inf
  }
  log.tau.u <- log(runif(1))
  
  # Decide to accept the pproposal
  got.draw.idx.tau <- (log.tau.u < log.tau.alpha)            #accepted proposal indicator 
  
  if(got.draw.idx.tau){
    tau.t <- tau.prop
  } else{
    tau.t <- tau.curr
  }
  
  return(list('tau.t'=tau.t,'idx'=got.draw.idx.tau))   
}

"met.alpha" <- function(alpha.curr, v.a, Yobs, n, N, tau.t, sd.tau, a, b, prior.type, g, 
                        sigma, const, use.const) {
  
  # Use Metropolis to sample a parameter alpha
  alpha.prop <- rnorm(1, mean=alpha.curr, sd=v.a)
  
  # evaluate pi with error
  if(use.const==FALSE) {
    pi.prop <- min(1, g(min(tau.t)) + rnorm(1,mean=0,sd=sigma))
    pi.curr <- min(1, g(min(tau.t)) + rnorm(1,mean=0,sd=sigma))
  } else {
    pi.prop <- min(1, g(min(tau.t)) + const)
    pi.curr <- min(1, g(min(tau.t)) + const)
  }
  
  #     cat("tau.t[k,]\n")
  #     print(tau.t[k,])
  #     cat("alpha.prop\n")
  #     print(alpha.prop)
  #     cat("pi.prop\n")
  #     print(pi.prop)
  #     cat("\n")
  
  # log-full-conditionals used in MH proposal ratio
  if (alpha.prop<=0) {
    p.prop <- -Inf
  } else {
    p.prop <- sum( dgamma(tau.t, shape=alpha.prop, rate=sd.tau, log=TRUE) ) + 
      log.dens.prior(alpha.prop, a=a, b=b) +
      (N-n)*log(1-pi.prop)
  }
  p.curr <- sum( dgamma(tau.t, shape=alpha.curr, rate=sd.tau, log=TRUE) ) + 
    log.dens.prior(alpha.curr, a=a, b=b) +
    (N-n)*log(1-pi.curr)
  
  log.alpha <- p.prop - p.curr
  if(is.nan(log.alpha)) {
    log.alpha <- -Inf
  }
  log.u <- log(runif(1))
  
  # Decide to accept the pproposal
  got.draw.idx <- (log.u < log.alpha)            #accepted proposal indicator 
  
  if(got.draw.idx){
    alpha.t <- alpha.prop
  } else{
    alpha.t <- alpha.curr
  }
  
  return(list('alpha.t'=alpha.t,'idx'=got.draw.idx))   
}


"mc" <- function(alpha.curr, tau.curr, niter, sigma, v.a, v.t, prop.ct.alpha, prop.ct.tau, prop.ct.state, 
                 tune.iter, stop.tune, sd.tau, prior.type, sd.y, g, const, use.const, Yobs, N ) {
  
  n <- length(Yobs)
  alpha.t    <- numeric(niter)
  alpha.t[1] <- alpha.curr
  tau.t      <- matrix(NA,nrow=niter,ncol=n)
  tau.t[1,]  <- tau.curr
  
  for (k in 2:niter) {
    alpha.curr <- alpha.t[k-1]
    tau.curr   <- tau.t[k-1,]
    
    # Use Metropolis to sample a parameter tau_i
    met <- met.tau(tau.curr=tau.curr, n=n, v.t=v.t, Yobs=Yobs, sd.y=sd.y, sd.tau=sd.tau, alpha.curr=alpha.curr)
    tau.t[k,] <- met$tau.t
    prop.ct.tau <- prop.ct.tau + met$idx
    
    # Use Metropolis to sample a parameter alpha
    met <- met.alpha(alpha.curr=alpha.curr, v.a=v.a, Yobs=Yobs, n=n, N=N, tau.t=met$tau.t, sd.tau=sd.tau, 
                     a=a, b=b, prior.type=prior.type, g=g, 
                     sigma=sigma, const=const, use.const=use.const)
    alpha.t[k] <- met$alpha.t
    prop.ct.alpha <- prop.ct.alpha + met$idx
    
    ############################
    #Tune v.S, v.theta parameters every tune.iter iterations
    if (k%%tune.iter == 0 & k <= stop.tune){
      tune <- tune.v(v.t=v.t, v.a=v.a, prop.ct.tau=prop.ct.tau, prop.ct.alpha=prop.ct.alpha,
                     prop.state=prop.ct.state, 
                     tune.iter=tune.iter)
      prop.ct.state <- tune$prop.state
      v.t <- tune$v.t
      v.a <- tune$v.a
      
#       cat("k:\n")
#       print(k)
#       cat("v.a, v.t:\n")
#       print(c(v.a, v.t))
#       cat("prop.state:\n")
#       print(prop.ct.state)
#       cat("\n")
    }    
    
    
  } # END for loop of k=iter
  
  return(list("alpha.t"=alpha.t,"tau.t"=tau.t))
} ## end mcmc step


########################################################################
####################### RUN ############################################
if (sinkit) {
  cat("Sinking to...\n")
  print(sink.file)
  sink(sink.file)
}

length.jobs <- length(sigma.vec)
temp.matrix <- matrix(NA, nrow=length.datasets, ncol=length.jobs)
res_b <- list("mean"=temp.matrix,"sd"=temp.matrix,"var"=temp.matrix,
              "median"=temp.matrix,
              "min"=temp.matrix,"max"=temp.matrix,
              "q025"=temp.matrix,"q975"=temp.matrix)
res <- rep(NA,length.jobs)

if (plot.to.file) {
  jpeg(output.plot.trace, quality=100,width=960,height=960)
}

for (data_num in 1:length.datasets) { # dataset is unique for each job.
  
  # Generate data: complete, observed and missing
  set.seed(762*data_num + 1330931) # set sid for producing the datasets
  
  tau <- rgamma(N, shape=alpha,rate=1)
  Y <- rnorm(N, mean=tau,sd=1)
  I <- rbinom(N, size=1, prob=g(min(tau)))
  Yobs <- Y[I==1]
  n <- length(Yobs)
  
  for (job_num in 1:length.jobs ) {
    
    sigma <- sigma.vec[job_num]
    res[job_num] <- c(alpha)
  
    prop.ct.alpha <- numeric(1)
    prop.ct.tau   <- numeric(1)
    prop.ct.state <- list("prop.ct.alpha"=prop.ct.alpha,"prop.ct.tau"=prop.ct.tau)
    
    tau.init   <- rep(mean(Yobs),n) # estimate starting state
    alpha.init <- mean(tau.init)
    
    mc_b <- mc(alpha.curr=alpha.init, tau.curr=tau.init, niter=niter, sigma=sigma, v.a=v.a, v.t=v.t,
               prop.ct.alpha=prop.ct.alpha,prop.ct.tau=prop.ct.tau, prop.ct.state=prop.ct.state, 
               tune.iter=tune.iter, stop.tune=stop.tune,
               sd.tau=sd.tau, sd.y=sd.y, prior.type=prior.type, g=g, const=const, use.const=use.const, Yobs=Yobs, N=N)
    alpha_b <- mc_b$alpha.t
    
    alpha_b <- alpha_b[(burnin+1):niter]
    
    # Store statistics
    res_b$mean[data_num,job_num] <- mean(alpha_b)
    res_b$sd[data_num,job_num]   <- sd(alpha_b)
    res_b$var[data_num,job_num]  <- var(alpha_b)
    res_b$median[data_num,job_num] <- median(alpha_b)
    res_b$min[data_num,job_num] <- min(alpha_b)
    res_b$max[data_num,job_num] <- max(alpha_b)
    res_b$q025[data_num,job_num] <- quantile(alpha_b,prob=0.025)
    res_b$q975[data_num,job_num] <- quantile(alpha_b,prob=0.975)
        
    # Save to file
    output.result.image <- file.path(output.dir,paste("Store_res_job_num_",job_num,".RData",sep=""))
    save(res_b, res, niter,burnin, a,b,sigma,prior.type,g,N,sd.y,  file=output.result.image)  
    
    if(data_num == length.datasets) {
      # Show results of the last fitted dataset
      par(mfrow=c(2,2))
      
      alpha_b <- mcmc(alpha_b)
      range.alpha <- range(c(alpha_b,alpha))
      traceplot(alpha_b,main=paste("Posterior: alpha_b. N=",N,", n=",n,", N-n=",N-n,", sigma=",sigma,sep=""))
      abline(h=alpha,col="red",lwd=1.8)
      densplot(alpha_b,main=paste("Posterior density: alpha_b"),xlim=range.alpha)
      abline(v=alpha,col="red",lwd=1.8)

      par(mfrow=c(1,1))
    }
    
  } ## end analyzing job (i.e. choice of sigma)
} ## end analyzing dataset

if (plot.to.file) {
  dev.off()
}


# Show results of the MSE of all datasets
# res_b is [length.datasets x n.stats] matrix
posterior_est_lower <- apply(res_b$q025,2,mean,na.rm=TRUE)
posterior_est_upper <- apply(res_b$q975,2,mean,na.rm=TRUE)
posterior_est_bar_var <- apply(res_b$var,2,mean,na.rm=TRUE)

"evaluate_results" <- function(posterior_est, est_title, model_title, res, N, sigma.vec, length.jobs){
  #Evaluate MSE statistics (average across datasets)
  # effective sample size is N, no NA cases.
  posterior_est_bar  <- apply(posterior_est,2,mean,na.rm=TRUE) # mean of means across datasets
  posterior_est_var  <- apply(posterior_est,2,var,na.rm=TRUE) # variance of means across datasets
  posterior_est_bias <- posterior_est_bar - res  # bias of mean of means
  posterior_est_mse  <- posterior_est_var*(N-1)/N + posterior_est_bias^2 # MSE of mean of means
  absbias <- sapply(1:length.jobs, function(x)abs(posterior_est[,x]-res[x])) # absolute bias of each dataset
  posterior_est_absbias <- apply(absbias,2,mean,na.rm=TRUE) # mean of absolute bias across datasets
  
  mse <- data.frame("mse"=round(posterior_est_mse,4),
                    "bias"=round(posterior_est_bias,4),
                    "absbias"=round(posterior_est_absbias,4),
                    "var"=round(posterior_est_var,4),
                    "mean"=round(posterior_est_bar,4),
                    "model"=factor(rep(model_title,length.jobs)),
                    "job"=factor(sigma.vec),
                    "est_title"=factor(rep(est_title,length.jobs)),
                    "N"=N)
  return(mse)
}

# compute summary statistics for posterior mean estimate
mse_b_mean <- evaluate_results(posterior_est=res_b$mean, est_title="Mean", model_title="Approx Pi", res, N=N, sigma.vec=sigma.vec, length.jobs=length.jobs)

# compute summary statistics for posterior median estimate
mse_b_median <- evaluate_results(posterior_est=res_b$median, est_title="Median", model_title="Approx Pi", res, N=N, sigma.vec=sigma.vec, length.jobs=length.jobs)

mse <- rbind(mse_b_mean)
mse2 <- rbind(mse_b_median)

cat("\n\nTrue value of alpha:\n")
print(alpha)
cat("N:\n")
print(N)
cat("sigma.vec:\n")
print(sigma.vec)
cat("Results for Mean:\n")
print(mse)
cat("Results for Median:\n")
print(mse2)

# Save to file
output.mse.image <- file.path(output.dir,paste("Store_mse.RData",sep=""))
save(mse, mse2, res, niter,burnin, a,b,sigma.vec,prior.type,g,N,sd.y,  file=output.mse.image)  

sink()

###############################
###############################

# after all results are done, make pretty pictures of MSE
if(N==N.vec[length(N.vec)]) {
  
  # read all other files and combine results
  main.mse <- main.mse2 <- numeric(0)
  
  for(N in N.vec) {
    output.dir.extra <- paste("/data_N_",N,"_alpha_",alpha,sep="")
    output.dir       <- paste(output.dir.outer,output.dir.extra,sep="")
    output.mse.image <- file.path(output.dir,paste("Store_mse.RData",sep=""))
    load(output.mse.image)
    main.mse  <- rbind(main.mse,  mse)
    main.mse2 <- rbind(main.mse2, mse2)
  }
  
  "mse_plot" <- function(mse, res, length.datasets, N, sigma) {
    par(mfrow=c(4,1))
    interaction.plot(x.factor=mse$job, trace.factor=mse$N, mse$mse, ylim=range(c(0,mse$mse)), 
                     type="b", pch=c(16,17,15), col=1:3, lwd=1.5,
                     trace.label="Population size: N",
                     ylab="MSE", xlab="Sigma",
                     main=paste("Posterior ",mse$est_title[1]," estimate of Alpha=",alpha,", averaged across ",length.datasets," datasets,",
                                " with error on incompleteness of size Sigma",sep=""),
                     fixed=TRUE, leg.bty="o", leg.bg="beige")
    abline(h=0, lty=1)
    
    interaction.plot(x.factor=mse$job, trace.factor=mse$N, mse$var, ylim=range(c(0,mse$var)), 
                     type="b", pch=c(16,17,15), col=1:3, lwd=1.5,
                     trace.label="Population size: N",
                     ylab="Var", xlab="Sigma",
                     fixed=TRUE, leg.bty="o", leg.bg="beige")
    abline(h=0, lty=1)
    interaction.plot(x.factor=mse$job, trace.factor=mse$N, mse$absbias, ylim=range(c(0,mse$absbias)), 
                     type="b", pch=c(16,17,15), col=1:3, lwd=1.5,
                     trace.label="Population size: N",
                     ylab="|Bias|", xlab="Sigma",
                     fixed=TRUE, leg.bty="o", leg.bg="beige")
    abline(h=0, lty=1)
    interaction.plot(x.factor=mse$job, trace.factor=mse$N, mse$mean, ylim=range(c(res[1],mse$mean)), 
                     type="b", pch=c(16,17,15), col=1:3, lwd=1.5,
                     trace.label="Population size: N",
                     ylab=mse$est_title[1], xlab="Sigma",
                     fixed=TRUE, leg.bty="o", leg.bg="beige")
    abline(h=res[1], lty=1, lwd=1, col="blue")
    text(x=10.8,y=res[1]*1.02,labels="True Alpha",col="blue",cex=1.3)
    par(mfrow=c(1,1))
  }
  
  if (plot.to.file) {
    jpeg(output.plot.mse, quality=100,width=960,height=960)
  }
  
  mse_plot(mse=main.mse, res=res, length.datasets=length.datasets, N=N, sigma=sigma)
  mse_plot(mse=main.mse2, res=res, length.datasets=length.datasets, N=N, sigma=sigma)
  
  if (plot.to.file) {
    dev.off()
  }
  
}