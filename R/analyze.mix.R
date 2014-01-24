"analyze.mix" <- function(Y, niter, burnin, v.so, v.th, v.sm=NULL, v.bp=NULL, 
                          alpha, beta, a, b, gamma, E, 
                          C=NULL, mu=NULL, am=NULL, bm=NULL,
                          nsamples=10000, sigma=0,
                          g=function(lambda){  return(rep(0.7,length(lambda)))  },
                          fixed.N=FALSE, fixed.theta=FALSE, fixed.S.obs=FALSE, fixed.S.mis=FALSE, 
                          fixed.Smin=FALSE, fixed.bp=NULL, 
                          tune.iter=100, stop.tune=1000, store_logPost=FALSE,
                          print.every=10000, save.every=1, verbose=FALSE, save.progress=FALSE, save.progress.dir="R_code" ){
  
  inputs <- list("g"=g,"nsamples"=nsamples, "sigma"=sigma, "E"=E,
                 "Y.obs.tot"=Y,"niter"=niter,"burnin"=burnin,"v.so"=v.so,"v.th"=v.th,"v.sm"=v.sm,"v.bp"=v.bp,
                 "a"=a,"b"=b,"C"=C,"mu"=mu,"alpha"=alpha,"beta"=beta,"am"=am,"bm"=bm,"gamma"=gamma,
                 "fixed.theta"=fixed.theta,"fixed.S.obs"=fixed.S.obs,"fixed.N"=fixed.N,
                 "fixed.S.mis"=fixed.S.mis,"fixed.Smin"=fixed.Smin,"fixed.bp"=fixed.bp,
                 "save.every"=save.every, "store_logPost"=store_logPost,
                 "tune.iter"=tune.iter, "stop.tune"=stop.tune, "verbose"=verbose)
  
  
  ####################################################################################################
  #               Methods: Gibbs sampler and Metropolis-Hastings algorithm
  #
  # Input: Y         = total observed photon counts
  #        niter     = total number of iterations per chain
  #        burnin    = burn-in number of iterations
  #        v.so       = value SD, tuning parameter for fluxes S to accept 20-60% of proposals
  #        v.th       = value SD, tuning parameter for slope theta to accept 20-60% of proposals
  #        v.sm       = value SD, tuning parameter for Smin to accept 20-60% of proposals
  #        v.bp       = value SD, tuning parameter for bp to accept 20-60% of proposals
  #        alpha     = target number of sccessful trials in negbinom prior for N
  #        beta      = dispersion parameter in negbinom prior for N
  #        a         = vector of shape hyper-parameter(s) in gamma prior for theta(s)
  #        b         = vector of rate hyper-parameter(s) in gamma prior for theta(s)
  #        am        = shape hyper-parameter(s) in gamma prior for Smin
  #        bm        = rate hyper-parameter(s) in gamma prior for Smin
  #        mu        = mean parameter in normal proposal for break-points
  #        C         = varianca-covariance matrix of parameters in normal proposal for break-points
  #        gamma     = constant of transformation, energy per photon
  #        pble      = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        E         = vector of exposureMap of observed sources
  #        fixed.N     = Value/FALSE: debugging for N - no mcmc
  #        fixed.theta = Value/FALSE: allow for user specified theta - no mcmc
  #        fixed.S.obs = Value/FALSE: debugging for S.obs - no mcmc
  #        fixed.S.mis = Value/FALSE: debugging for S.mis - no mcmc
  #        fixed.Smin  = Value/FALSE: allow for user specified Smin  - no mcmc
  #        fixed.bp    = either NULL (no fixed break-points), or, a vector of length p specifying the break-points
  #        nsamples    = number of MC samples to produce estimate of pi(theta,Smin) integral
  #        g           = function, probability of observing a source
  #        store_logPost = (T/F) whether to store estimate of log-posterior (without norm const)
  #        tune.iter   = How often to tune the proposal variance
  #        stop.tune   = when to stop tuning parameters
  #        save.every  = thin MCMC draws by this amount
  #        print.every = display code progress every kth iteration
  #        verbose     = (T/F) display progress of program
  #        save.progress      = manually save current draws to file
  #        save.progress.dir  = storage file for current draws
  #
  # Output:  list of :
  #             draws  = posterior samples of parameters (N, theta, Smin, bp, S.obs),
  #             prop.accept = accepted proportions for MH sampled parameters (S.obs, S.mis, theta, Smin, bp),
  #             draws.S.mis, inputs, eff.size, 
  #             loglik
  ####################################################################################################
  
  startTime <- proc.time()
  
  if (burnin>niter)
    stop("'burnin' must be less than 'niter'")
  if (niter<1)
    stop("'niter' must be positive")
  
  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  
  Y.obs.tot <- Y
  n <- length(Y.obs.tot)                         #number of observed sources
  m <- length(a)
  if (length(b) != m){
    stop("'a' and 'b' must be the same length")
  }
  
  if (verbose) {
    cat("Generating starting values.. \n") 
  } 
  
  ### Generate Smin, the minimum source flux, informative prior:
  if (any(!fixed.Smin)){     #some/all Smin values are not fixed and need to be generated
    # estimate
    S.obs <- ((Y.obs.tot)*gamma)/E     #vector of initial values: S.obs0
    subby <- (Y.obs.tot)>2 
    tmp.Sobs <- S.obs[subby]
    if (length(tmp.Sobs)==0){
      warning("Unable to produce valid starting state 'estimate' for Smin due to no fluxes (no temporary S.obs).")
      tmp.Smin <- max(c(min(0.5*gamma/E),min(S.obs)))
    } else {  
      tmp.Smin <- min(tmp.Sobs)
    }
    Smin.t <- tmp.Smin*0.95
  } else {
    Smin.t  <- fixed.Smin
  }
  if (all(!!fixed.bp)) {
    if (all(Smin.t>fixed.bp)) {
      warning("Initial Smin was estimated to be smaller than fixed.bp")
      Smin.t <- min(fixed.bp)*0.90
    }
  }
  if (verbose){
    cat("Smin.t=\n"); print(Smin.t)
  }
  
  
  ### Generate the breakpoints:
  # fixed.bp has 3 options:
  # Option 1 :: fixed.bp == NULL                 ==> No break-points
  # Option 2 :: fixed.bp == (vector of scalars)  ==> Break-points, fixed to user-specified values
  # Option 3 :: fixed.bp == FALSE                ==> Break-points, must be estimated
  
  if (any(!fixed.bp)){  
    # Option 3 :: Break-points, must be estimated
    bp.t <- cumsum(exp(rnorm(n=m-1, mean=mu, sd=C)))+Smin.t
    idx <- (fixed.bp != FALSE)          # index of fixed values
    bp.t[idx] <- fixed.bp[idx]       # overwrite with fixed values
  } else {
    # Option 2 :: Break-points, fixed to user-specified values
    bp.t <- fixed.bp
  }
  if (verbose){
    cat("bp.t=\n"); print(bp.t)
  }
  
  ### Generate total number of sources N:
  if (!fixed.N){
    prob.obs.source <- g(lambda=Smin.t*E/gamma)
    
    if(all(prob.obs.source>0.9)){
      #Probability of observing the source is very high, so start without missing sources
      N.t = n
    } else{
      N.t <- max(n,rnbinom(1,size=alpha,prob=beta/(1+beta)))     #initial value: N0
    }
    
  } else {
    if (fixed.N < n)
      stop("'fixed.n' is less than 'n'")
    N.t <- fixed.N
  }
  
  if (verbose) {
    cat("N.t = \n"); print(N.t)
    cat("n = \n"); print(n)
  } 
  
  ### Generate theta(s):
  if (any(!fixed.theta)){  
    #some/all theta values are not fixed and need to be generated
    theta.t <- rgamma(n=m,shape=a,rate=b)  #initial value: theta0
    idx <- (fixed.theta != FALSE)          # index of fixed values
    theta.t[idx] <- fixed.theta[idx]       # overwrite with fixed values
  } else {
    theta.t <- fixed.theta
  }
  
  ### Generate source flux S.obs:
  # Generalizing to handle both break-points and untruncated mixtures
  if (verbose){
    cat('Debugging S.obs generation\n')
    cat('!fixed.S.obs:'); print(!fixed.S.obs)
    cat('!fixed.S.mis:'); print(!fixed.S.mis)
    cat('Smin.t:'); print(Smin.t)
    cat('theta.t:'); print(theta.t)
  }
  if (any(!fixed.S.obs)){
    tmp <- (Y.obs.tot*gamma)/E     #vector of initial values: S.obs0
    tmp[tmp<Smin.t] <- Smin.t*1.0001
    S.obs.t <- tmp
    
    if (verbose){
      cat("Initial values of S.obs.t, before patching: \n")
      print(S.obs.t)
    }
    
    # Check starting states of the observed sources are all valid...
    S.obs.t <- patch.sobs(S.obs.t=S.obs.t, theta.t=theta.t, Smin=Smin.t, bp=bp.t,
                          E=E, gamma=gamma, g=g, verbose=verbose)
    if (verbose) {
      cat("Done with patching Sobs. \n") 
      print(S.obs.t)
    }
    
  } else {
    S.obs.t <- fixed.S.obs
  }
  
  ### Generate source flux S.mis:    
  if (any(!fixed.S.mis)){
    if (N.t==n){     
      S.mis.t <- numeric(0)  
      I.idx   <- numeric(0)
      
    } else if (N.t>n){
      
      # Generate from broken power-law mixture of Truncated Pareto's:
      S.mis.t <- rbrokenpareto(n=N.t-n,x_min=Smin.t,k=theta.t,bp=bp.t,verbose=verbose)     #vector of initial values: S.mis0
      # Check starting states of the missing sources are all valid...
      S.mis.t <- patch.smis(S.mis.t=S.mis.t, theta.t=theta.t, Smin=Smin.t, bp=bp.t,
                            E=E, gamma=gamma, g=g, verbose=verbose)
      if (verbose) {
        cat("Done with patching Smis. \n") 
      }
    }
  } else {
    S.mis.t <- fixed.S.mis
  }
  
  if (length(v.so)<1){
    stop("v.so (SD parameter for MH proposal for flux S) has not been specified as a scalar or is unknown. Please original check settings")
  } else if (length(v.so)==1){
    v.so <- rep(v.so,n)              #v.so = tuning parameter (vector) to accept 20-60% of proposals
  } else {
    if (length(v.so)!=n) {
      v.so <- rep(v.so[1],n)
    }
  } 
  if (length(v.th)<1){
    stop("v.th :theta (SD parameter for MH proposal for theta) has not been specified as a scalar or is unknown. Please original check settings")
  } else if (length(v.th)==1) {
    v.th <- rep(v.th,m)      #v.thheta = tuning parameter (vector) to accept 20-60% of proposals
  } else {
    if (length(v.th)!=m) {
      v.th <- rep(v.th[1],m)
    }
  } 
  if (length(v.bp)<1){
    stop("v.bp :bp.t (SD parameter for MH proposal for preak-point) has not been specified as a scalar or is unknown. Please original check settings")
  } else if (length(v.bp)==1){
    v.bp <- rep(v.bp,m-1)              #v.bp.t = tuning parameter (vector) to accept 20-60% of proposals
  } else {
    if (length(v.bp)!=m-1) {
      v.bp <- rep(v.bp[1],m-1)
    }
  } 
  
  prop.ct.S.obs <- numeric(n)    #count of accepted proposals for S.obs
  prop.ct.theta <- numeric(m)    #count of accepted proposals for theta
  prop.ct.Smin  <- numeric(1)    #count of accepted proposals for Smin
  prop.ct.bp   <- numeric(m-1)  #count of accepted proposals for bp
  prop.ct.S.mis <- numeric(N.t*2)        #count of accepted proposals for S.mis, larger than needed
  prop.state <- list("prop.ct.theta"=prop.ct.theta,
                     "prop.ct.Smin"=prop.ct.Smin,"prop.ct.bp"=prop.ct.bp,
                     "prop.ct.S.obs"=prop.ct.S.obs)   #current state of proposals
  
  ### Set-up storage matrix for 'draws'
  total.samples <- floor((niter-burnin)/save.every)
  len.draws <- 1+m # for N and theta
  name.draws <- c("N",paste("theta.",1:m,sep=""))
  # Add other parameters: tau.1=Smin, tau.j=bp, p.mix, S.obs.j
  if (!fixed.Smin) {
    len.draws <- len.draws+1 # for Smin
    name.draws <- c(name.draws, "tau.1")
  }
  if (any(!fixed.bp)) {
    len.draws <- len.draws+m-1 # for bp.t
    name.draws <- c(name.draws, paste("tau.",2:m,sep=""))
  }
  len.draws <- len.draws+n # for S.obs
  name.draws <- c(name.draws, paste("S.obs.",1:n,sep=""))
  
  draws <- mcmc(matrix(numeric(0),nrow=total.samples,ncol=len.draws))
  varnames(draws) <- name.draws
  
  # Add log-posterior results
  if (store_logPost) {
    draws <- mcmc(cbind(draws,"log.Poster"=NA))
  }
  draws.S.mis <- vector("list",total.samples)  
  
  loglik <- list("loglik"=numeric(total.samples))
  
  if (verbose) {
    cat("Generating draws.. \n") 
  } 
  
  save.i <- 1
  for (iter in 1:niter) {
    
    update <- update.lns(S.obs.t=S.obs.t, theta.t=theta.t, Smin.t=Smin.t, bp.t=bp.t, N.t=N.t, n=n, Y.obs.tot=Y.obs.tot, 
                         S.mis.t=S.mis.t, v.so=v.so, v.th=v.th, v.sm=v.sm, v.bp=v.bp, niter=niter,
                         prop.ct.S.obs=prop.ct.S.obs, prop.ct.S.mis=prop.ct.S.mis,
                         prop.ct.theta=prop.ct.theta, prop.ct.Smin=prop.ct.Smin, prop.ct.bp=prop.ct.bp,
                         a=a, b=b, C=C, mu=mu, alpha=alpha, beta=beta, am=am, bm=bm, gamma=gamma, 
                         E=E, nsamples=nsamples, g=g, sigma=sigma,
                         fixed.N=fixed.N, fixed.theta=fixed.theta, fixed.Smin=fixed.Smin, fixed.bp=fixed.bp,
                         fixed.S.obs=fixed.S.obs, fixed.S.mis=fixed.S.mis, 
                         store_logPost=store_logPost,
                         verbose=verbose)
    
    # Update all variables
    N.t     <- update$draws$N.t
    theta.t <- update$draws$theta.t
    S.obs.t <- update$draws$S.obs.t
    
    if (any(!fixed.Smin)){
      Smin.t <- update$draws$Smin.t
    }
    if(any(!fixed.bp)){
      bp.t   <- update$draws$bp.t 
    }  
    
    prop.ct.theta <- update$prop.ct.theta
    prop.ct.Smin  <- update$prop.ct.Smin
    prop.ct.bp   <- update$prop.ct.bp
    prop.ct.S.obs <- update$prop.ct.S.obs
    
    S.mis.t <- update$draws.S.mis
    prop.ct.S.mis <- update$prop.ct.S.mis
    
    #Tune v.so, v.thheta parameters every tune.iter iterations
    if (iter%%tune.iter == 0 & iter <= stop.tune){
      tune <- tune.v(v.th=v.th, v.sm=v.sm, v.bp=v.bp, v.so=v.so, prop.ct.theta=prop.ct.theta, 
                     prop.ct.Smin=prop.ct.Smin, prop.ct.bp=prop.ct.bp, 
                     prop.ct.S.obs=prop.ct.S.obs, prop.state=prop.state, 
                     tune.iter=tune.iter, verbose=verbose)
      prop.state <- tune$prop.state
      v.th <- tune$v.th
      v.sm <- tune$v.sm
      v.bp <- tune$v.bp
      v.so <- tune$v.so
    }    
    
    # Store draws after burnin (Automatically remove burnin) (Automatically thin draws)
    if (iter>burnin) {  
      if ((iter-burnin) %% save.every == 0) {
        draws[(iter-burnin)/save.every,] <- unlist(update$draws)
        draws.S.mis[[(iter-burnin)/save.every]] <- update$draws.S.mis
        
        # Evaluate the log( pr(Data|param(i)) ) in order to estimate of log of normalizing constant
        pi.value <- pi.theta.get(theta=theta.t, Smin=Smin.t, bp=bp.t, gamma=gamma, E=E, g=g,
                                 nsamples=nsamples, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
        lambda.obs <- S.obs.t*E/gamma
        pr.n <- dbinom(x=n, size=N.t, prob=pi.value, log=TRUE)
        pr.y <- dpois(x=Y.obs.tot, lambda=lambda.obs, log=TRUE)
        loglik$loglik[(iter-burnin)/save.every] <- pr.n + sum(pr.y) 
        
      }
    }
    
    if (iter%%print.every == 0){ 
      
      cat("N, theta, Smin, bp, and S[1] are:   ")
      print(c(update$draws$N.t,update$draws$theta.t,update$draws$Smin.t,update$draws$bp.t,update$draws$S.obs.t[1]))
      
      cat("______Show-Memory-Use_____\n")
      
      #############
      # Below is code for memory-use function, which MUST be used as a script in order to retrieve the objects 
      #   in the current function environment
      objectList <- ls(pos=-1)
      
      sort <- "size" 
      decreasing <- TRUE
      limit <- 10
      
      oneKB <- 1024
      oneMB <- 1048576
      oneGB <- 1073741824
      
      memoryUse <- sapply(objectList, function(x) as.numeric(object.size(eval(parse(text=x)))))
      
      memListing <- sapply(memoryUse, function(size) {
        if (size >= oneGB) return(paste(round(size/oneGB,2), "GB"))
        else if (size >= oneMB) return(paste(round(size/oneMB,2), "MB"))
        else if (size >= oneKB) return(paste(round(size/oneKB,2), "kB"))
        else return(paste(size, "bytes"))
      })
      
      memListing <- data.frame(objectName=names(memListing),memorySize=memListing,row.names=NULL)
      
      if (sort=="alphabetical") memListing <- memListing[order(memListing$objectName,decreasing=decreasing),] 
      else memListing <- memListing[order(memoryUse,decreasing=decreasing),] #will run if sort not specified or "size"
      
      if(!missing(limit)) memListing <- memListing[1:limit,]
      
      print(memListing, row.names=FALSE)
      
      rm(memoryUse,memListing)
      #############
      
      gc() # perform manual garbage collection
      
      cat(sprintf("MCMC computation time (%.1f seconds)\n", round((proc.time()-startTime)[3], digits=1)))    
    }
    
    if (save.progress){
      if (iter %% 50000 == 0){
        save( list = ls(all.names=TRUE,pos=-1), file = paste(save.progress.dir,"/dataset_",save.i,"_after_",iter,"_iterations.RData",sep=""))
        save.i <- save.i+1
      }
    }
    
    if (iter%%print.every == 0)
      cat(paste("Finished iteration ",iter,"\n",sep="")) 
    
    if (verbose){
      cat("________________________________\n")
    }
    
  } # END iter for-loop across iterations of MCMC
  
  ############################################
  # Compute effective sample size
  eff.size <- effectiveSize(draws)
  ############################################
  
  tpa <- round(prop.ct.theta/niter,2)*100
  mpa <- round(prop.ct.Smin/niter,2)*100
  epa <- round(prop.ct.bp/niter,2)*100
  spa <- round(prop.ct.S.obs/niter,2)*100
  prop.accept <- list("v"=c(v.so,v.th,v.sm,v.bp),
                      "S.prop.accept"=spa,"theta.prop.accept"=tpa,"Smin.prop.accept"=mpa,"eta.prop.accept"=epa,
                      "S.mis.ct.accept"=prop.ct.S.mis)
  
  return(list("draws"=draws,"prop.accept"=prop.accept,
              "draws.S.mis"=draws.S.mis,
              "inputs"=inputs,"eff.size"=eff.size,
              "loglik"=loglik))
}

analyze.mix <- cmpfun(analyze.mix)
  