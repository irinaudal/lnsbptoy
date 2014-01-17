"analyze.mix" <- function(Y, niter, burnin, v.S, v.theta, v.sm=NULL, v.bp=NULL, 
                          alpha, beta, a, b, gamma, pble, E.obs,L.obs,bg.obs,k.obs, 
                          C=NULL, mu=NULL, d=NULL, am=NULL, bm=NULL, prob.MH=0, start.state="random",
                          pi=NULL, fixed.S.pi=NULL, length.S=100, debug.pi=FALSE, nsamples=10000, g.type="step", 
                          g=function(lambda,bg,E,L,g.type){
                            g.compute(lambda=lambda,bg=bg,E=E,L=L,g.type=g.type)
                          }, area.conv.func=function(off){10.57508 * exp(0.5232*off)}, delta.vec=seq(0,0.1,by=0.01),
                          fixed.N=FALSE, fixed.theta=FALSE, fixed.S.obs=FALSE, fixed.S.mis=FALSE, fixed.Smin=FALSE, fixed.bp=NULL, fixed.p=NULL, 
                          model="regular", tune.iter=100, stop.tune=1000, compute.dic.extras=FALSE,
                          met.alg=list("Smin"="MH","bp"="MH"), store_logPost=FALSE,
                          print.every=10000, save.every=1, verbose=FALSE, save.progress=FALSE, save.progress.dir="R_code" ){
  
  inputs <- list("pi"=pi,"pble"=pble,"area.conv.func"=area.conv.func,"g"=g,"g.type"=g.type,"length.S"=length.S,"nsamples"=nsamples, 
                 "E.obs"=E.obs,"L.obs"=L.obs,"bg.obs"=bg.obs,"k.obs"=k.obs,
                 "Y.obs.tot"=Y,"niter"=niter,"burnin"=burnin,"v.S"=v.S,"v.theta"=v.theta, "v.sm"=v.sm, "v.bp"=v.bp,
                 "a"=a,"b"=b,"C"=C,"mu"=mu,"d"=d,"alpha"=alpha,"beta"=beta,"am"=am,"bm"=bm,"gamma"=gamma,
                 "fixed.theta"=fixed.theta,"fixed.S.obs"=fixed.S.obs,"fixed.N"=fixed.N,
                 "fixed.S.mis"=fixed.S.mis,"fixed.Smin"=fixed.Smin,"fixed.bp"=fixed.bp,
                 "fixed.p"=fixed.p,"delta.vec"=delta.vec,"model"=model,"met.alg"=met.alg,
                 "start.state"=start.state, "save.every"=save.every, "store_logPost"=store_logPost,
                 "tune.iter"=tune.iter, "stop.tune"=stop.tune, "verbose"=verbose)


  ####################################################################################################
  # analyze.mix   MCMC function draws from posterior for LogN-LogS project.
  #               Methods: Gibbs sampler and Metropolis-Hastings algorithm
  #
  # Input: Y         = total observed photon counts
  #        niter     = total number of iterations per chain
  #        burnin    = burn-in number of iterations
  #        v.S       = value SD, tuning parameter for fluxes S to accept 20-60% of proposals
  #        v.theta   = value SD, tuning parameter for slope theta to accept 20-60% of proposals
  #        v.sm      = value SD, tuning parameter for Smin to accept 20-60% of proposals
  #        v.bp      = value SD, tuning parameter for bp to accept 20-60% of proposals
  #        alpha     = target number of sccessful trials in negbinom prior for N
  #        beta      = dispersion parameter in negbinom prior for N
  #        a         = vector of shape hyper-parameter(s) in gamma prior for theta(s)
  #        b         = vector of rate hyper-parameter(s) in gamma prior for theta(s)
  #        d         = vector if hyper-parameter(s) in dirichlet prior for mixture proportions (p_0,...,p_m)) (broken power-law)
  #        am        = shape hyper-parameter(s) in gamma prior for Smin
  #        bm        = rate hyper-parameter(s) in gamma prior for Smin
  #        mu        = mean parameter in normal proposal for break-points
  #        C         = varianca-covariance matrix of parameters in normal proposal for break-points
  #        gamma     = constant of transformation, energy per photon
  #        pble      = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        E.obs     = vector of exposureMap of observed sources
  #        L.obs     = vector of location of observed sources
  #        bg.obs    = background counts/pixel
  #        k.obs     = intensity/expected background photon counts, k=z*E=B*E/gamma
  #        prob.MH   = probability of how often to perform Metropolis-Hastings step in sampling theta 
  #        start.state = how to generate fluxes {"random"/"estimate"/"fixed"list }
  #        fixed.N     = Value/FALSE: debugging for N - no mcmc
  #        fixed.theta = Value/FALSE: allow for user specified theta - no mcmc
  #        fixed.S.obs = Value/FALSE: debugging for S.obs - no mcmc
  #        fixed.S.mis = Value/FALSE: debugging for S.mis - no mcmc
  #        fixed.Smin  = Value/FALSE: allow for user specified Smin  - no mcmc
  #        fixed.bp    = either NULL (no fixed break-points), or, a vector of length p specifying the break-points
  #        fixed.p     = Value/FALSE, mixture Pareto probabilities: allow for user specified theta - no mcmc
  #        pi          = R.object: a list("pi","theta") with detection probalities: pi(theta)
  #        fixed.S.pi  = (override) fix the value of S (used mainly for debugging)
  #        length.S    = length of grid of S, for integration of g(theta)
  #        nsamples    = number of MC samples to produce estimate of pi(theta,Smin) integral
  #        debug.pi    = (T/F) (used mainly for debugging)
  #        g.type      = type of g-function {"step","smooth","table"}
  #        g           = function, probability of observing a source
  #        area.conv.func = function, convert from off-axis angle to area [pix] - not used in this function, only for storage purposes
  #        delta.vec   = probability vector for normalization constant calculation
  #        model       = computing model for fluxes {"pareto", "bp", "mix"}
  #        met.alg       = list("Smin","bp") with "MH"/"MTM"/("HitRun") Metropolis method variations, used in Smin, bp sampling
  #        store_logPost = (T/F) whether to store estimate of log-posterior (without norm const)
  #        tune.iter   = How often to tune the proposal variance
  #        stop.tune   = when to stop tuning parameters
  #        save.every  = thin MCMC draws by this amount
  #        print.every = display code progress every kth iteration
  #        verbose     = (T/F) display progress of program
  #        save.progress      = manually save current draws to file
  #        save.progress.dir  = storage file for current draws
  #        compute.dic.extras = (T/F) compute additional statistics 
  #
  # Output:  list of :
  #             draws  = posterior samples of parameters (N, theta, Smin, bp, S.obs),
  #             prop.accept = accepted proportions for MH sampled parameters (S.obs, S.mis, theta, Smin, bp),
  #             draws.S.mis, draws.mix.I, inputs, eff.size, 
  #             norm.const, loglik, DIC, D_compt, D_compt_hat
  ####################################################################################################
  
  startTime <- proc.time()
  
  # Create pble object, if needed
  if (class(pble)=="character") {
    pble <- new.pble.table(effects.file=pble)
  } else if (class(pble)!="pble.basic" && class(pble)!="pble.table"){
    stop("'pble' could not be created.")
  }  
  
  if (burnin>niter)
    stop("'burnin' must be less than 'niter'")
  if (niter<1)
    stop("'niter' must be positive")
  
  Y.obs.tot <- Y
  n <- length(Y.obs.tot)           	            #number of observed sources
  m <- length(a)
  if (length(b) != m){
    stop("'a' and 'b' must be the same length")
  }
  
  use.mix <- FALSE
  use.bp  <- FALSE
  if (model=="bp"){
    use.bp <- TRUE
  } else if (model=="mix"){
    use.mix <- TRUE
  }
  
  # sanity check
  if(use.bp & use.mix){
    stop("Check parameter specification! Cannot be both bp and mix case.")
  }
  
  if(use.bp) {
    if(mu>0){
      stop("'mu' parameter for break-point prior must be negative.")
    }
  }
  
  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  
  if (verbose) {
    cat("Generating starting values.. \n") 
  } 
  
  ### Generate Smin, the minimum source flux, informative prior:
  if (any(!fixed.Smin)){     #some/all Smin values are not fixed and need to be generated
    if (is.list(start.state)){
      if (start.state$method=="fixed"){
        Smin <- start.state$Smin.value
      } else {stop("Unable to produce valid starting state 'fixed' for Smin due to wrong list format.")}
    } else if (start.state=="estimate"){
      S.obs <- ((Y.obs.tot-k.obs)*gamma)/E.obs     #vector of initial values: S.obs0
      subby <- (Y.obs.tot-k.obs)>2 
      tmp.Sobs <- S.obs[subby]
      if (length(tmp.Sobs)==0){
        warning("Unable to produce valid starting state 'estimate' for Smin due to no fluxes (no temporary S.obs).")
        tmp.Smin <- max(c(min(0.5*gamma/E.obs),min(S.obs)))
      } else {  
        tmp.Smin <- min(tmp.Sobs)
      }
      Smin <- tmp.Smin*0.95
    } else if (start.state=="random"){
      Smin <- rgamma(n=1,shape=am,rate=bm)  #initial value: Smin0
    }
    idx <- (fixed.Smin != FALSE)          # index of fixed values
    Smin[idx] <- fixed.Smin[idx]          # overwrite with fixed values
  } else {
    Smin  <- fixed.Smin
  }
  if(use.bp) {
    if (all(!!fixed.bp)) {
      if (all(Smin>fixed.bp)) {
        Smin <- min(fixed.bp)*0.90
      }
    }
  }
  
  ### Generate the breakpoints:
  # fixed.bp has 3 options:
  # Option 1 :: fixed.bp == NULL                 ==> No break-points
  # Option 2 :: fixed.bp == (vector of scalars)  ==> Break-points, fixed to user-specified values
  # Option 3 :: fixed.bp == FALSE                ==> Break-points, must be estimated
  
  if (!use.bp){
    # Option 1 :: No break-points
    bp <- NULL
  } else {
    if (any(is.logical(fixed.bp))){
      if (is.list(start.state)){
        if (start.state$method=="fixed"){
          bp <- start.state$bp.value
        } else {stop("Unable to produce valid starting state 'fixed' for bp due to wrong list format.")}
      } else {
        if (any(!fixed.bp)){  
          # Option 3 :: Break-points, must be estimated
          if (verbose) {
            cat("Generating K from prior.\n")
          }
          if (!fixed.Smin && length(C)==m) {
            K <- cumsum(exp(rnorm(n=m, mean=mu, sd=C)))
            Smin <- K[1]
            bp <- K[-1]
          } else {
            bp <- cumsum(exp(rnorm(n=m-1, mean=mu, sd=C)))
            bp <- sapply(bp,max,Smin*1.01) # require that Smin<bp.
          }      
        } else {
          stop("When specifying a break-point, you must use 'fixed.bp=value_you_want_to_fix_to'")
          # Use the user-specified break-points
        }
      } # END start-state is estimate/random
    } else {
      # Option 2 :: Break-points, fixed to user-specified values
      if (length(fixed.bp) != (m-1)){
        stop("When 'fixed.bp' is a vector of scalars, it must be of length 'length(a)-1' and 'length(b)-1'")
      }
      bp <- fixed.bp
    }
    if (verbose){
      cat("bp=\n"); print(bp)
    }
  } # END generating bp
    
  ### Generate total number of sources N:
  if (!fixed.N){
    prob.obs.source <- g(lambda=Smin*pble["pars"]$Emin/gamma,bg=pble["pars"]$Bmin,E=pble["pars"]$Emin,L=pble["pars"]$Lmin,g.type=g.type)
    
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
    if (is.list(start.state)){
      if (start.state$method=="fixed"){ 
        theta.t <- start.state$theta.value
      } else {stop("Unable to produce valid starting state 'fixed' for theta due to wrong list format.")}
      
    } else {
      theta.t <- rgamma(n=m,shape=a,rate=b)  #initial value: theta0
    }
    idx <- (fixed.theta != FALSE)          # index of fixed values
    theta.t[idx] <- fixed.theta[idx]       # overwrite with fixed values
  } else {
    theta.t <- fixed.theta
  }
  
  ### Generate mixture components:
  # Mixture-pareto has 2 options:
  # Option 1 :: d == NULL                        ==> No mixtures
  # Option 2 :: d == (vector of scalars)         ==> Mixtures, user-specified
  
  if (!use.mix){
    # Option 1 :: No mixtures    
    p.t     <- NULL
    I.idx.t <- NULL
    
  } else {
    
    # Option 2 :: Mixtures, total number fixed to user-specified values
    if (length(d) != (m)){
      stop("When specifying mixtures, 'd' is a vector of scalars, it must be of length 'length(a)' and 'length(b)'")
    }
    
    ### Generate mixture proportion(s):
    idx <- (fixed.p==FALSE) #index of non-fixed values    
    if (all(idx)){
      p.t <- rdirichlet(n=1,alpha=d)                     #initial value: (p_0,...,p_m)
    } else { #some/all p values are fixed
      ## Suppose (X1,X2,X3,X4,X5)~Dirichlet(d1,d2,d3,d4,d5), s.t. sum(X_i)=1. 
      ## Then, (X1,X2,X3,X4+X5)~Dirichlet(d1,d2,d3,d4+d5), s.t. sum(X_i)=1.
      ## And, (X1,X2,X3)*1/(1-(X4+X5)) ~ Dirichlet(d1,d2,d3)  indep of  Y=(X4+X5) ~ Beta(d4+d5,d1+d2+d3) 
      ## So, X1,X2,X3 | X4, X5  ~  (1-(X4+X5))*Z, where Z~Dirichlet(d1,d2,d3)  
      p.t <- fixed.p
      if (any(idx)) { #some p values are fixed
        p.tilde <- rdirichlet(n=1,alpha=d[idx])
        p.t[idx] <- p.tilde*(1-sum(p.t[fixed.p != FALSE]))          #initial value: (p_0,...,p_m)
      }
    }
  } # END generating mixture probabilities
  
  
  ### Generate source flux S.obs:
  # Generalizing to handle both break-points and untruncated mixtures
  if (verbose){
    cat('Debugging S.obs generation\n')
    cat('!fixed.S.obs:'); print(!fixed.S.obs)
    cat('!fixed.S.mis:'); print(!fixed.S.mis)
    cat('Smin:'); print(Smin)
    cat('theta.t:'); print(theta.t)
    cat('p.t:'); print(p.t)
  }
  if (any(!fixed.S.obs)){
    if (is.list(start.state)){
      if (start.state$method=="fixed"){
        S.obs.t <- start.state$S.obs.value
      } else {stop("Unable to produce valid starting state for Smin due to wrong list format.")}
      
    } else if (start.state=="estimate"){
      
      tmp <- (Y.obs.tot*gamma)/E.obs     #vector of initial values: S.obs0
      tmp[tmp<Smin] <- Smin*1.0001
      S.obs.t <- tmp
      
      if (verbose){
        cat("Initial values of S.obs.t, before patching:")
        print(S.obs.t)
      }
      
    } else if (start.state=="random"){
      
      if (use.bp){
        # Generate from broken power-law mixture of Truncated Pareto's:
        S.obs.t <- rbrokenpareto(n=n,x_min=Smin,k=theta.t,bp=bp,verbose=verbose)     #vector of initial values: S.obs0
        I.idx   <- NULL
      } else if (use.mix){
        # Generate from mixture of Pareto's:
        tmp     <- rmixpareto(n=n,p=p.t,k=theta.t,x_min=Smin,verbose=verbose)
        I.idx.t <- tmp$I.idx         #vector of initial values: I.idx(1,...,n)
        S.obs.t <- tmp$S             #vector of initial values: S.obs0
      } else {
        # Generate from single Pareto:
        S.obs.t <- rpareto(n=n,k=theta.t,x_min=Smin) #vector of initial values: S.obs0
        I.idx   <- NULL
      }
    }
    
    # Check starting states of the observed sources are all valid...
    S.obs.t <- patch.sobs(S.obs.t=S.obs.t, theta.t=theta.t, I.idx.t=I.idx.t[1:n], Smin=Smin, 
                          E.obs=E.obs, L.obs=L.obs, bg.obs=bg.obs, gamma=gamma, 
                          g.type=g.type, g=g, bp=bp, use.bp=use.bp, use.mix=use.mix, verbose=verbose)
    if (verbose) {
      cat("Done with patching Sobs. \n") 
      print(S.obs.t)
    }
    
  } else {
    S.obs.t <- fixed.S.obs
  }
  
  ### Generate standard background, effective area etc.:
  par <- sample.pble(pble,N.t-n)
  bg.mis <- par$B  #vector if initial values: bg.mis0, not saved during update
  L.mis <- par$L  #vector if initial values: L.mis0, not saved during update
  E.mis <- par$E  #vector if initial values: E.mis0, not saved during update
  if (verbose){
    cat("E.mis:",E.mis,"\n")
    cat("L.mis:",L.mis,"\n")
    cat("bg.mis:",bg.mis,"\n")
  }
  
  
  ### Generate source flux S.mis:    
  if (any(!fixed.S.mis)){
    if (N.t==n){     
      S.mis.t <- numeric(0)  
      I.idx   <- numeric(0)
      
    } else if (N.t>n){
      
      if (use.bp){
        # Generate from broken power-law mixture of Truncated Pareto's:
        S.mis.t <- rbrokenpareto(n=N.t-n,x_min=Smin,k=theta.t,bp=bp,verbose=verbose)     #vector of initial values: S.mis0
        I.idx   <- NULL
      } else if (use.mix){
        # Generate from mixture of Pareto's:
        tmp     <- rmixpareto(n=N.t-n,p=p.t,k=theta.t,x_min=Smin)
        I.idx.t <- c(I.idx.t, tmp$I.idx)         #UPDATE/ADD vector of initial values: I.idx(n+1,...,N) for S.mis component
        S.mis.t <- tmp$S             #vector of initial values: S.mis0
      } else {
        # Generate from single Pareto:
        S.mis.t <- rpareto(n=N.t-n,k=theta.t,x_min=Smin)  #vector of initial values: S.mis0
        I.idx   <- NULL
      }
      # Check starting states of the missing sources are all valid...
      S.mis.t <- patch.smis(S.mis.t=S.mis.t, theta.t=theta.t, I.idx.t=I.idx.t[(n+1):N.t], Smin=Smin, 
                            E.mis=E.mis, L.mis=L.mis, bg.mis=bg.mis, gamma=gamma, pble=pble,
                            g.type=g.type, g=g, bp=bp, use.bp=use.bp, use.mix=use.mix, verbose=verbose)
      if (verbose) {
        cat("Done with patching Smis. \n") 
      }
    }
  } else {
    S.mis.t <- fixed.S.mis
  }
  
  if (length(v.S)<1){
    stop("v.S (SD parameter for MH proposal for flux S) has not been specified as a scalar or is unknown. Please original check settings")
  } else if (length(v.S)==1){
    v.S <- rep(v.S,n)              #v.S = tuning parameter (vector) to accept 20-60% of proposals
  } else {
    if (length(v.S)!=n) {
      v.S <- rep(v.S[1],n)
    }
  } 
  if (length(v.theta)<1){
    stop("v.theta (SD parameter for MH proposal for theta) has not been specified as a scalar or is unknown. Please original check settings")
  } else if (length(v.theta)==1) {
    v.theta <- rep(v.theta,m)      #v.theta = tuning parameter (vector) to accept 20-60% of proposals
  } else {
    if (length(v.theta)!=m) {
      v.theta <- rep(v.theta[1],m)
    }
  } 
  if (length(v.bp)<1){
    stop("v.bp (SD parameter for MH proposal for preak-point) has not been specified as a scalar or is unknown. Please original check settings")
  } else if (length(v.bp)==1){
    v.bp <- rep(v.bp,m-1)              #v.bp = tuning parameter (vector) to accept 20-60% of proposals
  } else {
    if (length(v.bp)!=m-1) {
      v.bp <- rep(v.bp[1],m-1)
    }
  } 
  
  prop.ct.S.obs <- numeric(n)    #count of accepted proposals for S.obs
  prop.ct.theta <- numeric(m)    #count of accepted proposals for theta
  prop.ct.bp    <- numeric(m-1)    #count of accepted proposals for bp
  prop.ct.Smin  <- numeric(1)    #count of accepted proposals for Smin
  prop.ct.S.mis <- numeric(N.t*2)        #count of accepted proposals for S.mis, larger than needed
  prop.state <- list("prop.ct.S.obs"=prop.ct.S.obs,"prop.ct.theta"=prop.ct.theta,
                     "prop.ct.bp"=prop.ct.bp, "prop.ct.Smin"=prop.ct.Smin)   #current state of proposals
  
  ### Set-up storage matrix for 'draws'
  total.samples <- floor((niter-burnin)/save.every)
  len.draws <- 1+m # for N and theta
  name.draws <- c("N",ppaste("theta.",1:m))
  # Add other parameters: tau.1=Smin, tau.j=bp, p.mix, S.obs.j
  if (use.mix){
    len.draws <- len.draws+m # for p.t
    name.draws <- c(name.draws, ppaste("p.",1:m))
  }
  if (!fixed.Smin) {
    len.draws <- len.draws+1 # for Smin
    name.draws <- c(name.draws, "tau.1")
  }
  if (use.bp){
    if (any(!fixed.bp)) {
      len.draws <- len.draws+m-1 # for p.t
      name.draws <- c(name.draws, ppaste("tau.",2:m))
    }
  }
  len.draws <- len.draws+n # for S.obs
  name.draws <- c(name.draws, ppaste("S.obs.",1:n))
  
  draws <- mcmc(matrix(numeric(0),nrow=total.samples,ncol=len.draws))
  varnames(draws) <- name.draws
 
  # Add log-posterior results
  if (store_logPost) {
    draws <- mcmc(cbind(draws,"log.Poster"=NA))
  }
  draws.S.mis <- vector("list",total.samples)  
  draws.mix.I <- vector("list",total.samples) 
  
  loglik <- list("loglik"=numeric(total.samples))
  components <- matrix(NA,ncol=(total.samples),nrow=n+1)
  rownames(components) <- c("Bin",ppaste("Pois_y",1:n))
  colnames(components) <- ppaste("Pars_",1:(total.samples))
  
  if (verbose) {
    cat("Generating draws.. \n") 
  } 
  
  save.i <- 1
  for (iter in 1:niter) {
    
    update <- update.lns(S.obs.t=S.obs.t, theta.t=theta.t, N.t=N.t, I.idx.t=I.idx.t, p.t=p.t, n=n, Y.obs.tot=Y.obs.tot, 
                         S.mis.t=S.mis.t, v.S=v.S, v.theta=v.theta, v.bp=v.bp, v.sm=v.sm, niter=niter,
                         prop.ct.S.obs=prop.ct.S.obs, prop.ct.theta=prop.ct.theta, prop.ct.S.mis=prop.ct.S.mis, prop.ct.bp=prop.ct.bp, prop.ct.Smin=prop.ct.Smin,
                         a=a, b=b, C=C, mu=mu, d=d, alpha=alpha, beta=beta, Smin=Smin, gamma=gamma, pble=pble,
                         E.obs=E.obs, L.obs=L.obs, bg.obs=bg.obs, k.obs=k.obs, length.S=length.S, 
                         fixed.S.pi=fixed.S.pi, debug.pi=debug.pi, nsamples=nsamples,
                         g.type=g.type, g=g, pi=pi, bp=bp, prob.MH=prob.MH,am=am,bm=bm,
                         fixed.N=fixed.N, fixed.theta=fixed.theta, fixed.p=fixed.p, fixed.bp=fixed.bp,
                         fixed.S.obs=fixed.S.obs, fixed.S.mis=fixed.S.mis, fixed.Smin=fixed.Smin, 
                         met.alg=met.alg, store_logPost=store_logPost,
                         use.bp=use.bp, use.mix=use.mix, verbose=verbose)
    
    # Update all variables
    N.t     <- update$draws$N.t
    theta.t <- update$draws$theta.t
    S.obs.t <- update$draws$S.obs.t
    
    if (any(!fixed.Smin)){
      Smin <- update$draws$Smin
    }
    if (use.bp){
      if(any(!fixed.bp)){ # bp=random
        bp   <- update$draws$bp  
      }  
    } else if(use.mix){
      p.t    <- update$draws$p.t 
    } 
    
    prop.ct.S.obs <- update$prop.ct.S.obs
    prop.ct.theta <- update$prop.ct.theta
    prop.ct.bp    <- update$prop.ct.bp
    prop.ct.Smin  <- update$prop.ct.Smin
    
    S.mis.t <- update$draws.S.mis
    I.idx.t <- update$draws.mix.I
    prop.ct.S.mis <- update$prop.ct.S.mis
    
    #Tune v.S, v.theta parameters every tune.iter iterations
    if (iter%%tune.iter == 0 & iter <= stop.tune){
      tune <- tune.v(v.S=v.S, v.theta=v.theta, v.bp=v.bp, v.sm=v.sm, prop.ct.S.obs=prop.ct.S.obs, prop.ct.theta=prop.ct.theta, 
                     prop.ct.bp=prop.ct.bp, prop.ct.Smin=prop.ct.Smin, prop.state=prop.state, 
                     tune.iter=tune.iter, N.t=N.t, n=n, verbose=verbose)
      prop.state <- tune$prop.state
      v.S     <- tune$v.S
      v.theta <- tune$v.theta
      v.bp    <- tune$v.bp
      v.sm    <- tune$v.sm
    }    
    
    # Store draws after burnin (Automatically remove burnin) (Automatically thin draws)
    if (iter>burnin) {  
      if ((iter-burnin) %% save.every == 0) {
        draws[(iter-burnin)/save.every,] <- unlist(update$draws)
        draws.S.mis[[(iter-burnin)/save.every]] <- update$draws.S.mis
        draws.mix.I[[(iter-burnin)/save.every]] <- update$draws.mix.I
        
        # Evaluate the log( pr(Data|param(i)) ) in order to estimate of log of normalizing constant
        pi.value <- pi.theta.get(pi=pi, theta=theta.t, p.t=p.t, bp=bp, Smin=Smin, gamma=gamma, 
                                 pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi, nsamples=nsamples,
                                 g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources  
        lambda.obs <- S.obs.t*E.obs/gamma
        pr.n <- dbinom(x=n, size=N.t, prob=pi.value, log=TRUE)
        pr.y <- dpois(x=Y.obs.tot, lambda=lambda.obs+k.obs, log=TRUE)
        components[1,(iter-burnin)/save.every] <- pr.n
        components[1+(1:n),(iter-burnin)/save.every] <- pr.y
        loglik$loglik[(iter-burnin)/save.every] <- pr.n + sum(pr.y) 
        
      }
    }
    
    if (iter%%print.every == 0){ 
      
      cat("N, theta, Smin, and S[1] are:   ")
      print(c(update$draws$N.t,update$draws$theta.t,update$draws$Smin,update$draws$S.obs.t[1]))
      
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
  
  if (compute.dic.extras){
    
    # Compute log of normalizing constant (predictive score)
    cat("Computing log(normalizing constant)...\n")
    nn <- total.samples
    adj.const <- min(loglik$loglik)   #constant to help numerical balance
    temp <- exp(adj.const - loglik$loglik)
    norm.const <- (adj.const) + log(nn) - log(sum(temp))
    
    #   # Compute log of normalizing constant (predictive score) for iter=[1:Iter], ending at every iteration
    #   nc.vec <- rep(NA,nn)
    #   for (i in 1:nn){
    #     tmp.ll <- loglik[1:i]
    #     adj.const <- min(tmp.ll)   #constant to help numerical balance
    #     temp <- exp(adj.const - tmp.ll)
    #     nc.vec[i] <- (adj.const) + log(i) - log(sum(temp))
    #     if (i %% 1000 == 0){
    #       cat(i,"\n")
    #     }
    #   }
    #   nc.vec.iters <- nc.vec
    #   norm.const <- nc.vec[nn]
    
    #   # Method #3
    #   TOL <- exp(-400)
    #   n.tries.max <- 1000
    #   
    #   nc.delta3 <- rep(NA, length(delta.vec))
    #   for(j in 1:length(delta.vec)){
    #     delta <- delta.vec[j]
    #     # Iteratively compute the normalizing constant
    #     tmp.l <- exp(loglik)
    #     n.tries <- 0
    #     diff <- 1
    #     # start state
    #     nc0 <- exp(nc.vec[nn])
    #     nc <- rep(NA,n.tries.max)
    #     nc[1] <- nc0
    #     c <- delta*nn/(1-delta)
    #     while (diff>TOL & n.tries<n.tries.max){
    #       n.tries <- n.tries+1
    #       # update equation   
    #       denom <- (delta*nc0 + (1-delta)*tmp.l)
    #       nc[n.tries+1] <- (c+sum( tmp.l/denom ))/(c*nc0+sum( 1/denom ))
    #       diff <- abs(nc[n.tries+1]-nc0)
    #       # update state
    #       nc0 <- nc[n.tries+1]     
    #     }  
    #     nc <- log(nc)
    #     if(verbose){
    #       cat("log(norm.const.)=\n")
    #       print(nc[n.tries+1])
    #     }
    #     nc.delta3[j] <- nc[n.tries+1]
    #   }
    #   nc.vec <- nc.delta3
    
    ############################################
    ### Compute DIC: Deviance Information Criteria
    cat("Computing DIC...\n")
    
    stats <- summary(draws,quantiles=0.5)
    draws.mean <- stats$statistics[,"Mean"]
    draws.median <- stats$quantiles
    #mode:
    draws.mode <- rep(NA,ncol(draws))
    theta.id <- grep("theta",names(draws.mean))
    sobs.id <- grep("S.obs",names(draws.mean))
    p.id <- grep("p.",names(draws.mean))   #either a vector or integer(0) 
    for(mm in 1:ncol(draws)){
      if (any(mm==sobs.id)){
        rounding <- (-1)*floor(log(min(fixed.Smin)/100,base=10))
      } else {
        rounding <- 3
      }
      draws.mode[mm] <- as.numeric(names(sort(-table(round(draws[,mm],rounding))))[1])
    }
    
    components_hat <- matrix(NA,ncol=3,nrow=1+n)
    rownames(components_hat) <- rownames(components)
    colnames(components_hat) <- ppaste("Pars_",c("mean","median","mode"))
    # Evaluate loglik.at.mean
    N.m <- round(draws.mean[1])
    theta.m <- draws.mean[theta.id]
    S.obs.m <- draws.mean[sobs.id]
    if (model=="mix"){ p.m <- draws.mean[p.id] } else {p.m <- NULL}
    pi.value <- pi.theta.get(pi=pi, theta=theta.m, p.t=p.m, bp=bp, Smin=Smin, gamma=gamma,  pble=pble, 
                             length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi, g.type=g.type, g=g, 
                             nsamples=nsamples, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources  
    lambda.obs <- S.obs.m*E.obs/gamma
    pr.n <- dbinom(x=n, size=N.m, prob=pi.value, log=TRUE)
    pr.y <- dpois(x=Y.obs.tot, lambda=lambda.obs+k.obs, log=TRUE)
    components_hat[,"Pars_mean"] <- c(pr.n,pr.y)
    loglik$loglik.at.mean <- pr.n + sum(pr.y)   
    # Evaluate loglik.at.median
    N.m <- round(draws.median[1])
    theta.m <- draws.median[theta.id]
    S.obs.m <- draws.median[sobs.id]
    if (model=="mix"){ p.m <- draws.median[p.id] } else {p.m <- NULL}
    pi.value <- pi.theta.get(pi=pi, theta=theta.m, p.t=p.m, bp=bp, Smin=Smin, gamma=gamma,  pble=pble, 
                             length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi, g.type=g.type, g=g, 
                             nsamples=nsamples, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources  
    lambda.obs <- S.obs.m*E.obs/gamma
    pr.n <- dbinom(x=n, size=N.m, prob=pi.value, log=TRUE)
    pr.y <- dpois(x=Y.obs.tot, lambda=lambda.obs+k.obs, log=TRUE)
    components_hat[,"Pars_median"] <- c(pr.n,pr.y)
    loglik$loglik.at.median <- pr.n + sum(pr.y) 
    # Evaluate loglik.at.mode
    N.m <- round(draws.mode[1])
    theta.m <- draws.mode[theta.id]
    S.obs.m <- draws.mode[sobs.id]
    if (model=="mix"){ p.m <- draws.mode[p.id] } else {p.m <- NULL}
    pi.value <- pi.theta.get(pi=pi, theta=theta.m, p.t=p.m, bp=bp, Smin=Smin, gamma=gamma,  pble=pble, 
                             length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi, g.type=g.type, g=g, 
                             nsamples=nsamples, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources  
    lambda.obs <- S.obs.m*E.obs/gamma
    pr.n <- dbinom(x=n, size=N.m, prob=pi.value, log=TRUE)
    pr.y <- dpois(x=Y.obs.tot, lambda=lambda.obs+k.obs, log=TRUE)
    components_hat[,"Pars_mode"] <- c(pr.n,pr.y)
    loglik$loglik.at.mode <- pr.n + sum(pr.y)  
    
    D_bar_k <- (-2)*mean(loglik$loglik)
    D_theta_bar <- (-2)*loglik$loglik.at.mean
    D_theta_median <- (-2)*loglik$loglik.at.median
    D_theta_mode <- (-2)*loglik$loglik.at.mode
    
    DIC_mean <- 2*D_bar_k - D_theta_bar
    DIC_median <- 2*D_bar_k - D_theta_median
    DIC_mode <- 2*D_bar_k - D_theta_mode
    DIC <- list("D"     =list("D_bar_k"=D_bar_k,"D_all"=(-2)*loglik$loglik),
                "mean"  =list("DIC"=DIC_mean,  "penalty"=D_bar_k-D_theta_bar,   "D_theta_hat"=D_theta_bar),
                "median"=list("DIC"=DIC_median,"penalty"=D_bar_k-D_theta_median,"D_theta_hat"=D_theta_median),
                "mode"  =list("DIC"=DIC_mode,  "penalty"=D_bar_k-D_theta_mode,  "D_theta_hat"=D_theta_mode))
  } else {
    DIC <- NULL
    norm.const <- NULL
    #loglik <- NULL
    D_compt <- NULL
    D_compt_hat <- NULL
    components_hat <- NULL
  } # END compute DIC and extras
  
  ############################################
  # Compute effective sample size
  eff.size <- effectiveSize(draws)
  ############################################
  
  spa <- round(prop.ct.S.obs/niter,2)*100
  tpa <- round(prop.ct.theta/niter,2)*100
  bpa <- round(prop.ct.bp/niter,2)*100
  mpa <- round(prop.ct.Smin/niter,2)*100
  prop.accept <- list("v"=c(v.S,v.theta,v.sm,v.bp),
                      "S.prop.accept"=spa,"theta.prop.accept"=tpa,"Smin.prop.accept"=mpa,"bp.prop.accept"=bpa,
                      "S.mis.ct.accept"=prop.ct.S.mis)
  
  return(list("draws"=draws,"prop.accept"=prop.accept,
              "draws.S.mis"=draws.S.mis,"draws.mix.I"=draws.mix.I,
              "inputs"=inputs,"eff.size"=eff.size,"norm.const"=norm.const,
              "loglik"=loglik,
              "DIC"=DIC,"D_compt"=(-2)*components,"D_compt_hat"=(-2)*components_hat))
}

analyze.mix <- cmpfun(analyze.mix)
