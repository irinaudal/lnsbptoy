"update.lns" <- function(S.obs.t, theta.t, Smin.t, bp.t, N.t, n, Y.obs.tot, S.mis.t, v.so, v.th, v.sm, v.bp, niter,
                         prop.ct.S.obs, prop.ct.S.mis,  prop.ct.theta, prop.ct.Smin, prop.ct.bp,
                         alpha, beta, a, b, am, bm, C, mu, gamma, E, g, nsamples, sigma,
                         fixed.N=FALSE, fixed.theta=FALSE, fixed.S.obs=FALSE, fixed.S.mis=FALSE, 
                         fixed.Smin=FALSE, fixed.bp=NULL, store_logPost=FALSE,
                         verbose=FALSE){
  
  
  ####################################################################################################
  # update.lns   MCMC function draws from posterior for LogN-LogS project.
  #                 Generate Y.obs.src.t, N.t, theta.t, S.obs.t; auto-tune SD parameters
  #                 Methods: Gibbs sampler and Metropolis-Hastings algorithm
  #
  # Input: S.obs.t = previous value of observed flux S.obs.t
  #        theta.t = previous value of theta.t
  #        Smin.t  = minimum flux the sources can be detected to, hyper-parameter in pareto prior for S
  #        bp.t    = either NULL (no fixed break-points), or, a vector of length p specifying the break-points
  #        N.t     = total number of sources, N.com
  #        n       = number of observed sources
  #        Y.obs.tot = total observed photon counts
  #        S.mis.t = NULL/flux of missing sources to be updated
  #        v.so     = value SD, tuning parameter to accept 20-60% of proposals of S.obs
  #        v.th     = value SD, tuning parameter to accept 20-60% of proposals of theta
  #        v.sm     = value SD, tuning parameter to accept 20-60% of proposals of Smin
  #        v.bp     = value SD, tuning parameter to accept 20-60% of proposals of eta(bp)
  #        niter   = total number of iterations per chain
  #        prop.ct.S.obs = number of accepted proposals for S.obs
  #        prop.ct.theta = number of accepted proposals for theta
  #        prop.ct.Smin = number of accepted proposals for Smin
  #        prop.ct.bp = number of accepted proposals for bp
  #        prop.ct.S.mis = number of accepted proposals for S.mis
  #        alpha   = target number of sccessful trials in negbinom prior for N
  #        beta    = dispersion parameter in negbinom prior for N
  #        a       = shape hyper-parameter in gamma prior for theta
  #        b       = rate hyper-parameter in gamma prior for theta
  #        am      = shape hyper-parameter(s) in gamma prior for Smin
  #        bm      = rate hyper-parameter(s) in gamma prior for Smin
  #        mu        = mean parameter in normal proposal for break-points
  #        C         = varianca-covariance matrix of parameters in normal proposal for break-points
  #        gamma   = constant of transformation, energy per photon
  #        E       = vector of exposureMap of observed sources
  #        nsamples   = number of MC samples to produce estimate of pi(theta,Smin) integral
  #        fixed.N     = Value/FALSE: debugging for N - no mcmc
  #        fixed.theta = Value/FALSE: allow for user specified theta - no mcmc
  #        fixed.S.obs = Value/FALSE: debugging for S.obs - no mcmc
  #        fixed.S.mis = Value/FALSE: debugging for S.mis - no mcmc
  #        fixed.Smin  = Value/FALSE: allow for user specified Smin  - no mcmc
  #        fixed.bp    = NULL/Value/FALSE, break-points: allow for user specified bp - no mcmc  
  #        g       = function, probability of observing source
  #        store_logPost = (T/F) whether to store estimate of log-posterior (without norm const)
  #        verbose = (T/F) display progress of program
  #
  # Output: draws  = posterior samples of parameters (N, theta, S.obs, Smin, bp, p)
  #         draws.S.mis = posterior samples of parameters S.mis
  #         prop.ct.S.obs = number of accepted proposals for S.obs
  #         prop.ct.theta = number of accepted proposals for theta
  #         prop.ct.bp    = number of accepted proposals for bp
  #         prop.ct.Smin  = number of accepted proposals for Smin
  #         prop.ct.S.mis = number of accepted proposals for S.mis
  ####################################################################################################
  
  if (verbose) {
    cat("Update step of N.t, theta.t, S.obs.t, Y.obs.src.t, and S.mis.t, I.idx.t, p.t, Smin.t \n") 
  } 
    
  m <- length(a)
  
  ### Sample S.obs  
  if (any(!fixed.S.obs)){
    
    # Metropolis-Hastings, sample S.obs               
    if (verbose) {
      cat("_____S.obs.t updating step (MH)_____\n") 
    }        
    met <- met.sobs(S.obs.t, theta.t=theta.t, n=n, Y.obs.tot=Y.obs.tot, Smin.t=Smin.t, bp.t=bp.t,
                    v.so=v.so, E=E, gamma=gamma, g=g, verbose=verbose)
    S.obs.t <- met$S.obs.t                        #next value of S.obs at iteration t=iter+1
    prop.ct.S.obs[met$idx] <- prop.ct.S.obs[met$idx]+1    #count accepted proposals 
    
  } # END S.obs update
  
  
  ### Rejection Sampling, sample S.mis
  if (any(!fixed.S.mis)){
    
    if (verbose) {
      cat("_____S.mis.t updating step (Rej.Sampl)_____\n") 
    }  
    if(!is.null(S.mis.t)){
      rej <- rej.smis(S.mis.t=S.mis.t, theta.t=theta.t, N.t=N.t, Smin.t=Smin.t, bp.t=bp.t,
                      n=n, E=E, gamma=gamma, g=g, verbose=verbose)            
      S.mis.t <- rej$S.mis.t                        #next value of S.mis at iteration t=iter+1
      prop.ct.S.mis[rej$idx] <- prop.ct.S.mis[rej$idx]+1    #count accepted proposals 
      
    } else{
      S.mis.t <- NULL
      prop.ct.S.mis <- NULL
      if (verbose)
        cat("S.mis.t = NULL, because g(S,E,L)=1 for all.\n")
    }
    
  } # END S.mis update

  ### Sample theta(s)
  if (any(!fixed.theta)){   #some/all theta values are not fixed and need to be generated
    
    #Rejection Sampling, &/or Metropolis-Hastings, sample theta.t
    if (verbose) {
      cat("_____theta.t updating step (Rej.Sampl/MH)_____\n") 
    }  
    rej <- met.theta(theta.t=theta.t, S.obs.t=S.obs.t, N.t=N.t, n=n, Smin.t=Smin.t, bp.t=bp.t,
                     a=a, b=b, v.th=v.th, fixed.theta=fixed.theta, sigma=sigma,
                     gamma=gamma, E=E, g=g, nsamples=nsamples, verbose=verbose)  
    theta.t <- rej$theta.t                                #next value of theta at iteration t=iter+1
    prop.ct.theta[rej$idx] <- prop.ct.theta[rej$idx]+1    #count accepted proposals 
    
    if (verbose) {
      cat("theta.t after update:\n"); print(theta.t) 
    }
    
  } # END theta update
  
  
  ### Exact, sample Smin.t, minimum threshold of pareto-s
  # if use.bp model, generate Smin.t separately from bp, i.e. conditional on bp.
  if (!fixed.Smin) {    
    if (verbose) {
      cat("_____Smin.t updating step (Rej.Sampl/MH/MTM)_____\n") 
    }
    smp <- met.smin(Smin.t=Smin.t, N.t=N.t, n=n, theta.t=theta.t, bp.t=bp.t, S.obs.t=S.obs.t,
                    am=am, bm=bm, v.sm=v.sm, sigma=sigma,
                    gamma=gamma, E=E, g=g, nsamples=nsamples, verbose=verbose)      
    Smin.t <- smp$Smin.t                        #next value of Smin.t at iteration t=iter+1
    prop.ct.Smin[smp$idx] <- prop.ct.Smin[smp$idx]+1 # count accepted proporals
    
    if (verbose) {
      cat("Smin.t after update:\n"); print(Smin.t) 
    }
    
  } # END Smin.t update
  
  ### Sample bp.t 
  if (!is.null(fixed.bp)) {
    if (any(!fixed.bp)) { 
      # Metropolis-Hastings, sample K.t bpeakpoints, without Smin.t.              
      if (verbose) {
        cat("_____bp updating step (MH/MTM)_____\n") 
      }       
      # Generate bp.t only: Smin.t is fixed at the current state
      met <- met.bp(bp.t=bp.t, Smin.t=Smin.t, S.obs.t=S.obs.t,  N.t=N.t, n=n, theta.t=theta.t, 
                    C=C, mu=mu, v.bp=v.bp, fixed.bp=fixed.bp, sigma=sigma,
                    gamma=gamma, E=E, g=g, nsamples=nsamples, verbose=verbose)      
      bp.t   <- met$bp.t
      prop.ct.bp[met$idx] <- prop.ct.bp[met$idx]+1    #count accepted proposals      
      
      if (verbose) {
        cat("bp.t after update:\n"); print(bp.t) 
      }
    }
  } # END fixed.bp update
  
  ### Exact via numerical integration, sample N
  if (!fixed.N){
    
    if (verbose) {
      cat("_____N.t updating step (Num.Integr)_____\n") 
    }   
    #      ni <- numint.N(N.t, theta.t=theta.t, n=n, alpha=alpha, beta=beta, pi=pi, verbose=verbose)
    ni <- numint.N(N.t, theta.t=theta.t, n=n, Smin.t=Smin.t, bp.t=bp.t,
                   alpha=alpha, beta=beta, sigma=sigma,
                   gamma=gamma, E=E, g=g, nsamples=nsamples, verbose=verbose)      
    N.t <- ni$N.t                        #next value of N at iteration t=iter+1
    
    if (verbose) {
      cat("N.t after update:\n"); print(N.t) 
    }
    
  } # END N update
  
  
  ### Store the current samples to 'draws'    
  if (any(!fixed.Smin)){ 
    if (verbose){
      cat("Bp case: Smin.t about to be stored to draws:\n"); print(Smin.t)
    }
    if (is.null(fixed.bp)) {
      draws <- list("N.t"=N.t, "theta.t"=theta.t, "Smin.t"=Smin.t, "S.obs.t"=S.obs.t)      
    } else if (any(!fixed.bp)){
      if (verbose){
        cat("Bp case: bp.t about to be stored to draws:\n"); print(bp.t)
      }
      draws <- list("N.t"=N.t, "theta.t"=theta.t, "Smin.t"=Smin.t, "bp.t"=bp.t, "S.obs.t"=S.obs.t)      
    } else {
      draws <- list("N.t"=N.t, "theta.t"=theta.t, "Smin.t"=Smin.t, "S.obs.t"=S.obs.t)      
    }
  } else{
    if (is.null(fixed.bp)) {
      draws <- list("N.t"=N.t, "theta.t"=theta.t, "S.obs.t"=S.obs.t)      
    } else if (any(!fixed.bp)){
      draws <- list("N.t"=N.t, "theta.t"=theta.t, "bp.t"=bp.t, "S.obs.t"=S.obs.t) 
    } else {
      draws <- list("N.t"=N.t, "theta.t"=theta.t, "S.obs.t"=S.obs.t) 
    }
  }

  if (store_logPost) {
    # Evaluate log(posterior)
    # works for cage: g=1, pi=1, dim(bp)=1, Smin~gamma, g.type="smooth"
    pi.value <- pi.theta.get(theta=theta.t, Smin=Smin.t, bp=bp.t, gamma=gamma, E=E,
                             nsamples=nsamples, sigma=sigma, g=g, verbose=verbose>1)     #marginal prob. of observing sources  
    pi.const <- ifelse(N.t==n,0,(N.t-n)*log(1.0-pi.value))
    p.N   <- lgamma(N.t+alpha)-lgamma(N.t-n+1) + pi.const - N.t*log(1+beta)
    p.th  <- sum(dgamma(x=theta.t, shape=a, rate=b, log=TRUE)) 
    p.Sm  <- dgamma(x=Smin.t, shape=am, rate=bm, log=TRUE)
    p.bp  <- dnorm(x=log(bp.t-Smin.t), mean=mu, sd=C, log=TRUE)
    lambda <- S.obs.t*E/gamma
    gp    <- log(g(lambda=lambda))
    p.S   <- dbrokenpareto(S.obs.t, x_min=Smin.t, k=theta.t, bp=bp.t, log=TRUE)
    p.Yt  <- dpois(x=Y.obs.tot, lambda=lambda, log=TRUE)
    
    terms <- list(p.N=p.N,p.th=p.th,p.Sm=p.Sm,p.bp=p.bp,p.gp=gp,p.S=p.S,p.Yt=p.Yt)
    draws$logPost <- sum(unlist(terms))

    if(verbose) {
      cat("log.Poster.values"); print(terms)
    }
  }
  
  draws.S.mis <- S.mis.t
  
  return(list("draws"=draws,"prop.ct.S.obs"=prop.ct.S.obs,"prop.ct.theta"=prop.ct.theta,
              "prop.ct.Smin"=prop.ct.Smin,"prop.ct.bp"=prop.ct.bp,
              "draws.S.mis"=draws.S.mis,"prop.ct.S.mis"=prop.ct.S.mis
              ))
}

update.lns <- cmpfun(update.lns)