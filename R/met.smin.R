"met.smin" <- function(Smin.t, bp.t, N.t, n, theta.t, S.obs.t, am,bm, v.sm,
                       gamma, E, g, nsamples, sigma, 
                       verbose=FALSE){
  
  ####################################################################################################
  # met.smin   Metropolis-Hasting step to draw sample of Smin, minimum threshold for S, in mixture Pareto case.
  #
  # Input: Smin.t  = previous iteration of Smin, minimum flux
  #        N.t     = total number of sources, current iteration
  #        theta.t = current iteration of theta, lns slope
  #        S.obs.t = current vector of observed sources of flux
  #        S.mis.t = current vector of missing sources of flux
  #        I.idx.t = current iteration of mixture neighborhoods 
  #        n.t     = number of observed sources
  #        am      = shape hyper-parameter(s) in gamma prior for Smin
  #        bm      = rate hyper-parameter(s) in gamma prior for Smin
  #        v.sm    = value SD, tuning parameter for Smin to accept 20-60% of proposals
  #        gamma   = constant of transformation, energy per photon
  #        pble    = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        length.S = length of grid of S, for integration of g(theta)
  #        g.type   = type of g-function {"step","smooth","table"}
  #        g        = function, probability of observing a source
  #        bp.t     = current state of bp
  #        pi         = R.object: a list("pi","theta") with detection probalities: pi(theta)
  #        fixed.S.pi = (override) fix the value of S (used mainly for debugging)
  #        nsamples   = number of MC samples to produce estimate of pi(theta,Smin) integral
  #        debug.pi = (T/F) (used mainly for debugging)
  #        use.bp   = (T/F) break-point pareto version
  #        use.mix  = (T/F) mixture pareto version
  #        algorithm.type = value:(1)/(2)  #1=uses Scom w/o binomial term,  2=uses Sobs
  #        met.alg  = "MH"/"MTM"/("HitRun") Metropolis method variations, used in Smin, bp sampling
  #        kk       = number of tries to propose new draws in Multiple-Try-Metropolis
  #        prob.MH     = probability of how often to perform Metropolis-Hastings step in sampling Smin  
  #        n.tries.max = max number of tries to propose new vectors of S.mis.t
  #        verbose = (T/F) display progress of program
  #
  # Output: Smin.t  = sample draw of vector of Smin
  #         idx     = indices of accepted proposals
  ####################################################################################################
  
  
  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  
  upperBound <- min(S.obs.t)
  S <- S.obs.t
  
  # current state
  eta.t <- log(Smin.t)
  pi.value.curr <- pi.theta.get(theta=theta.t, Smin=Smin.t, bp=bp.t, gamma=gamma, 
                                E=E, nsamples=nsamples, g=g, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
  jacobian.curr <- eta.t
  p1.curr <- jacobian.curr + dgamma(x=Smin.t, shape=n*theta.t[1] + am, rate=bm, log=TRUE)
  p2.curr <- (N.t-n)*log(1-pi.value.curr)
  if(N.t==n){
    p2.curr <- 0
  }
  p.curr <- p1.curr + p2.curr
  
  # proposal state
  eta.prop <- rnorm(n=1, mean=eta.t, sd=v.sm)    #v.sm = tuning parameter (vector) to accept 20-60% of proposals
  Smin.prop <- exp(eta.prop)
  
  if (0 < Smin.prop && all(Smin.prop < bp.t) && Smin.prop < upperBound) {
    pi.value.prop <- pi.theta.get(theta=theta.t, Smin=Smin.prop, bp=bp.t, gamma=gamma, 
                                  E=E, nsamples=nsamples, g=g, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
    jacobian.prop <- eta.prop
    p1.prop <- jacobian.prop + sum(dgamma(x=Smin.prop, shape=n*theta.t[1] + am, rate=bm, log=TRUE))
    p2.prop <- (N.t-n)*log(1-pi.value.prop)
    if(N.t==n){
      p2.prop <- 0
    }
    p.prop <- p1.prop + p2.prop
  } else {
    p.prop <- -Inf
    p1.prop <- p2.prop <- jacobian.prop <- NA
  }
  
  if(verbose){
    cat("--- Inside updating step for S.obs.t: --- \n")
    cat("Smin.t:     #before update\n"); print(Smin.t)
    cat("Smin.prop:\n"); print(Smin.prop)
  }
    
  # Compute ratio of densities for MH
  log.alpha <- p.prop - p.curr            #compare proposal to current
  log.u <- log(runif(1))                #random vector U(0,1)
  got.draw.idx <- (log.u < log.alpha)            #accepted proposal indicator 
  got.draw.idx[is.na(got.draw.idx)] <- FALSE
  Smin.t[got.draw.idx] <- Smin.prop[got.draw.idx]       #update accepted proposals, iteration t=iter+1
  
  if(verbose){
    # Print above results (error checking)
    cat("Smin.t:     #after update\n"); print(Smin.t)
    cat("Smin.prop:\n"); print(Smin.prop)
    cat(paste("v.sm :Smin = ")); print(v.sm)
    cat("dgamma.curr:\n"); print(p1.curr)
    cat("dgamma.prop:\n"); print(p1.prop)
    cat("log.pi.curr:\n"); print(p2.curr)
    cat("log.pi.prop:\n"); print(p2.prop)
    cat("p.curr:\n"); print(p.curr)
    cat("p.prop:\n"); print(p.prop)
    cat("log.alpha:\n"); print(log.alpha)
    cat("log.u:\n"); print(log.u)
    cat("got.draw.idx:\n"); print(got.draw.idx)
  }
  
  return(list('Smin.t'=Smin.t,"idx"=got.draw.idx))
}

met.smin <- cmpfun(met.smin)
