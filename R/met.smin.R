"met.smin" <- function(Smin.t, bp.t, N.t, n, theta.t, S.obs.t, am,bm, v.sm,
                       gamma, E, g, nsamples, pi, sigma, 
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
                                E=E, nsamples=nsamples, g=g, pi=pi, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
  # try to reassemble the full conditional out of all parts: prior, p(S|pars), and pi-comp
  p0.c <- dgamma(x=Smin.t, shape=am, rate=bm, log=TRUE)
  if (is.null(bp.t)) {
    p1.c <- sum(dpareto(S, x_min=Smin.t, k=theta.t, log=TRUE)) # likelihood
#     p1.curr <- p0.c + p1.c
    p1.curr <- dgamma(x=Smin.t, shape=am + n*theta.t[1], rate=bm, log=TRUE)
  } else {
    p1.c <- sum(dbrokenpareto(S, x_min=Smin.t, k=theta.t, bp=bp.t, log=TRUE, verbose=verbose2)) # likelihood
    p1.curr <- p0.c + p1.c
  }
  
  p2.curr <- (N.t-n)*log(1-pi.value.curr)
  if(N.t==n){
    p2.curr <- 0
  }
  p.curr <- p1.curr + p2.curr
  
  # proposal state: Smin ~ log-Normal
  eta.prop <- rnorm(n=1, mean=eta.t, sd=v.sm)    #v.sm = tuning parameter (vector) to accept 20-60% of proposals
  Smin.prop <- exp(eta.prop)
  
  if (0 < Smin.prop && all(Smin.prop < bp.t) && Smin.prop < upperBound) {
    pi.value.prop <- pi.theta.get(theta=theta.t, Smin=Smin.prop, bp=bp.t, gamma=gamma, 
                                  E=E, nsamples=nsamples, g=g, pi=pi, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
    
    # try to reassemble the full conditional out of all parts: prior, p(S|pars), and pi-comp
    p0.p <- dgamma(x=Smin.prop, shape=am, rate=bm, log=TRUE)
    if (is.null(bp.t)) {
      p1.p <- sum(dpareto(S, x_min=Smin.prop, k=theta.t, log=TRUE)) # likelihood
#       p1.prop <- p0.p + p1.p
      p1.prop <- dgamma(x=Smin.prop, shape=am + n*theta.t[1], rate=bm, log=TRUE)
    } else {
      p1.p <- sum(dbrokenpareto(S, x_min=Smin.prop, k=theta.t, bp=bp.t, log=TRUE, verbose=verbose2)) # likelihood
      p1.prop <- p0.p + p1.p
    }
    p2.prop <- (N.t-n)*log(1-pi.value.prop)
    if(N.t==n){
      p2.prop <- 0
    }
    p.prop <- p1.prop + p2.prop

    # transition density: log-normal
    q.curr.to.prop <- dlnorm(x=Smin.prop, meanlog=eta.t,    sdlog=v.sm, log=TRUE)
    q.prop.to.curr <- dlnorm(x=Smin.t,    meanlog=eta.prop, sdlog=v.sm, log=TRUE)
    
  } else {
    p.prop <- q.prop.to.curr <- -Inf
    q.curr.to.prop <- 0
    p1.prop <- p2.prop <- jacobian.prop <- NA
  }
  
  if(verbose){
    cat("--- Inside updating step for S.obs.t: --- \n")
    cat("Smin.t:     #before update\n"); print(Smin.t)
    cat("Smin.prop:\n"); print(Smin.prop)
  }
        
  # Compute ratio of densities for MH
  proposal.ratio <- q.prop.to.curr - q.curr.to.prop
  dens.ratio <- p.prop - p.curr
  log.alpha <- proposal.ratio + dens.ratio           #compare proposal to current
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
