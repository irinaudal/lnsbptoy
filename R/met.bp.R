"met.bp" <- function(bp.t, Smin.t, S.obs.t, theta.t, N.t, n, v.bp,
                     C, mu, gamma, E, g, nsamples, pi, sigma,
                     fixed.bp=NULL, 
                     TOL=1e-10,
                     verbose=FALSE ){
  
  ####################################################################################################
  # met.bp   Metropolis-Hasting step for single draw from posterior for bp
  #
  # Input: bp.t    = previous iteration of break points bp=[bp1,bp2,bp3,...,bp_m]
  #        Smin.t  = minimum flux the sources can be detected to according to E.obs and g.
  #        S.obs.t = current vector of observed sources of flux
  #        S.mis.t = current vector of missing sources of flux
  #        theta.t = current iteration of theta
  #        N.t     = total number of sources, current iteration
  #        n.t     = number of observed sources
  #        v.bp    = tuning parameter of SD to accept 20-60% of proposals of bp.t
  #        z.e     = hit-run Metropolis specified parameter of generating non-random direction
  #        mu      = mean parameter in normal proposal for break-points, must be of correct dimensions
  #        C       = varianca-covariance matrix of parameters in normal proposal for break-points
  #        gamma   = constant of transformation, energy per photon
  #        pble    = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        length.S = length of grid of S, for integration of g(theta)
  #        g.type   = type of g-function {"step","smooth","table"}
  #        g        = function, probability of observing source
  #        pi       = R.object: a list("pi","theta") with detection probalities: pi(theta)
  #        fixed.bp = Value/FALSE: allow for user specified bp
  #        fixed.S.pi = (override) fix the value of S (used mainly for debugging)
  #        debug.pi   = (T/F) (used mainly for debugging)  
  #        met.alg  = "MH"/"MTM"/("HitRun") Metropolis method variations, used in Smin, bp sampling
  #        kk       = number of tries to propose new draws in Multiple-Try-Metropolis
  #        use.bp   = (T/F) break-point pareto version
  #        use.mix  = (T/F) mixture pareto version
  #        verbose  = (T/F) display progress of program
  #
  # Output: bp.t    = posterior samples of parameter: bp
  #         idx     = indices of accepted proposals
  ####################################################################################################
  
  ### NOTE: bp MUST be sampled separately from Smin
  
  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  
  S <- S.obs.t
  m.bp <- length(theta.t)-1
  
  # current state
  tau <- c(Smin.t, bp.t)
  eta.t <- log(diff(tau))
  
  pi.value.curr <- pi.theta.get(theta=theta.t, Smin=Smin.t, bp=bp.t, gamma=gamma, 
                                E=E, nsamples=nsamples, g=g, pi=pi, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
  # pi.value.curr <- ifelse(any(pi.value.curr>-Inf), pi.value, NA) #sanity check
  
  p1.curr <- sum(dnorm(x=eta.t, mean=mu, sd=C, log=TRUE)) # prior, no Jacobian
  p2.curr <- sum(dbrokenpareto(S, x_min=Smin.t, k=theta.t, bp=bp.t, log=TRUE, verbose=verbose2)) # likelihood
  p3.curr <- (N.t-n)*log(1-pi.value.curr)
  if(N.t==n){
    p3.curr <- 0
  }
  p.curr <- p1.curr + p2.curr + p3.curr
  
  
  # proposal state
  eta.prop <- rnorm(n=m.bp, mean=eta.t, sd=v.bp)    #v.bp = tuning parameter (vector) to accept 20-60% of proposals  
  bp.prop <- cumsum(exp(eta.prop)) + Smin.t
  tau <- c(Smin.t,bp.prop)
  eta.after <- log(diff(tau))
  
  # Overwrite/take care of fixed.bp
  fixed.idx  <- (fixed.bp != FALSE)
  # Check for numerical accuracy
  if (abs(eta.prop-eta.after)>TOL) {
    warning("Proposal of z(bp,Smin) is too small, numerical error incured with transformation.")
    bp.prop <- -Inf
  } else if (any(fixed.idx)) {
    bp.prop[fixed.idx] <- fixed.bp[fixed.idx]  # overwrite with fixed values 
    tau <- c(Smin.t,bp.prop)
    eta.prop <- log(diff(tau))
  }
  
  pi.value.prop <- pi.theta.get(theta=theta.t, Smin=Smin.t, bp=bp.t, gamma=gamma, 
                                E=E, nsamples=nsamples, g=g, pi=pi, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
  # pi.value.prop <- ifelse(any(pi.value.prop>-Inf), pi.value, NA) #sanity check
  
  if (any(Smin.t > bp.prop)) {
    p.prop <- -Inf
    p1.prop <- p2.prop <- p3.prop <- NA
  } else {
    p1.prop <- sum(dnorm(x=eta.prop, mean=mu, sd=C, log=TRUE)) # prior, no Jacobian
    p2.prop <- sum(dbrokenpareto(S, x_min=Smin.t, k=theta.t, bp=bp.prop, log=TRUE, verbose=verbose2)) # likelihood
    p3.prop <- (N.t-n)*log(1-pi.value.prop)
    if(N.t==n){
      p3.prop <- 0
    }
    p.prop <- p1.prop + p2.prop + p3.prop
  }
  
  if(verbose){
    cat("--- Inside updating step for eta/bp.t: --- \n")
    cat("eta.t:     #before update\n"); print(eta.t)
    cat("eta.prop:\n"); print(eta.prop)
    cat("bp.t:     #before update\n"); print(bp.t)
    cat("bp.prop:\n"); print(bp.prop)
    cat(paste("v.bp :eta/bp = ")); print(v.bp)
    if(verbose2){
      cat("Debugging of samples of observed flux in met.bp():\n")
      cat("theta:\n"); print(theta.t)
      cat("Smin.t:\n"); print(Smin.t)
      cat("N.t:\n"); print(N.t)
    }
  }
  
  # Compute ratio of densities for MH
  log.alpha <- p.prop - p.curr          #compare proposal to current
  log.u <- log(runif(1))                #random vector U(0,1)
  got.draw.idx <- (log.u < log.alpha)            #accepted proposal indicator 
  got.draw.idx[is.na(got.draw.idx)] <- FALSE
  bp.t[got.draw.idx] <- bp.prop[got.draw.idx]       #update accepted proposals, iteration t=iter+1
  
  
  if(verbose){
    # Print above results (error checking)
    cat("log of pi.curr:\n"); print(log(pi.value.curr))
    cat("log of pi.prop:\n"); print(log(pi.value.prop))
    cat("eta.t:     #after update\n"); print(eta.t)          #after update
    cat("eta.prop:\n"); print(eta.prop)
    cat("bp.t:     #after update\n"); print(bp.t)          #after update
    cat("bp.prop:\n"); print(bp.prop)
    cat("dpareto.curr:\n"); print(p2.curr)
    cat("dpareto.prop:\n"); print(p2.prop)
    cat("normal.curr:\n"); print(p1.curr)
    cat("normal.prop:\n"); print(p1.prop)
    cat("p.curr:\n"); print(p.curr)
    cat("p.prop:\n"); print(p.prop)
    cat("log.alpha:\n"); print(log.alpha)
    cat("log.u:\n"); print(log.u)
    cat("got.draw.idx:\n"); print(got.draw.idx)
  }
  
 
  return(list('bp.t'=bp.t,"idx"=got.draw.idx))
}

met.bp <- cmpfun(met.bp)