"met.sobs" <- function(S.obs.t, theta.t, Y.obs.tot, n, Smin.t, bp.t, v.so, E, gamma, g, sigma,
                       proposal="truncnorm", verbose=FALSE ){
  
  ####################################################################################################
  # met.sobs   Metropolis step for single draw from posterior for S.obs
  #
  # Input: S.obs.t = previous iteration of S.obs
  #        theta   = theta of current iteration
  #        n       = number of observed sources
  #        Y.obs.tot = total (bkg+src) observed photon counts
  #        v.so     = vector of tuning parameters of SD to accept 20-60% of proposals
  #        E       = exposure for observed sources
  #        Smin.t  = minimum flux the sources can be detected to based on E.obs and g
  #        bp.t    = breakpint flux for the sources 
  #        gamma   = constant of transformation, energy per photon
  #        g       = function, probability of observing a source
  #        proposal = "truncnorm"/"normal" proposal distribution
  #        verbose = (T/F) display progress of program
  #
  # Output: S.obs.t = posterior samples of parameter: flux
  #         idx     = indices of accepted proposals
  ####################################################################################################

  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  
  # current state
  lambda.curr <- S.obs.t * E/gamma
  g.curr <- g(lambda=lambda.curr)

  p1.curr <- dbrokenpareto(S.obs.t, x_min=Smin.t, k=theta.t, bp=bp.t, log=TRUE, verbose=verbose2)
  p2.curr <- dpois(x=Y.obs.tot, lambda=lambda.curr, log=TRUE) #density of Y.obs.tot, with k=0. (Pois+Binom term)
  p.curr <- p1.curr + p2.curr + log(g.curr)

  # proposal state
  if (proposal=="truncnorm") {
    # Propose truncated-normal parameters for Sobs
    S.obs.prop <- rtnorm(n=n, mean=S.obs.t, sd=v.so, lower=Smin.t)
    bad.prop <- NULL
    # Compute asymmetric proposal density
    q.prop.to.curr <- dtnorm(x=S.obs.t,    mean=S.obs.prop, sd=v.so, lower=Smin.t, log=TRUE)
    q.curr.to.prop <- dtnorm(x=S.obs.prop, mean=S.obs.t,    sd=v.so, lower=Smin.t, log=TRUE)

  } else if (proposal=="normal") {
    # Propose from normal distribution
    S.obs.prop <- rnorm(n, mean=S.obs.t, sd=v.so)    #v.so = tuning parameter (vector) to accept 20-60% of proposals  
    bad.prop <- (S.obs.prop <= Smin.t)    # identify impossible proposal values
    # Compute asymmetric proposal density
    q.prop.to.curr <- dnorm(x=S.obs.t,    mean=S.obs.prop, sd=v.so, log=TRUE)
    q.curr.to.prop <- dnorm(x=S.obs.prop, mean=S.obs.t,    sd=v.so, log=TRUE)
  }

  lambda.prop <- S.obs.prop * E/gamma           #vector of proposed lambda for observed cources, not saved
  lambda.prop[bad.prop] <- NA             #set the negative lambda proposals to NA
  g.prop <- g(lambda=lambda.prop)   #vector of proposed prob. of observing sources, not saved  
  
  p1.prop <- dbrokenpareto(S.obs.prop, x_min=Smin.t, k=theta.t, bp=bp.t, log=TRUE, verbose=verbose2)
  p2.prop <- dpois(x=Y.obs.tot, lambda=lambda.prop, log=TRUE) #density of Y.obs.tot, with k=0.
  p.prop <- p1.prop + p2.prop + log(g.prop)
  p.prop[bad.prop] <- -Inf
  
  if(verbose){
    cat("--- Inside updating step for S.obs.t: --- \n")
    cat("S.obs.t:     #before update\n"); print(S.obs.t)
    cat("S.obs.prop:\n"); print(S.obs.prop)
    if(verbose2){
      cat("Debugging of samples of observed flux in met.sobs():\n")
      cat("Moment estimate of S.obs:\n")
      print(Y.obs.tot*gamma/E)
      cat("g.curr:\n"); print(g.curr)
      cat("g.prop:\n"); print(g.prop)      
      cat("theta:\n"); print(theta.t)
      cat("Smin.t:\n"); print(Smin.t)
      cat("bp.t:\n"); print(bp.t)
      cat("Y.obs.tot:\n"); print(Y.obs.tot)
    }
  }
  
  # Compute ratio of densities for MH
  prob.ratio <- p.prop - p.curr
  proposal.ratio <- q.prop.to.curr - q.curr.to.prop

  log.alpha <- prob.ratio + proposal.ratio            #compare proposal to current
  log.u <- log(runif(n))                #random vector U(0,1)
  idx <- (log.u < log.alpha)            #accepted proposal indicator 
  idx[is.na(idx)] <- FALSE
  S.obs.t[idx] <- S.obs.prop[idx]       #update accepted proposals, iteration t=iter+1


  if(verbose){
    # Print above results (error checking)
    cat("lambda.curr:\n"); print(lambda.curr)
    cat("lambda.prop:\n"); print(lambda.prop)
    cat("log of g.curr:\n"); print(log(g.curr))
    cat("log of g.prop:\n"); print(log(g.prop))
    cat("S.obs.t:     #after update\n"); print(S.obs.t)          #after update
    cat("S.obs.prop:\n"); print(S.obs.prop)
    cat("dpareto.curr:\n"); print(p1.curr)
    cat("dpareto.prop:\n"); print(p1.prop)
    cat("dpois.curr:\n"); print(p2.curr)
    cat("dpois.prop:\n"); print(p2.prop)
    cat("p.curr:\n"); print(p.curr)
    cat("p.prop:\n"); print(p.prop)
    cat("prob.ratio:\n"); print(prob.ratio)
    cat("proposal:\n"); print(proposal)
    cat("proposal.ratio:\n"); print(proposal.ratio)
    cat("log.alpha:\n"); print(log.alpha)
    cat("log.u:\n"); print(log.u)
    cat("idx:\n"); print(idx)
  }
  
  return(list('S.obs.t'=S.obs.t,'idx'=idx))  
}

met.sobs <- cmpfun(met.sobs)
