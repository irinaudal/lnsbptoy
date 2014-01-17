"met.sobs" <- function(S.obs.t, theta.t, n, Y.obs.tot, Y.obs.src.t, v.S, E.obs, L.obs, bg.obs, k.obs,
                       Smin, gamma, g.type, g, bp=NULL, p.t=NULL, use.bp=FALSE, use.mix=FALSE, 
                       proposal="truncnorm", verbose=FALSE ){
  
  ####################################################################################################
  # met.sobs   Metropolis step for single draw from posterior for S.obs
  #
  # Input: S.obs.t = previous iteration of S.obs
  #        theta.t = theta of current iteration
  #        n       = number of observed sources
  #        Y.obs.tot   = total (bkg+src) observed photon counts
  #        Y.obs.src.t = observed source photon counts of current iteration
  #        p.t     = previous vector of mixture component probabilities
  #        v.S     = vector of tuning parameters of SD to accept 20-60% of proposals
  #        E.obs   = exposure for observed sources
  #        L.obs   = location for observed sources
  #        bg.obs  = background counts/pixel
  #        k.obs   = intensity/expected background counts for observed sources
  #        Smin    = minimum flux the sources can be detected to based on E.obs and g. 
  #        gamma   = constant of transformation, energy per photon
  #        g       = function, probability of observing a source
  #        g.type  = type of g-function {"step","smooth","table"}
  #        use.bp  = (T/F) break-point pareto version
  #        use.mix = (T/F) mixture pareto version
  #        proposal = "truncnorm"/"normal" proposal distribution
  #        verbose = (T/F) display progress of program
  #
  # Output: S.obs.t = posterior samples of parameter: flux
  #         idx     = indices of accepted proposals
  ####################################################################################################
  
  if (verbose) {
    cat("S.obs.t:     #before update\n"); print(S.obs.t)          #before update
  }   
  
  # Metropolis algorithm for sampling from p(S.obs.t| .) for whole vector of S.obs, i=1,...,n 
  lambda.curr <- S.obs.t*E.obs/gamma         #vector of lambda for observed cources, not saved  
  g.curr <- g(lambda=lambda.curr, bg=bg.obs, E=E.obs, L=L.obs, g.type=g.type)    #vector of prob. of observing sources, not saved  
  
  if (proposal=="truncnorm") {
    # Propose truncated-normal parameters for Sobs
    S.obs.prop <- rtnorm(n=n, mean=S.obs.t, sd=v.S, lower=Smin)
  } else if (proposal=="normal") {
    # Propose from normal distribution
    S.obs.prop <- rnorm(n, mean=S.obs.t, sd=v.S)    #v.S = tuning parameter (vector) to accept 20-60% of proposals   
  }
  lambda.prop <- S.obs.prop*E.obs/gamma           #vector of proposed lambda for observed cources, not saved
  lambda.prop[lambda.prop <= 0] <- NA             #set the negative lambda proposals to zero
  
  g.prop <- g(lambda=lambda.prop, bg=bg.obs, E=E.obs, L=L.obs, g.type=g.type)   #vector of proposed prob. of observing sources, not saved  
  
  if(verbose){
    cat("Debugging met.sobs():\n")
    cat("g.curr:\n")
    print(g.curr)
    cat("g.prop:\n")
    print(g.prop)      
    cat("S.obs.t:\n")
    print(S.obs.t)
    cat("S.obs.prop:\n")
    print(S.obs.prop)
    cat("theta.t:\n")
    print(theta.t)
    cat("Smin:\n")
    print(Smin)
    cat("bg.obs:\n")
    print(bg.obs)
    cat("Moment estimate of S.obs:\n")
    print(Y.obs.tot*gamma/E.obs)
  }
  
  if (use.bp){
    
    pareto.curr <- dbrokenpareto(S.obs.t, x_min=Smin, k=theta.t, bp=bp, log=TRUE)        + log(g.curr)
    pareto.prop <- dbrokenpareto(S.obs.prop, x_min=Smin, k=theta.t, bp=bp, log=TRUE)     + log(g.prop)
        
  } else if (use.mix){
    
    pareto.curr <- dmixpareto(S.obs.t, p=p.t, x_min=Smin, k=theta.t, log=TRUE)        + log(g.curr)
    pareto.prop <- dmixpareto(S.obs.prop, p=p.t, x_min=Smin, k=theta.t, log=TRUE)     + log(g.prop)    
    
  } else { #No mixtures case
    
    pareto.curr <- dpareto(S.obs.t, x_min=Smin, k=theta.t, log=TRUE)        + log(g.curr)
    pareto.prop <- dpareto(S.obs.prop, x_min=Smin, k=theta.t, log=TRUE)     + log(g.prop)    
    
  }
  
  if (verbose){
    cat("pareto.curr:\n")
    print(pareto.curr)      
    cat("pareto.prop:\n")
    print(pareto.prop)
  }
  
  p.curr <- pareto.curr + 
    dbinom(Y.obs.src.t, size=Y.obs.tot, prob=lambda.curr/(k.obs+lambda.curr), log=TRUE) +
    dpois(Y.obs.tot, lambda=k.obs+lambda.curr, log=TRUE)    #vector of log(target-distr) evaluated at S.obs.t 
  p.prop <- pareto.prop +
    dbinom(Y.obs.src.t, size=Y.obs.tot, prob=lambda.prop/(k.obs+lambda.prop), log=TRUE) +
    dpois(Y.obs.tot, lambda=k.obs+lambda.prop, log=TRUE)    #vector of log(target-distr) evaluated at S.obs.prop   
  
  # Compute asymmetric proposal ratio
  if (proposal=="truncnorm") {
    # Propose truncated-normal parameters for Sobs
    q.prop.to.curr <- dtnorm(x=S.obs.t,   mean=S.obs.prop,sd=v.S,lower=Smin,log=TRUE)
    q.curr.to.prop <- dtnorm(x=S.obs.prop,mean=S.obs.t,   sd=v.S,lower=Smin,log=TRUE)
  } else if (proposal=="normal") {
    # Propose from normal distribution
    q.prop.to.curr <- dnorm(x=S.obs.t,   mean=S.obs.prop,sd=v.S,log=TRUE)
    q.curr.to.prop <- dnorm(x=S.obs.prop,mean=S.obs.t,   sd=v.S,log=TRUE)
  }
    
  # Compute ratio of densities for MH
  prob.ratio <- p.prop - p.curr
  proposal.ratio <- q.prop.to.curr - q.curr.to.prop
  
  log.alpha <- prob.ratio + proposal.ratio            #compare proposal to current
  log.u <- log(runif(n))                #random vector U(0,1)
  idx <- (log.u < log.alpha)            #accepted proposal indicator 
  idx[is.na(idx)] <- FALSE
  S.obs.t[idx] <- S.obs.prop[idx]       #update accepted proposals, iteration t=iter+1
  
  if (verbose==TRUE) {
    # Print above results (error checking)
    cat("Y.obs.tot:\n"); print(Y.obs.tot)
    cat("Y.obs.src.t:\n"); print(Y.obs.src.t)
    cat("k.obs:\n"); print(k.obs)
    cat("bg.obs:\n"); print(bg.obs)
    cat("theta.t:\n"); print(theta.t)
    cat("lambda.curr:\n"); print(lambda.curr)
    cat("lambda.prop:\n"); print(lambda.prop)
    cat("log of g.curr:\n"); print(log(g.curr))
    cat("log of g.prop:\n"); print(log(g.prop))
    cat("S.obs.t:     #after update\n"); print(S.obs.t)          #after update
    cat("S.obs.prop:\n"); print(S.obs.prop)  
    if (!use.bp & !use.mix){
      cat("dpatero.curr: (-Inf means that current.Smin > S.t)\n"); print(dpareto(S.obs.t, x_min=Smin, k=theta.t, log=TRUE))
      cat("dpatero.prop:\n"); print(dpareto(S.obs.prop, x_min=Smin, k=theta.t, log=TRUE))
    }
    cat("dbinom.curr:\n"); print(dbinom(Y.obs.src.t, size=Y.obs.tot, prob=lambda.curr/(k.obs+lambda.curr), log=TRUE))
    cat("dbinom.prop:\n"); print(dbinom(Y.obs.src.t, size=Y.obs.tot, prob=lambda.prop/(k.obs+lambda.prop), log=TRUE))  
    cat("dpois.curr:\n"); print(dpois(Y.obs.tot, lambda=k.obs+lambda.curr, log=TRUE))
    cat("dpois.prop:\n"); print(dpois(Y.obs.tot, lambda=k.obs+lambda.prop, log=TRUE)) 
    cat("p.curr:\n"); print(p.curr)
    cat("p.prop:\n"); print(p.prop)
    cat("prob.ratio:\n"); print(prob.ratio)
    cat("proposal.ratio:\n"); print(proposal.ratio)
    cat("log.alpha:\n"); print(log.alpha)
    cat("log.u:\n"); print(log.u)
    cat("idx:\n"); print(idx)
  }
  
  return(list('S.obs.t'=S.obs.t,'idx'=idx))  
}

met.sobs <- cmpfun(met.sobs)
