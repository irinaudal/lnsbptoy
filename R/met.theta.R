"met.theta" <- function(theta.t, Smin.t, bp.t, S.obs.t, N.t, n, v.th, 
                        a, b, gamma, E, g, nsamples, sigma,
                        fixed.theta=NULL, verbose=FALSE ){ 
  
  ####################################################################################################
  # met.theta   Metropolis and Rejection Sampling step for single draw from posterior for theta
  #
  # Input: theta.t = previous iteration of theta, lns slope
  #        S.obs.t = S.obs of current iteration
  #        S.mis.t = S.mis of current iteration
  #        N.t     = total number of sources, current iteration
  #        I.idx.t = previous vector of mixture component indices
  #        p.t     = previous vector of mixture component probabilities
  #        n       = number of observed sources
  #        v.thheta = tuning parameter of SD to accept 20-60% of proposals
  #        a       = shape hyper-parameter in gamma prior for theta
  #        b       = rate hyper-parameter in gamma prior for theta
  #        Smin    = minimum flux the sources can be detected to according to E.obs and g.
  #        gamma   = constant of transformation, energy per photon
  #        pble    = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        length.S = length of grid of S, for integration of g(theta)
  #        pi       = R.object: a list("pi","theta") with detection probalities: pi(theta)
  #        bp      = either NULL (no fixed break-points), or, a vector of length p specifying the break-points
  #        prob.MH    = probability of how often to perform Metropolis-Hastings step in sampling theta  
  #        fixed.S.pi = (override) fix the value of S (used mainly for debugging)
  #        nsamples   = number of MC samples to produce estimate of pi(theta,Smin) integral
  #        debug.pi   = (T/F) (used mainly for debugging)
  #        fixed.theta = Value/FALSE: allow for user specified theta - no mcmc
  #        g.type  = type of g-function {"step","smooth","table"}
  #        g       = function, probability of observing source
  #        use.bp   = (T/F) break-point pareto version
  #        use.mix  = (T/F) mixture pareto version
  #        algorithm.type = value:(1)/(2)  #1=uses Scom w/o binomial term,  2=uses Sobs
  #        verbose  = (T/F) display progress of program
  #
  # Output: theta.t = posterior samples of parameter: theta
  #         idx     = indices of accepted proposals
  ####################################################################################################
  
  
  if (verbose) {
    cat("--- Inside updating step for theta.t: --- \n")
    cat("theta.t:     #before update\n"); print(theta.t)          #before update
  } 
  
  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  m <- length(theta.t)
  
  # Prepare gamma parameters for theta for portion of theta posterior
  if (is.null(bp.t)) {
    a.post <- a + n
    b.post <- b + sum(log(S.obs.t/Smin))   
  } else {
    post.pars <- brokenpareto.posterior(a=a,b=b,S.obs=S.obs.t,S.mis=NULL,Smin=Smin.t,bp=bp.t,verbose=verbose)
    a.post <- post.pars$a
    b.post <- post.pars$b
  }
  
  if(verbose){
    cat("Param: a in post.theta:\n"); print(a.post)
    cat("Param: b in post.theta:\n"); print(b.post)
    cat(paste("Smin.t =",Smin.t,"\n"))
    cat(paste("bp.t =",bp.t,"\n"))
  }
  
  # current state
  pi.value.curr <- pi.theta.get(theta=theta.t, Smin=Smin.t, bp=bp.t, gamma=gamma, 
                                E=E, nsamples=nsamples, g=g, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
  p1.curr <- sum(dgamma(x=theta.t, shape=a.post, rate=b.post, log=TRUE))
  p2.curr <- (N.t-n)*log(1-pi.value.curr)
  if(N.t==n){
    p2.curr <- 0
  }
  p.curr <- p1.curr + p2.curr
  
  # proposed state
  theta.prop <- rnorm(n=m, mean=theta.t, sd=v.th)    #v.th = tuning parameter to accept 20-60% of proposals 
  # Overwrite/take care of fixed.theta
  fixed.idx <- (fixed.theta != FALSE)            # index of fixed values
  theta.prop[fixed.idx] <- fixed.theta[fixed.idx]  # overwrite with fixed values 
  
  if (all(theta.prop>0)){
    pi.value.prop <- pi.theta.get(theta=theta.prop, Smin=Smin.t, bp=bp.t, gamma=gamma, 
                                  E=E, nsamples=nsamples, g=g, sigma=sigma, verbose=verbose2)     #marginal prob. of observing sources  
    p1.prop <- sum(dgamma(x=theta.prop, shape=a.post, rate=b.post, log=TRUE))
    p2.prop <- (N.t-n)*log(1-pi.value.prop)
    if(N.t==n){
      p2.prop <- 0
    }
    p.prop <- p1.prop + p2.prop
    
  } else {
    p.prop <- -Inf
    p1.prop <- p2.prop <- NA
  }
  
  if(verbose){
    cat(paste("pi.value.curr =",pi.value.curr,"\n"))          #before update
  }
  
  
  log.alpha <- p.prop-p.curr            #compare proposal to current
  log.u <- log(runif(1))                #random scalar U(0,1)
  got.draw.idx <- (log.u < log.alpha)            #accepted proposal indicator 
  got.draw.idx[is.na(got.draw.idx)] <- FALSE
  theta.t[got.draw.idx] <- theta.prop[got.draw.idx]       #update accepted proposals, iteration t=iter+1
  
  if (verbose>2){ 
    cat(paste("Debugging MH sampling for theta:\n"))
    cat(paste("posterior a = \n")); print(a.post)
    cat(paste("posterior b = \n")); print(b.post)
    cat(paste("posterior 2.5th percentiles  =",qgamma(p=0.025,shape=a.post,rate=b.post),"\n")) 
    cat(paste("posterior 50th percentiles   =",qgamma(p=0.500,shape=a.post,rate=b.post),"\n")) 
    cat(paste("posterior 97.5th percentiles =",qgamma(p=0.975,shape=a.post,rate=b.post),"\n"))
    cat(paste("proposed theta =\n")); print(theta.prop)
    cat(paste("selected theta =\n")); print(theta.t)
    cat(paste("pi(theta)=",pi.value.prop,"\n"))
    cat(paste("N.t      =",N.t,"\n"))
    cat(paste("E[Binom] =",N.t*pi.value.prop,"\n"))
    cat(paste("n        =",n,"\n"))
    #cat(paste("dbinom   =",exp(log.r.alpha),"\n"))
  } 
  
  if (verbose>0) {
    #print above results (error checking)
    cat(paste("--- MH sampling for theta: results:\n"))
    if (verbose>1) {
      cat(paste("m = ",m,"\n")); print(m)
      cat(paste("v.th :theta = ")); print(v.th)
      cat(paste("log of pi.value.curr     =",log(pi.value.curr),"\n"))
      cat(paste("log of pi.value.prop     =",log(pi.value.prop),"\n"))
      cat(paste("log of (1-pi.curr)^(N-n) =",(N.t-n)*log(1-pi.value.curr),"\n"))
      cat(paste("log of (1-pi.prop)^(N-n) =",ifelse(all(theta.prop>0), (N.t-n)*log(1-pi.value.prop),-Inf),"\n"))
      cat(paste("sum of log of dgamma.curr =",p1.curr,"\n"))
      cat(paste("sum of log of dgamma.prop =",p1.prop,"\n"))
    }
    cat(paste("p.curr    = ",p.curr,"\n"))
    cat(paste("p.prop    = ",p.prop,"\n"))
    cat(paste("log.alpha = ",log.alpha,"\n"))
    cat(paste("log.u     = ",log.u,"\n"))
    cat(paste("idx       = ",got.draw.idx,"\n"))
    cat(paste("proposed theta = \n")); print(theta.prop)
    cat(paste("selected theta = \n")); print(theta.t)
  }     
  
  
  return(list('theta.t'=theta.t,"idx"=got.draw.idx))  
}

met.theta <- cmpfun(met.theta)
