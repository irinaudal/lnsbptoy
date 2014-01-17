"met.theta" <- function(theta.t, S.obs.t, S.mis.t, N.t, n, v.theta, 
                        a, b, Smin, gamma, pble, length.S, g.type, g, pi, nsamples,
                        bp=NULL, I.idx.t=NULL, p.t=NULL, prob.MH = 0, fixed.theta=NULL,
                        fixed.S.pi=NULL, debug.pi=FALSE, use.bp=FALSE, use.mix=FALSE, 
                        algorithm.type=2, #1=uses Scom w/o binomial term,  2=uses Sobs
                        verbose=FALSE ){ 
  
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
  #        v.theta = tuning parameter of SD to accept 20-60% of proposals
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
  
  pi.value.curr <- pi.theta.get(pi=pi, theta=theta.t, p.t=p.t, bp=bp, Smin=Smin, gamma=gamma, 
                                pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi,
                                nsamples=nsamples, g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources
  
  if (verbose) {
    cat(ppaste("pi.value.curr = ",pi.value.curr,"\n"))          #before update
  }   
  
  m <- length(theta.t)
  
  use.MH <- (runif(1) < prob.MH) # rbinom(n=1,size=1,prob=prob.MH)  
  
  # Prepare gamma parameters for theta for portion of theta posterior
  if (use.bp){
    ## Broken power-law version:
    if (algorithm.type==1) {
      post.pars <- brokenpareto.posterior(a=a,b=b,S.obs=S.obs.t,S.mis=S.mis.t,Smin=Smin,bp=bp,verbose=verbose)
      a.post <- post.pars$a
      b.post <- post.pars$b
    } else if (algorithm.type==2) {
      post.pars <- brokenpareto.posterior(a=a,b=b,S.obs=S.obs.t,S.mis=NULL,Smin=Smin,bp=bp,verbose=verbose)
      a.post <- post.pars$a
      b.post <- post.pars$b
    }
    
  } else if (use.mix) {   
    ## Mixture version:
    stop("Theta sampling for Mixture model must be checekd again in derivation - posterior in the case of integrating S.mis.\n")
    post.pars <- mixpareto.posterior(a=a,b=b,S.obs=S.obs.t,S.mis=S.mis.t,Smin=Smin,I.idx=I.idx.t,verbose=verbose)
    a.post <- post.pars$a
    b.post <- post.pars$b
    
  } else {
    ## Regular version:
    if (algorithm.type==1) {
      a.post <- a + N.t
      b.post <- b + sum(log(S.obs.t/Smin)) + sum(log(S.mis.t/Smin))   
    } else if (algorithm.type==2) {
      a.post <- a + n
      b.post <- b + sum(log(S.obs.t/Smin))   
    } 
  }
  if(verbose){
    cat("Param: a in post.theta:\n"); print(a.post)
    cat("Param: b in post.theta:\n"); print(b.post)
    cat(ppaste("Smin = ",Smin,"\n"))
    if(!use.mix & !use.bp) {
      cat("sum(log(S.obs.t/Smin))\n"); print(sum(log(S.obs.t/Smin)))
      cat("sum(log(S.mis.t/Smin))\n"); print(sum(log(S.mis.t/Smin)))
      cat("log(S.obs.t)\n"); print(log(S.obs.t))
    }
  }
  
  if (!use.MH){ # Rejection Sampling for getting a draw of theta
    
    if (verbose){
      stime <- proc.time()
    }
    
    got.draw.idx <- FALSE
    max.tries <- 1000
    num.tries <- 0
    
    while (!got.draw.idx && num.tries<max.tries){
      
      theta.prop <- rgamma(n=m, shape=a.post, rate=b.post)
      
      # Overwrite/take care of fixed.theta
      fixed.idx <- (fixed.theta != FALSE)            # index of fixed values
      theta.prop[fixed.idx] <- fixed.theta[fixed.idx]  # overwrite with fixed values
      
      pi.value.prop <- pi.theta.get(pi=pi, theta=theta.prop, p.t=p.t, bp=bp, Smin=Smin, gamma=gamma, 
                                    pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi,
                                    nsamples=nsamples, g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources
      
      pi.value.prop <- ifelse(any(pi.value.prop>-Inf), pi.value.prop, NA) #sanity check
      
      if (algorithm.type==1) {
        log.r.alpha <- dbinom(n, size=N.t, prob=pi.value.prop, log=TRUE)
      } else if (algorithm.type==2) {
        log.r.alpha <- (N.t-n)*log(1-pi.value.prop)
        if(N.t==n){
          log.r.alpha <- 0
        }
      }
      if (log(runif(1)) < log.r.alpha){
        theta.t <- theta.prop
        got.draw.idx <- TRUE #returned index of theta
      }
      num.tries <- num.tries+1
      
      if (verbose>2){ #TRUE)
        cat(ppaste("Debugging rejection sampling attempt ",num.tries,":\n"))
        cat(ppaste("posterior a = \n")); print(a.post)
        cat(ppaste("posterior b = \n")); print(b.post)
        cat(ppaste("posterior 2.5th percentiles  = ",qgamma(p=0.025,shape=a.post,rate=b.post),"\n")) 
        cat(ppaste("posterior 50th percentiles   = ",qgamma(p=0.500,shape=a.post,rate=b.post),"\n")) 
        cat(ppaste("posterior 97.5th percentiles = ",qgamma(p=0.975,shape=a.post,rate=b.post),"\n"))
        cat(ppaste("proposed theta = \n")); print(theta.prop)
        cat(ppaste("selected theta = \n")); print(theta.t)
        cat(ppaste("log.r.alpha    = ",log.r.alpha,"\n"))
        cat(ppaste("log.u          = ",log.u,"\n"))
        cat(ppaste("idx            = ",got.draw.idx,"\n"))
        cat(ppaste("pi(theta)= ",pi.value.prop,"\n"))
        cat(ppaste("N.t      = ",N.t,"\n"))
        cat(ppaste("E[Binom] = ",N.t*pi.value.prop,"\n"))
        cat(ppaste("n        = ",n,"\n"))
        cat(ppaste("dbinom   = ",exp(log.r.alpha),"\n"))
      } 
    } # END while-loop
    
    if (verbose>0){   
      etime <- proc.time()
      cat(ppaste("--- Rejection sampling for theta took ",num.tries," attempts (",(etime-stime)["elapsed"],")...\n"))
      cat(ppaste("selected theta = \n")); print(theta.t)
      if (verbose>1) {
        cat(ppaste("Last step of rejection sampling: attempt ",num.tries,":\n"))
        cat(ppaste("posterior a = \n")); print(a.post)
        cat(ppaste("posterior b = \n")); print(b.post)
        cat(ppaste("posterior 2.5th percentiles  = ",qgamma(p=0.025,shape=a.post,rate=b.post),"\n")) 
        cat(ppaste("posterior 50th percentiles   = ",qgamma(p=0.500,shape=a.post,rate=b.post),"\n")) 
        cat(ppaste("posterior 97.5th percentiles = ",qgamma(p=0.975,shape=a.post,rate=b.post),"\n"))
        cat(ppaste("pi(theta)= ",pi.value.prop,"\n"))
        cat(ppaste("N.t      = ",N.t,"\n"))
        cat(ppaste("E[Binom] = ",N.t*pi.value.prop,"\n"))
        cat(ppaste("n        = ",n,"\n"))
        cat(ppaste("dbinom   = ",exp(log.r.alpha),"\n"))
        cat("...........................................\n")
      }
    }
    
  } ## END rejection sampling piece...
  
  if (use.MH){
    
    #Metropolis algorithm for sampling from p(theta.t| .) for whole vector of theta, i=1,...,m  
    
    if(use.mix){
      stop("Mixture model should be checked first.")        
    } else {
      # Regular OR bp model
      p.curr.1 <- sum(dgamma(x=theta.t, shape=a.post, rate=b.post, log=TRUE)) 
      
      if (algorithm.type==1) {
        p.curr.2 <- dbinom(n, size=N.t, prob=pi.value.curr, log=TRUE)
      } else if (algorithm.type==2) {
        p.curr.2 <- (N.t-n)*log(1-pi.value.curr)
        if(N.t==n){
          p.curr.2 <- 0
        }
      }
      p.curr <- p.curr.1 + p.curr.2
    }
    
    theta.prop <- rnorm(n=m, mean=theta.t, sd=v.theta)    #v.theta = tuning parameter to accept 20-60% of proposals 
    
    # Overwrite/take care of fixed.theta
    fixed.idx <- (fixed.theta != FALSE)            # index of fixed values
    theta.prop[fixed.idx] <- fixed.theta[fixed.idx]  # overwrite with fixed values 
    
    if (all(theta.prop>0)){
      
      pi.value.prop <- pi.theta.get(pi=pi, theta=theta.prop, p.t=p.t, bp=bp, Smin=Smin, gamma=gamma, 
                                    pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi,
                                    nsamples=nsamples, g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources
      
      pi.value.prop <- ifelse(any(pi.value.prop>-Inf), pi.value.prop, NA) #sanity check
      
      p.prop.1 <- sum(dgamma(x=theta.prop, shape=a.post, rate=b.post, log=TRUE)) 
      
      if (algorithm.type==1) {
        p.prop.2 <- dbinom(n, size=N.t, prob=pi.value.prop, log=TRUE)
      } else if (algorithm.type==2) {
        p.prop.2 <- (N.t-n)*log(1-pi.value.prop)
        if(N.t==n){
          p.prop.2 <- 0
        }
      }
      p.prop <- p.prop.1 + p.prop.2
      
    } else {
      # Avoid bunch of annoying dpareto warnings...
      pi.value.prop <- 0
      p.prop <- -Inf
    }
    
    log.alpha <- p.prop-p.curr            #compare proposal to current
    log.u <- log(runif(1))                #random scalar U(0,1)
    got.draw.idx <- (log.u < log.alpha)            #accepted proposal indicator 
    got.draw.idx[is.na(got.draw.idx)] <- FALSE
    theta.t[got.draw.idx] <- theta.prop[got.draw.idx]       #update accepted proposals, iteration t=iter+1
    
    if (verbose>2){ 
      cat(ppaste("Debugging MH sampling for theta:\n"))
      cat(ppaste("posterior a = \n")); print(a.post)
      cat(ppaste("posterior b = \n")); print(b.post)
      cat(ppaste("posterior 2.5th percentiles  = ",qgamma(p=0.025,shape=a.post,rate=b.post),"\n")) 
      cat(ppaste("posterior 50th percentiles   = ",qgamma(p=0.500,shape=a.post,rate=b.post),"\n")) 
      cat(ppaste("posterior 97.5th percentiles = ",qgamma(p=0.975,shape=a.post,rate=b.post),"\n"))
      cat(ppaste("proposed theta = \n")); print(theta.prop)
      cat(ppaste("selected theta = \n")); print(theta.t)
      cat(ppaste("pi(theta)= ",pi.value.prop,"\n"))
      cat(ppaste("N.t      = ",N.t,"\n"))
      cat(ppaste("E[Binom] = ",N.t*pi.value.prop,"\n"))
      cat(ppaste("n        = ",n,"\n"))
      cat(ppaste("dbinom   = ",exp(log.r.alpha),"\n"))
    } 
    
    if (verbose>0 && use.MH) {
      #print above results (error checking)
      cat(ppaste("--- MH sampling for theta: results:\n"))
      if (verbose>1) {
        cat(ppaste("m = ",m,"\n")); print(m)
        cat(ppaste("v.theta = ")); print(v.theta)
        cat(ppaste("log of pi.value.curr     = ",log(pi.value.curr),"\n"))
        cat(ppaste("log of pi.value.prop     = ",log(pi.value.prop),"\n"))
        cat(ppaste("log of (1-pi.curr)^(N-n) = ",(N.t-n)*log(1-pi.value.curr),"\n"))
        cat(ppaste("log of (1-pi.prop)^(N-n) = ",ifelse(all(theta.prop>0), (N.t-n)*log(1-pi.value.prop),-Inf),"\n"))
        cat(ppaste("sum of log of dgamma.curr = ",sum(dgamma(x=theta.t, shape=a.post, rate=b.post, log=TRUE)),"\n"))
        cat(ppaste("sum of log of dgamma.prop = ",sum(dgamma(x=theta.prop, shape=a.post, rate=b.post, log=TRUE)),"\n"))
      }
      cat(ppaste("p.curr    = ",p.curr,"\n"))
      cat(ppaste("p.prop    = ",p.prop,"\n"))
      cat(ppaste("log.alpha = ",log.alpha,"\n"))
      cat(ppaste("log.u     = ",log.u,"\n"))
      cat(ppaste("idx       = ",got.draw.idx,"\n"))
      cat(ppaste("proposed theta = \n")); print(theta.prop)
      cat(ppaste("selected theta = \n")); print(theta.t)
    }     
    
  } ## END if (use.MH){...}
  
  return(list('theta.t'=theta.t,"idx"=got.draw.idx))  
}

met.theta <- cmpfun(met.theta)
