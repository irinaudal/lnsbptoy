"met.smin" <- function(Smin.t, N.t, n.t, theta.t, S.obs.t, S.mis.t, I.idx.t, am,bm, v.sm,
                       gamma, pble, length.S, g.type, g, pi, nsamples,
                       bp.t=NULL,
                       algorithm.type=2, #1=uses Scom,  2=uses Sobs,
                       met.alg="MH", kk=5, symmetry=TRUE, # "MTM" or "MH" or "HitRun"
                       prob.MH=1, n.tries.max=5000,
                       stepsize=10^-16, n.steps=25,
                       fixed.S.pi=NULL, debug.pi=FALSE, use.bp=FALSE, use.mix=FALSE, verbose=FALSE){
  
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
  
  if (use.mix) {
    stop("Cannot sample Smin for 'mix' model - to be implemented.")
  }
  
  if (met.alg=="HitRun") {
    met.alg="MTM"  # 1-D parameter case gives MTM method
  }
  
  m <- length(theta.t)
  if (algorithm.type==1) {
    S     <- c(S.obs.t,S.mis.t)
    minSs <- min(S)
  } else if (algorithm.type==2) {
    minSs <- min(S.obs.t)
  }
  
  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  
#   use.MH <- (runif(1) < prob.MH) # rbinom(n=1,size=1,prob=prob.MH)  
#   use.HMC <- prob.MH == -1
#   
#   
#   if (use.HMC){
#     # Hamiltonian Monte Carlo
#     
#     #####################################
#     ### Define utility functions ###
#     #####################################
#     "U.fun" <- function(q, n=n.t, theta=theta.t[1], bp=bp.t, am.=am, bm.=bm, upperBound=minSs) {
#       # U(q) = -log( prob(q | else) )
#       return(ifelse(q<min(upperBound,bp), -(am. + n*theta - 1)*log(q) + bm.*q - 0, +Inf))  # TODO: omitted is the log(1-pi) term.
#     }
#     "grad.U.fun" <- function(q, n=n.t, theta=theta.t[1], bp=bp.t, am.=am, bm.=bm, upperBound=minSs) {
#       return(ifelse(q<min(upperBound,bp), -(am. + n*theta - 1)/q + bm. + 0, 0))  # TODO: omitted is the derivative of log(1-pi) term.
#     }
#     
#     tmp <- alg.hmc(q.curr=Smin.t, U=U.fun, grad.U=grad.U.fun, eps=stepsize, L=n.steps)
#     
#     got.draw.idx <- tmp$got.draw.idx
#     if (got.draw.idx) {
#       # Retrieve results
#       Smin.prop <- tmp$q.prop
#       if(verbose) {
#         cat("Selected Smin.t:\n"); print(Smin.prop)          #before update
#       }
#       Smin.t <- Smin.prop    #Smin update   
#     }    
#     if (verbose) {
#       cat("Smin.t after update:\n"); print(Smin.t)          #before update
#     } 
#     
#   } else if (!use.MH) {
#     # Rejection-Smapling algorithm
#     
#     # propose from p.1 <- dgamma(x=Smin.t, shape=N.t*theta.t[1] + am, rate=bm, log=TRUE)
#     # or
#     # propose from p.1 <- dgamma(x=Smin.t, shape=n.t*theta.t[1] + am, rate=bm, log=TRUE)
#     if (algorithm.type==1) {
#       Smin.prop <- rgamma(n=1, shape=am, rate=bm)
#     } else if (algorithm.type==2) {
#       Smin.prop <- rgamma(n=1, shape=am, rate=bm)
#     } 
#     
#     # Check boundary conditions or Try rejecting based on value of (1-pi)^N-n
#     if (all(Smin.prop < bp.t) && Smin.prop < minSs) {
#       pi.value <- pi.theta.get(pi=pi, theta=theta.t, p.t=NULL, bp=bp.t, Smin=Smin.prop, gamma=gamma, 
#                                pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi,
#                                nsamples=nsamples, g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources
#       pi.value <- ifelse(any(pi.value>-Inf), pi.value, NA) #sanity check
#       if (algorithm.type==1) {
#         p.2 <- dbinom(n.t, size=N.t, prob=pi.value,log=TRUE) +  #vector of log(target-distr) evaluated at theta.t,S.t,N.t
#                 (N.t*theta.t[1])*log(Smin.prop)
#       } else if (algorithm.type==2) {
#         p.2 <- (N.t-n.t)*log(1-pi.value)
#         if(N.t==n.t){
#           p.2 <- 0
#         }
#         p.2 <- p.2 + (n.t*theta.t[1])*log(Smin.prop)
#       }
#       log.alpha <- p.2             #compare proposal to current
#       log.u <- log(runif(1))            #random vector U(0,1)
#       idx <- (log.u < log.alpha)            #accepted proposal indicator 
#     } else {
#       idx <- FALSE
#     }
#     
#     #fail safe mechanism
#     failed.indices <- (!idx)
#     n.tries <- 0
#     
#     while (failed.indices && n.tries<n.tries.max){
#       
#       n.tries <- n.tries+1     
#       
#       # Propose new value
#       if (algorithm.type==1) {
#         #Smin.prop <- rgamma(n=1, shape=N.t*theta.t[1] + am, rate=bm)
#         Smin.prop <- rgamma(n=1, shape=am, rate=bm)
#       } else if (algorithm.type==2) {
#         #Smin.prop <- rgamma(n=1, shape=n.t*theta.t[1] + am, rate=bm)
#         Smin.prop <- rgamma(n=1, shape=am, rate=bm)
#       } 
#       
#       # Check boundary conditions or Try rejecting based on value of (1-pi)^N-n
#       if (all(Smin.prop < bp.t) && Smin.prop < minSs) {
#         pi.value <- pi.theta.get(pi=pi, theta=theta.t, p.t=NULL, bp=bp.t, Smin=Smin.prop, gamma=gamma, 
#                                  pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi,
#                                  nsamples=nsamples, g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources
#         pi.value <- ifelse(any(pi.value>-Inf), pi.value, NA) #sanity check
#         if (algorithm.type==1) {
#           p.2 <- dbinom(n.t, size=N.t, prob=pi.value,log=TRUE) +  #vector of log(target-distr) evaluated at theta.t,S.t,N.t
#                    (N.t*theta.t[1])*log(Smin.prop)
#         } else if (algorithm.type==2) {
#           p.2 <- (N.t-n.t)*log(1-pi.value)
#           if(N.t==n.t){
#             p.2 <- 0
#           }
#           p.2 <- p.2 + (n.t*theta.t[1])*log(Smin.prop)
#         }
#         log.alpha <- p.2              #compare proposal to current
#         log.u <- log(runif(1))            #random vector U(0,1)
#         idx <- (log.u < log.alpha)            #accepted proposal indicator 
#       } else {
#         idx <- FALSE
#       }
# 
#       # Check if Smin is still invalid...
#       failed.indices <- (!idx)
#       
#       # Error catch after fail safe mechanism
#       if(n.tries>=n.tries.max){
#         stop(paste("Error: Cannot draw acceptible Smin values. Check [new proposals of Smin]. 
#                    The current value of bp = ",bp.t, ", min(Si) = ",minSs,",\n 
#                    theta[1] = ",theta.t[1],", N.t=",N.t, ", n=",n.t,"\n
#                    log.alpha = ", log.alpha, "log.u = ", log.u,"\n", sep=""))
#       } else {
#         Smin.t <- Smin.prop       #update accepted proposals, iteration t=iter+1, NOTE: size will change!
#       }
#       
#     }
#     
#     got.draw.idx <- TRUE
#     if (verbose) {
#       cat("Smin.t after update:\n"); print(Smin.t)          #before update
#     }
#     
#     
#   } else {
    
  
    # Use Metropolis algorithm
  
  #####################################
  ### Define utility functions ###
  #####################################
  "r.q.prop.fun.mtm" <- function(x.curr, v=v.sm, m=1, upperBound=log(minSs), symm=symmetry) {
    if (symm) {
      # Symmetric proposal
      y.prop <- rnorm(n=m, mean=x.curr, sd=v)
    } else {  
      # Non-symmetric proposal
      y.prop <- rtnorm(n=m, mean=x.curr, sd=v, upper=upperBound)
    }
    return(y.prop)
  }

  "d.q.prop.fun.mtm" <- function(x.curr, y.prop, v=v.sm, m=1, upperBound=log(minSs), symm=symmetry) {
    if (symm) {
      # Symmetric log-density
      dq <- 0
    } else {  
      # Non-symmetric log-density
      dq <- dtnorm(x=y.prop, mean=x.curr, sd=v, upper=upperBound, log=TRUE)
    }
    return(dq)
  }
  
  "log.target.fun" <- function(x, bp=bp.t, upperBound=minSs, 
                               pi.vec=pi, theta=theta.t, N=N.t, n=n.t, am.=am, bm.=bm,...) {
    
    # Transform variables
    z <- x
    Smin <- exp(z)
    
    if (0 < Smin && all(Smin < bp) && Smin < upperBound) { # check boundary conditions
      
      pi.value <- pi.theta.get(pi=pi.vec, theta=theta, p.t=NULL, bp=bp, Smin=Smin, gamma=gamma, 
                               pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi,
                               nsamples=nsamples, g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources
      pi.value <- ifelse(any(pi.value>-Inf), pi.value, NA) #sanity check
      if (algorithm.type==1) {
        p.1 <- dgamma(x=Smin, shape=N*theta[1] + am., rate=bm., log=TRUE)
        p.2 <- dbinom(n, size=N, prob=pi.value, log=TRUE)  #vector of log(target-distr) evaluated at theta.t,S.t,N.t
      } else if (algorithm.type==2) {
        p.1 <- dgamma(x=Smin, shape=n*theta[1] + am., rate=bm., log=TRUE)
        p.2 <- (N-n)*log(1-pi.value)  
        if(N==n){
          p.2 <- 0
        }
      }
      p <- p.1 + p.2
      
    } else {
      # Avoid bunch of annoying dpareto warnings...
      pi.value <- 0
      p <- -Inf
      p.1 <- p.2 <- NA
    }
    p.jacobian <- z
    
    log.target <- p + p.jacobian
    
    if (verbose>1) {
      #print above results (error checking)
      cat("log of pi.value:\n"); print(log(pi.value))
      cat("p.1:\n"); print(p.1)
      cat("p.2:\n"); print(p.2)
      cat("p.jacobian:\n"); print(p.jacobian)
      cat("log target density, p:\n"); print(p)
      cat("z:\n"); print(z)
      cat("Smin:\n"); print(Smin) 
      cat("min(S_i):\n"); print(upperBound)
      cat("bp:\n"); print(bp)
      if (algorithm.type==2 && is.na(p.1)) {
        cat("p1:gamma:\n"); print( dgamma(x=Smin, shape=n*theta[1] + am., rate=bm., log=TRUE) )
        cat("p2:log(1-pi)^(N-n):\n"); print( (N-n)*log(1-pi.value) )
      }
    }
    
    return(log.target)
    
  } # END log.target.fun
  #####################################
  
  Smin.curr <- Smin.t
  z.curr <- log(Smin.curr)
  
  if (verbose) {
    cat("Inside updating step for Smin.t(bp.t)     #before update\n")
    cat("Smin.t:\n"); print(Smin.t)          #before update
    cat("z.curr:\n"); print(z.curr)          #before update
    cat("v.sm:\n"); print(v.sm)
    cat("min{S}:\n"); print(minSs)
    cat("theta.t:\n"); print(theta.t)
    cat(ppaste("N.t: ",N.t,", n: ",n.t,"\n"))
    cat(ppaste("am: ",am,", bm: ",bm,"\n"))
    
  }  
  
  # Run various sampling algorithms of Metropolis alternatives
  if (met.alg=="MH") {
    # Metropolis algorithm for sampling from p(z.t| .) for z(Smin).
    tmp <- alg.met(x.curr=z.curr, r.q.prop.fun=r.q.prop.fun.mtm, 
                   log.target.fun=log.target.fun, verbose=verbose)
    
  } else if (met.alg=="MTM") {
    # Multiple-Try Metropolis algorithm for sampling from p(z.t| .) for whole vector of z(bp), i=0,1,...,m  with bp_0=Smin.
    tmp <- alg.mtm(x.curr=z.curr, k=kk, r.q.prop.fun=r.q.prop.fun.mtm, d.q.prop.fun=d.q.prop.fun.mtm,
                   log.target.fun=log.target.fun, verbose=verbose)
    
  } # END Metropolis-algorithm-alternatives
  
  got.draw.idx <- tmp$got.draw.idx
  if (got.draw.idx) {
    # Retrieve results
    z.prop <- tmp$y.prop

    if(verbose) {
      cat("Selected Smin.t:\n"); print(exp(z.prop))          #before update
    }
    
    # Update
    z.t <- z.prop     #update accepted proposals, iteration t=iter+1
    Smin.t <- exp(z.t)    #Smin transformation   
  }    
  if (verbose) {
    cat("Smin.t after update:\n"); print(Smin.t)          #before update
  } 
  
#   }
  
  return(list('Smin.t'=Smin.t,"idx"=got.draw.idx))
}

met.smin <- cmpfun(met.smin)
