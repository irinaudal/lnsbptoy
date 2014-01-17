"met.bp" <- function(bp.t, Smin.t, S.obs.t, S.mis.t, theta.t, N.t, n.t, v.bp, z.e,
                     C, mu, gamma, pble, length.S, g.type, g, pi, nsamples,
                     fixed.bp=NULL,  
                     algorithm.type=2,
                     met.alg="MH", kk=5, # "MTM" or "MH" or "HitRun"
                     TOL=1e-10,
                     fixed.S.pi=NULL, debug.pi=FALSE, use.bp=TRUE, use.mix=FALSE, 
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
  
  if (met.alg=="HitRun") {
    stop("HitRun is not implemented yet. Also, cannot use this alg. to sample 1-D parameter.")
  }
  if (met.alg=="MTM" && kk==1) {
    met.alg <- "MH"
  }
  
  if (verbose) {
    cat("bp.t:     #before update\n"); print(bp.t)          #before update
  }  
  
  n <- n.t
  K.t <- c(Smin.t,bp.t)
  m.bp <- length(bp.t)
  
  if (algorithm.type==1) {
    S.t <- c(S.obs.t,S.mis.t)
  } else if (algorithm.type==2) {
    S.t <- S.obs.t
  }
  
  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  
  #####################################
  ### Define utility functions ###
  #####################################
  
  
#   "r.q.prop.fun.hitrun" <- function(x.curr, e, v=v.par[1], m=m.par) {
#     rj <- rnorm(n=m,0,v) # Note: must be a univariate proposal
#     y.prop <- x.curr + rj*e
#     return(y.prop)
#   }
#   "d.q.prop.fun.hitrun" <- function(x.curr, y.prop, v=v.par[1], m=m.par) {
#     # q.curr.to.prop: Normal density in 'r'-space: Te(x,y): prob of y given x.
#     d.prop <- sqrt(sum((x.curr-y.prop)^2))
#     return(dnorm(x=d.prop, mean=0, sd=v, log = FALSE))
#   }
#   "r.vec.fun" <- function(z=z.e, m=m.par, method=2) {
#     if (is.null(z)) {
#       if (method==1 && m==2) { # omega method
#         Pi <- acos(0)*2
#         omega <- runif(1, 0, 2*Pi)
#         e <- c(cos(omega), sin(omega))
#       } else { # U(0,1) method
#         z <- runif(n=1)
#         e <- c(1-z,rep(z/m,times=m-1))
#       }
#     } else {
#       e <- c(1-z,rep(z/m,times=m-1))
#     }
#     
#     # normalize to unit vector
#     e <- e/sqrt(sum(e^2))
#     return(e)
#   }
  
  "r.q.prop.fun.mtm" <- function(x.curr, v=v.bp, m=m.bp) {
    return(rnorm(n=m, mean=x.curr, sd=v))
  }
  "d.q.prop.fun.mtm" <- function(x.curr, y.prop) {
    return(0)
  }

  "transform.update.variables" <- function(z, Smin.t, fixed.bp, TOL){

    # Transform variables
    bp   <- cumsum(exp(z)) + Smin.t
    K    <- c(Smin.t,bp)
    z.after <- log(diff(K))
    
    # Overwrite/take care of fixed.bp
    fixed.idx  <- (fixed.bp != FALSE)
    # Check for numerical accuracy
    if (abs(z-z.after)>TOL) {
      warning("Proposal of z(bp,Smin) is too small, numerical error incured with transformation.")
      bp <- -Inf
    } else if (any(fixed.idx)) {
      bp[fixed.idx] <- fixed.bp[fixed.idx]  # overwrite with fixed values 
      K <- c(Smin.t,bp)
      z <- log(diff(K))
    }
    return(list("z"=z, "K"=K, "bp"=bp))
  } # END transform
  
  "log.target.fun" <- function(x, Smin=Smin.t, fixed.Kj=fixed.bp, TOL.=TOL, 
                               pi.vec=pi, theta=theta.t, N=N.t, n=n.t, S=S.t, mu.=mu, C.=C,...) {
    
    # Transform variables
    tmp <- transform.update.variables(z=x, Smin.t=Smin, fixed.bp=fixed.Kj, TOL=TOL.)
    z    <- tmp$z
    K    <- tmp$K
    bp   <- tmp$bp
    
    if (all(Smin < bp)) { # check boundary conditions
      
      pi.value <- pi.theta.get(pi=pi.vec, theta=theta, p.t=NULL, bp=bp, Smin=Smin, gamma=gamma, 
                               pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi,
                               nsamples=nsamples, g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources
      pi.value <- ifelse(any(pi.value>-Inf), pi.value, NA) #sanity check
      
      p.1 <- sum(dnorm(x=z, mean=mu., sd=C., log=TRUE)) # prior, no Jacobian
      p.2 <- sum(dbrokenpareto(S, x_min=Smin, k=theta, bp=bp, log=TRUE)) # likelihood
      if (algorithm.type==1) {
        p.3 <- dbinom(n, size=N, prob=pi.value, log=TRUE)
      } else if (algorithm.type==2) {
        p.3 <- (N-n)*log(1-pi.value)
        if(N.t==n){
          p.3 <- 0
        }
      }
      log.target <- p.1 + p.2 + p.3
    } else {
      # Avoid bunch of annoying dpareto warnings...
      pi.value <- 0
      log.target <- -Inf
      p.1 <- p.2 <- p.3 <- NA
    }
    
    if (verbose>1) {
      #print above results (error checking)
      cat("log of pi.value:\n"); print(log(pi.value))
      cat("z:\n"); print(z)
      cat("p.1:\n"); print(p.1)
      cat("p.2:\n"); print(p.2)
      cat("p.3:\n"); print(p.3)
      cat("log.target:\n"); print(log.target)
    }
    
    return(log.target)
    
  } # END log.target.fun
  #####################################
  
  if (verbose) {
    cat("Inside updating step for z.t(bp.t)     #before update\n")
    cat("m.bp:\n"); print(m.bp)
    cat("v.bp:\n"); print(v.bp)
    cat("K.t:\n"); print(K.t)
    cat("theta.t:\n"); print(theta.t)
  }  
  
  # Evaluate current state
  z.curr <- log(diff(K.t))
  
  # Run various sampling algorithms of Metropolis alternatives
  if (met.alg=="MH") {
    # Metropolis algorithm for sampling from p(z.t| .) for whole vector of z(bp), i=0,1,...,m  with bp_0=Smin.
    # prior parameters: mu, C. Assume diagonal form or not? C[1,1] or C[1]... see p.curr: dnorm(.)
    tmp <- alg.met(x.curr=z.curr, r.q.prop.fun=r.q.prop.fun.mtm, 
                   log.target.fun=log.target.fun, verbose=verbose)
    
  } else if (met.alg=="MTM") {
    # Multiple-Try Metropolis algorithm for sampling from p(z.t| .) for whole vector of z(bp), i=0,1,...,m  with bp_0=Smin.
    tmp <- alg.mtm(x.curr=z.curr, k=kk, r.q.prop.fun=r.q.prop.fun.mtm, d.q.prop.fun=d.q.prop.fun.mtm,
                   log.target.fun=log.target.fun, verbose=verbose)
    
#   } else if (met.alg=="HitRun") {
#     # Hit-and-Run Metropolis algorithm for sampling from p(z.t| .) for whole vector of z(bp), i=0,1,...,m  with bp_0=Smin.
#     tmp <- alg.hitrun(x.curr=z.curr, k=kk, r.vec.fun=r.vec.fun,
#                       r.q.prop.fun=r.q.prop.fun.hitrun, 
#                       d.q.prop.fun=d.q.prop.fun.hitrun,
#                       log.target.fun=log.target.fun, verbose=verbose)
  } # END Metropolis-algorithm-alternatives
 
  got.draw.idx <- tmp$got.draw.idx
  if (got.draw.idx) {
    # Retrieve results
    z.prop <- tmp$y.prop
    
    # Update
    z.t  <- z.prop       #update accepted proposals, iteration t=iter+1
    bp.t <- cumsum(exp(z.t))+Smin.t    # whole bp vector transformation
  }
  
  if (bp.t>1) {
    stop(ppaste("Proposed BP was accepted unreasonably large:\n bp.new = ",bp.t))
  }
  
  return(list('Smin.t'=Smin.t, 'bp.t'=bp.t,"idx"=got.draw.idx))
}

met.bp <- cmpfun(met.bp)