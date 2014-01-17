"numint.N" <- function(N.t, theta.t, n, alpha, beta, pi, 
                       Smin, gamma, pble, length.S, fixed.S.pi, debug.pi, g.type, g, nsamples,
                       bp=NULL, p.t=NULL, tol=10^-10, use.bp=FALSE, use.mix=FALSE, verbose=FALSE){
  
  ####################################################################################################
  # numint.N   Numerical Integration step for single draw from posterior for N
  #
  # Input: N.t     = previous iteration of N
  #        theta.t = current iteration of theta
  #        n       = observed number of sources 
  #        alpha   = target number of sccessful trials in negbinom prior for N
  #        beta    = dispersion parameter in negbinom prior for N
  #        Smin    = minimum flux the sources can be detected to according to E.obs and g.
  #        gamma   = constant of transformation, energy per photon
  #        pble    = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        bp      = either NULL (no fixed break-points), or, a vector of length p specifying the break-points
  #        p.t     = previous vector of mixture component probabilities
  #        length.S = length of grid of S, for integration of g(theta)
  #        pi       = R.object: a list("pi","theta") with detection probalities: pi(theta)
  #        fixed.S.pi = (override) fix the value of S (used mainly for debugging)
  #        nsamples   = number of MC samples to produce estimate of pi(theta,Smin) integral
  #        debug.pi = (T/F) (used mainly for debugging)
  #        g.type   = type of g-function {"step","smooth","table"}
  #        g        = function, probability of observing source
  #        use.bp   = (T/F) break-point pareto version
  #        use.mix  = (T/F) mixture pareto version
  #        verbose  = (T/F) display progress of program
  #
  # Output: N.t  = posterior sample of parameter N
  ####################################################################################################
  
  verbose2 <- 0
  if (verbose>1) {
    verbose2 <- 1
  }
  
  pi.value <- pi.theta.get(pi=pi, theta=theta.t, p.t=p.t, bp=bp, Smin=Smin, gamma=gamma, 
                           pble=pble, length.S=length.S, fixed.S=fixed.S.pi, debug=debug.pi,
                           nsamples=nsamples, g.type=g.type, g=g, use.bp=use.bp, use.mix=use.mix, verbose=verbose2)     #marginal prob. of observing sources  
  
  "log.f" <- function(N,n,a,b,theta,pi.value){
    #vectorized log-pmf, unnormalized:
    return( lgamma(N+a)-lgamma(N-n+1) + (N-n)*log(1.0-pi.value) - N*log(1+b) )  
  }
  
  if (verbose) {
    cat("N.t:     #before update\n"); print(N.t)          #before update
  }     
  
  #check for special case where pi(theta) = 1,
  #in this case, there can be no missing data...
  if (1.0-pi.value <= 0.0){
    if (verbose){
      cat("pi(theta) = 1.0 => Returning N.t = n\n")
    }
    return(list("N.t"=n))
  }
  
  #find range for N for function evaluation
  #note, N.max cannot be smaller than n
  if (pi.value==0) {
    N.max <- n*10000  # make sure we get the far tail of the distribution  
  } else {
    N.max <- max(n,ceiling(n/pi.value)) + 50  
  }
  N.grid <- n:max(n,N.max)
  
  #compute proportional pmf
  f.grid <- log.f(N=N.grid,n=n,a=alpha,b=beta,theta=theta.t,pi.value=pi.value)
  f.unnorm.grid <- f.grid
  
  adj.const <- max(f.grid)   #constant to help numerical balance
  f.grid <- exp( f.grid - adj.const )
  
  #compute normalized pmf  
  pmf <- f.grid/sum(f.grid)
  
  if (is.na(pmf[length(N.grid)])){
    if (verbose){
      cat("\n### NA's found in CDF of N -- BEGIN Debugging Information ###\n")
      cat("pi.value:\n")
      print(pi.value)
      cat("N.grid:\n")
      print(N.grid)
      cat(paste("(n=",n,"theta=",theta.t,")\n",sep=""))
      cat(paste("pi(theta) = ",pi.value,"\n",sep=""))
      cat("f.unnorm.grid:\n")
      print(f.unnorm.grid)
      cat(paste("adjustment constant = ",adj.const,"\n",sep=""))
      cat("pmf:\n")
      print(pmf)
      cat("### NA's found in CDF of N -- END Debugging Information ###\n\n")
    }
    stop("'NA' found in CDF of N computed using numerical integration")
  }
  
  # Check that the last value has a small enough probability:
  j <- 0
  while (pmf[length(N.grid)] > tol) {
    j <- j+1
    if(j>100) {
      stop("Warning: cannot find NegBin pdf < Tolerance")
    } 
    cur.N.max <- N.grid[length(N.grid)]
    extra.N.grid <- (cur.N.max+1):(cur.N.max+50)
    #compute proportional pmf, append extra function evaluations to the tail of f.grid      
    extra.f.grid <- exp( log.f(N=extra.N.grid,n=n,a=alpha,b=beta,theta=theta.t,pi.value=pi.value) - adj.const )
    f.grid <- c(f.grid,extra.f.grid)
    N.grid <- c(N.grid,extra.N.grid)      
    #compute normalized pmf  
    pmf <- f.grid/sum(f.grid)
    if (verbose){
      cat("sum(f.grid), last.value.of.pmf, n, theta.t: \n");
      print(c(sum(f.grid), pmf[length(N.grid)],n,theta.t)) 
    }
  }
  
  #get cdf of N.t
  F <- c(0, cumsum(pmf) )
  
  #take a random sample from pmf of N
  u <- runif(1)
  #find first cdf value that exceeds u, and subtract 1 to find largest less than it
  idx <- which(u <= F)[1]-1 
  
  #first value of F is zero, so idx should never be zero! check just in case...
  #although, the while loop should never really be needed...
  j <- 0
  while (idx==0) {
    j <- j+1
    if(j>100) {
      stop("Warning: NegBin pdf is too narrow!")
    }
    u <- runif(1)
    idx <- which(u <= F)[1]-1 
  } 
  
  #assume was well, now pick off the actual draw.
  #recall that the lowest possible idx is 1, in which case N.t=n
  N.t <- N.grid[idx]
  
  if (verbose){
    cat("_Num.Intgr. for N.t_  ::    \n");cat("log(f.unnorm.grid): \n");
    print(round(f.unnorm.grid,2))
    cat("pi.value:\n")
    print(pi.value)
    cat("f.grid: \n");
    print(round(f.grid,2))
    cat("N.grid: \n");
    print(N.grid)
    cat("pmf: \n");
    print(round(pmf,2))
    cat("cdf: \n");
    print(round(F,4))
    cat(paste("Unif(0,1) = ",u,"\n",sep=""))
    cat(paste("Draw      = ",N.t,"\n",sep=""))
    cat("N.t:      #after update\n"); print(N.t)        #after update
    cat("N.max:\n"); print(N.max)
    cat("corrected.N.max:\n"); print(N.grid[length(N.grid)])
    cat("unif:\n"); print(u)
    cat("idx:\n"); print(idx)
  } 
  
  return(list("N.t"=N.t))  
}

numint.N <- cmpfun(numint.N)
