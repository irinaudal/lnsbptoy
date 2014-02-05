"patch.sobs" <- function(S.obs.t, theta.t, Smin, bp, E, gamma, g, 
                         n.tries.max=1000, verbose=FALSE){
  
  ####################################################################################################
  # patch.sobs   Patch S.observed starting values, if needed
  #
  # Input: S.obs.t = initial value of S.obs to be patched
  #        theta.t = initial value of theta
  #        I.idx.t = initial value of mixture indices I
  #        Smin    = minimum flux the sources can be detected to, hyper-parameter in pareto prior for S
  #        E.obs   = exposure for observed sources
  #        L.obs   = location for observed sources
  #        bg.obs  = background counts/pixel for observed sources
  #        gamma   = constant of transformation, energy per photon
  #        pble    = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        bp      = vector of break-points for BrokenPareto model
  #        g.type  = type of g-function {"step","smooth","table"}
  #        g       = function, probability of observing a source
  #        n.tries.max = max number patches, of while-loop iterations
  #        use.bp  = (T/F) break-point pareto version
  #        use.mix = (T/F) mixture pareto version
  #        verbose = (T/F) display progress of program
  #
  # Output: S.obs.t = patched observed starting values
  ####################################################################################################
  
  
  # Check starting states of S.obs are all valid...
  if (any(S.obs.t<Smin)) {
    S.obs.t[S.obs.t<Smin] <- Smin*1.0001
  }
  
  lambda.curr <- S.obs.t*E/gamma         #vector of lambda for observed cources, not saved  
  g.curr <- g(lambda=lambda.curr)    #vector of prob. of observing sources, not saved  
  
  # If breakpoints are specified then the minimum value can't be increased above the lowest breakpoint:
  S.min.range <- Smin
  if(!is.null(bp)) {
    S.min.range <- min(bp)-Smin
  }
  if (any(S.min.range<0)){
    stop("'bp' less than 'Smin' detected in 'patch.sobs'")
  }
  # Increase 10% of the way to the maximum value it can be:
  S.star <- Smin + S.min.range*0.1
  
  if (verbose){
    cat("Patching up observed starting states...\n")
    cat("Failed indices vector for starting states of S.obs (first attempt):\n")
    print(g.curr==0)
  }
  
  n.tries <- 0
  
  while (any(g.curr==0) & n.tries<n.tries.max){
    # Patch up the impossible starting states...
    failed.indices <- (g.curr==0)
    N.failed.indices <- sum(failed.indices)
    
    if (is.null(bp)) {
      # Find the value of s_* s.t. g(lambda(s_*),E) > 0.1
      S.star <- S.star + S.min.range*0.1
    } else {
      # Find the value of s_* s.t. g(lambda(s_*),E) > 0.1
      if (S.star + S.min.range*0.1 < min(bp)){
        S.star <- S.star + S.min.range*0.1
      }
    }
    # Then generate from a pareto with lower limit s_* instead of S_min...
    S.obs.t[failed.indices] <- rbrokenpareto(n=N.failed.indices,x_min=S.star,k=theta.t,bp=bp,verbose=verbose) 
    lambda.curr <- S.obs.t[failed.indices]*E/gamma
    g.curr[failed.indices] <- g(lambda=lambda.curr)
    if (verbose){
      cat("Failed indices vector for starting states of S.obs:\n")
      print(failed.indices)
      cat("Adjusted (lambda,bg,E,L):\n")
      print(cbind(lambda.curr,E))
    } 
    n.tries <- n.tries+1
  }
  
  if(n.tries==n.tries.max){
    stop("Error: Reached maximum number of tries to give valid proposals for S.obs.")
  }
  
  return(S.obs.t)
}

patch.sobs <- cmpfun(patch.sobs)
