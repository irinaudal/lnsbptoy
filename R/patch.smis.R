"patch.smis" <- function(S.mis.t, theta.t, Smin, bp, E, gamma, g,
                         n.tries.max=5000, verbose=FALSE){
  
  ####################################################################################################
  # patch.smis   Patch S.missing starting values, if needed
  #
  # Input: S.mis.t = initial value of S.mis to be patched
  #        theta.t = initial value of theta
  #        I.idx.t = initial value of mixture indices I
  #        Smin    = minimum flux the sources can be detected to, hyper-parameter in pareto prior for S
  #        E.mis   = exposure for missing sources
  #        L.mis   = location for missing sources
  #        bg.mis  = background counts/pixel for missing sources
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
  # Output: S.mis.t = patched missing starting values
  ####################################################################################################
  
  # Check starting states of S.mis are all valid...
  lambda.curr <- S.mis.t*E/gamma         #vector of lambda for missing cources, not saved  
  g.curr <- g(lambda=lambda.curr)   #vector of prob. of observing sources, not saved  

  n.tries <- 0
  
  while (any(g.curr==1) & n.tries<n.tries.max){
    # Patch up the impossible starting states...
    failed.indices <- (g.curr==1)
    N.failed.indices <- sum(failed.indices)
    
    # Find the value of S_{upper} s.t. g(lambda(S_{upper}),E) < 0.9
    # Then generate from a pareto with lower limit S_min, but now upper-truncated at S.upper...
    S.mis.t[failed.indices] <- rbrokenpareto(n=N.failed.indices,x_min=Smin,k=theta.t,bp=bp,verbose=verbose)     #vector of initial values: S.obs0
    lambda.curr <- S.mis.t[failed.indices]*E/gamma
    g.curr[failed.indices] <- g(lambda=lambda.curr)
    
    if (verbose){
      cat("Debugging patch.smis ... \n")
      cat(paste("This is try #",n.tries,"\n"))        
      cat(paste("There were",N.failed.indices,"failed indices...\n"))
      cat("g.curr(all) = \n")
      print(g.curr)
      cat("good g.curr (slightly <1): residual is = \n")
      print(1-g.curr[g.curr<1])
      cat("failed.indices (rejected) = \n"); print(failed.indices)        
      cat("New [failed] proposals for S.mis:\n")
      print(S.mis.t[failed.indices])
      cat("Smin = \n")
      print(Smin)
      cat("theta.t = \n")
      print(theta.t)
      cat("E = \n")
      print(E)
      cat("lambda.curr[failed] = \n")
      print(lambda.curr)
      cat("g.curr[failed] = \n")
      print(g.curr[failed.indices])
    }
    n.tries <- n.tries+1
  }
  
  if(n.tries==n.tries.max & all(g.curr==1)){
    S.mis.t <- NULL
  } else if(n.tries==n.tries.max){
    stop("Error: Reached maximum number of tries to give valid proposals for S.mis.")
  } 
  
  return(S.mis.t)
}

patch.smis <- cmpfun(patch.smis)
