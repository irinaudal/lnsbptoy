"patch.smis" <- function(S.mis.t, theta.t, Smin, E.mis, L.mis, bg.mis, gamma, pble, g.type, g, I.idx.t=NULL, bp=NULL, 
                         n.tries.max=5000, use.bp=FALSE, use.mix=FALSE, verbose=FALSE){
  
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
  lambda.curr <- S.mis.t*E.mis/gamma         #vector of lambda for missing cources, not saved  
  
  g.curr <- g(lambda=lambda.curr,bg=bg.mis,E=E.mis,L=L.mis,g.type=g.type)   #vector of prob. of observing sources, not saved  
  if (verbose){
    cat("Patching up missing starting states...\n")
    cat("E.mis = \n")
    print(E.mis)
    cat("L.mis = \n")
    print(L.mis)
    cat("bg.mis = \n")
    print(bg.mis)
  }
  
  n.tries <- 0
  
  while (any(g.curr==1) & n.tries<n.tries.max){
    # Patch up the impossible starting states...
    failed.indices <- (g.curr==1)
    N.failed.indices <- sum(failed.indices)
    par <- sample.pble(pble,N.failed.indices)
    bg.mis <- par$B  #vector of new values of B.mis, not saved
    L.mis <- par$L  #vector of new values of L.mis, not saved
    E.mis <- par$E  #vector of new values of E.mis, not saved
    
    # Find the value of S_{upper} s.t. g(lambda(S_{upper}),E) < 0.9
    # Then generate from a pareto with lower limit S_min, but now upper-truncated at S.upper...
    if(use.bp){
      S.mis.t[failed.indices] <- rbrokenpareto(n=N.failed.indices,x_min=Smin,k=theta.t,bp=bp,verbose=verbose)     #vector of initial values: S.obs0
    } else if(use.mix){        
      S.mis.t[failed.indices] <- rpareto(n=N.failed.indices,x_min=Smin[I.idx.t[failed.indices]],k=theta.t[I.idx.t[failed.indices]])
    } else{ #No mixtures case
      S.mis.t[failed.indices] <- rpareto(n=N.failed.indices,x_min=Smin,k=theta.t)
    }
    lambda.curr <- S.mis.t[failed.indices]*E.mis/gamma
    g.curr[failed.indices] <- g(lambda=lambda.curr, bg=bg.mis, E=E.mis, L=L.mis, g.type=g.type)
    if (verbose){
      cat("Debugging patch.smis ... \n")
      cat(ppaste("This is try #",n.tries,"\n"))        
      cat(ppaste("There were ",N.failed.indices," failed indices...\n"))
      cat("g.curr(all) = \n")
      print(g.curr)
      cat("good g.curr (slightly <1): residual is = \n")
      print(1-g.curr[g.curr<1])
      cat("failed.indices (rejected) = \n"); print(failed.indices)        
      cat("New [failed] proposals for S.mis:\n")
      print(S.mis.t[failed.indices])
      cat("I.idx.t[failed] = \n")
      print(I.idx.t[failed.indices])
      cat("Smin[I.idx.t[failed]] = \n")
      print(Smin[I.idx.t[failed.indices]])
      cat("theta.t[I.idx.t[failed]] = \n")
      print(theta.t[I.idx.t[failed.indices]])
      cat("E.mis[failed] = \n")
      print(E.mis)
      cat("bg.mis[failed] = \n")
      print(bg.mis)
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
