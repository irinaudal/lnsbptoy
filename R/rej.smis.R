"rej.smis" <- function(S.mis.t, theta.t, N.t, n, Smin.t, bp.t, gamma, E, g,
                      n.tries.max=5000, verbose=FALSE ){
  
  ####################################################################################################
  # rej.smis   Rejection sampling step for single draw from posterior for S.mis
  #
  # Input: S.mis.t = previous iteration of S.mis
  #        theta.t = theta of current iteration
  #        N.t     = N of current iteration
  #        n       = number of missing sources
  #        Smin    = minimum flux the sources can be detected to, hyper-parameter in pareto prior for S
  #        p.t     = previous vector of mixture component probabilities
  #        I.idx.t = previous vector of mixture component indices
  #        pble    = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        gamma   = constant of transformation, energy per photon
  #        g       = function, probability of observing a source
  #        g.type  = type of g-function {"step","smooth","table"}
  #        use.bp  = (T/F) break-point pareto version
  #        use.mix = (T/F) mixture pareto version
  #        n.tries.max = max number of tries to propose new vectors of S.mis.t
  #        verbose = (T/F) display progress of program
  #
  # Output: S.mis.t = posterior samples of parameter: flux of missing sources
  #         idx     = indices of accepted proposals
  ####################################################################################################
  
  # Pr(obs source) via select lowest possible setting of parameters based on pble
  prob.obs.source <- 1-g(lambda=Smin.t*E/gamma)
  if (verbose) {
    cat("S.mis.t:     #before update\n"); print(S.mis.t)          #before update
    cat(paste("N.t = ",N.t,", n = ",n,"\n",sep=""))
    cat(paste("E = ",E,"\n",sep=""))
  }   
  
  if (N.t==n) {
    S.mis.t <- numeric(0)
    idx <- 0
  } else {
    
    S.mis.prop <- rbrokenpareto(n=N.t-n,x_min=Smin.t,k=theta.t,bp=bp.t)   #vector of proposals of S.mis    
    
    lambda.prop <- S.mis.prop*E/gamma           #vector of proposed lambda for missing cources, not saved
    lambda.prop[lambda.prop <= 0] <- NA             #set the negative lambda proposals to NA
    g.prop <- 1-g(lambda=lambda.prop)   #vector of proposed prob. of missing sources, not saved  
    
    log.alpha <- log(g.prop)              #compare proposal to current
    log.u <- log(runif(N.t-n))            #random vector U(0,1)
    idx <- (log.u < log.alpha)            #accepted proposal indicator 
    idx[is.na(idx)] <- FALSE
    
    #NOTE: MUST return completely new vector of S.mis.t
    #      so that all S.mis proposals are eventually excepted
    
    #fail safe mechanism
    failed.indices <- (!idx)
    n.tries <- 0
    
    if (verbose) {
      cat("S.mis.prop = \n"); print(S.mis.prop)
      cat("failed.indices = \n"); print(failed.indices)
    }
    
    while (any(failed.indices) && n.tries<n.tries.max) {
      
      if (verbose) {
        cat("Printing TEST in rej.smis.R ...\n")
      }      
      
      n.tries <- n.tries+1     
      # In this case, we have to propose more S.mis:
      N.fail <- sum(failed.indices)
      S.mis.prop[failed.indices] <- rbrokenpareto(n=N.fail,x_min=Smin.t,k=theta.t,bp=bp.t)  #vector of proposals of S.mis    
      
      lambda.prop <- S.mis.prop[failed.indices]*E/gamma           #vector of proposed lambda for missing cources, not saved
      lambda.prop[lambda.prop <= 0] <- NA             #set the negative lambda proposals to zero
      g.prop <- 1-g(lambda=lambda.prop)   #vector of proposed prob. of missing sources, not saved  
      log.alpha <- log(g.prop)              #compare proposal to current
      log.u <- log(runif(N.fail))            #random vector U(0,1)
      idx[failed.indices] <- (log.u < log.alpha)            #accepted proposal indicator 
      idx[is.na(idx)] <- FALSE
            
      if (verbose){
        cat(ppaste("This is try #",n.tries,"\n"))
        cat(ppaste("There were ",N.fail," failed indices...\n"))
        cat("New proposals for S.mis:\n"); print(S.mis.prop[failed.indices])
        cat("Smin = \n"); print(Smin.t)
        cat("theta.t = \n"); print(theta.t)
        cat("E.mis = \n"); print(E)
        cat("lambda.prop = \n"); print(lambda.prop)
        cat("Pr(observe source) = \n"); print(1-g.prop)
        cat("g.prop = \n"); print(g.prop)
        cat("log.alpha = \n"); print(log.alpha)
        cat("log.u = \n"); print(log.u)
        cat("idx (accepted proposals)  = \n"); print(idx)
        cat("failed.indices (rejected) = \n"); print(!idx)
      }
      
      # Check if S.mis is still invalid...
      failed.indices <- (!idx)
      
    }
    
    # Error catch after fail safe mechanism
    if(n.tries>=n.tries.max){
      stop(paste("Error: Cannot draw acceptible S.mis values. Check [new proposals of S.mis]. 
                 The current value of Smin = ",Smin.t, ", N.t=",N.t, ", n=",n,"\n", sep=""))
    } else {
      S.mis.t <- S.mis.prop       #update accepted proposals, iteration t=iter+1, NOTE: size will change!
    }
    
    if (verbose==TRUE) {
      # Print above results (error checking)
      cat("Rej.Sampling for S.mis took this many iterations:\n"); print(n.tries)
      cat("theta.t:\n"); print(theta.t)
      cat("N.t:\n"); print(N.t)
      cat("lambda.prop:\n"); print(lambda.prop)
      cat("log of g.prop:\n"); print(log(g.prop))
      cat("S.mis.t:     #after update\n"); print(S.mis.t)          #after update
      cat("S.mis.prop:\n"); print(S.mis.prop)  
      cat("log.alpha:\n"); print(log.alpha)
      cat("log.u:\n"); print(log.u)
      cat("idx (accepted proposals):\n"); print(idx)
    }
  }    
  return(list('S.mis.t'=S.mis.t,'idx'=idx))  
}

rej.smis <- cmpfun(rej.smis)
