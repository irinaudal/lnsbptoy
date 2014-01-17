"pi.theta.eval.mc" <- function(nsamples=10000, theta.grid=seq(0.1,2.1,by=0.01), Smin.grid=10^-17, bp.grid=10^-16, p.t.grid=1,
                               m=1, gamma, pble, g.type="smooth", model="regular",
                               g=function(lambda,bg,E,L,g.type){
                                 g.compute(lambda=lambda,bg=bg,E=E,L=L,g.type=g.type)
                               }, use.bp=TRUE, use.mix=FALSE, verbose=FALSE ){
  
  ####################################################################################################
  # pi.theta.eval.mc   Function evaluates pi(theta,Smin) integral via importance sampling / Monte Carlo integration
  #
  # Input: nsamples   = number of samples to draw for MC
  #        theta.grid = theta values (fixed) for which pi(theta,Smin) is to be evaluated
  #        Smin.grid  = Smin values (fixed) for which pi(theta,Smin) is to be evaluated
  #        bp.grid    = either NULL (no fixed break-points), or, a vector of length p specifying the break-points
  #        p.t.grid   = first mixture proportion of p=(p1,1-p1) for mixture model (fixed) for which pi(theta,Smin) is to be evaluated
  #        m     = dimension of theta, number of Pareto populations
  #        d     = vector if hyper-parameter(s) in dirichlet prior for mixture proportions (p_0,...,p_m)) (broken power-law)
  #        gamma      = constant of transformation, energy per photon
  #        pble       = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        g.type     = type of g-function {"step","smooth","table"}
  #        g          = function, probability of observing a source
  #        verbose    = (T/F) display progress of program
  #
  # Output: pi.theta.approx = MC approximation to pi(theta,Smin) integral
  ####################################################################################################
  
  ### TODO: add code to support bp and mix models
  if (use.mix) {
    stop("Cannot compute pi(theta,Smin) via Monte Carlo for 'mix' model - to be implemented.")
  }
  
  if (verbose) {
    cat("Computing pi(theta,Smin) integral via Monte Carlo.. \n") 
  }
  
  verbose3 <- verbose>2
  if (verbose) {
    cat("class(pble): \n")
    print(class(pble))
    cat("g.type: \n")
    print(g.type)
    cat("model: \n")
    if (use.bp) { cat("bp \n")} else if (use.mix) { cat("mix \n")} else { cat("regular \n")}
    if (verbose3) {
      cat("bp.grid:\n"); print(bp.grid) 
    }
  }
  
  # Create pble object, if needed
  if (class(pble)=="character") {
    pble <- new.pble.table(effects.file=pble)
  } else if (class(pble)!="pble.basic" && class(pble)!="pble.table"){
    stop("'pble' could not be created.")
  }  
  
  # Generate the total number of sources:
  N <- nsamples
  
  if(is.vector(theta.grid)){
    theta.grid <- matrix(theta.grid,ncol=m)
  }
  if(is.vector(bp.grid)){
    bp.grid <- matrix(bp.grid,ncol=m-1)
  }
  
  # Generate Smin, the minimum source flux ? NO. Use Smin.grid
  length.theta <- nrow(theta.grid)
  
  if (length(Smin.grid)==1) {
    Smin.grid <- rep(Smin.grid,length.theta)
  } 
  
  if (!use.bp){
    bp.grid <- NULL
  } 
  
  # Mixture-pareto:
  # Generate mixture proportion(s) ? NO.  Use p.t.grid  
  
  ################## ################## ###############
  # Create storage vector
  pi <- rep(NA,times=length.theta)
  
  for (i in 1:length.theta) {
    
    theta <- theta.grid[i,]
    Smin  <- Smin.grid[i]
    if (use.bp) {
      bp    <- bp.grid[i,]
    }
    if (use.mix) {
      p   <- p.t.grid[i,1]
    }
    # sanity check
    if (use.bp) {
      if (any(Smin > bp)) {
        pi[i] <- NA
        next
      }
    }
    
    # Generate source flux S:
    ## Generalizing to handle both break-points and untruncated mixtures
    if (use.bp){      
      # Generate from broken power-law mixture of Truncated Pareto's:
      S     <- rbrokenpareto(n=N,x_min=Smin,k=theta,bp=bp,verbose=verbose3)
      #I.idx <- NULL
    } else if (use.mix){
      # Generate from mixture of Pareto's:
      tmp   <- rmixpareto(n=N,p=p,k=theta,x_min=Smin) #this function is already compatible with single Pareto case
      #I.idx <- tmp$I.idx
      S     <- tmp$S
    } else {
      # Generate from single Pareto:
      S     <- rpareto(n=N,k=theta,x_min=Smin)
      #I.idx <- NULL
    }
    idx <- S==Inf
    if(any(idx)){
      S[idx] <- runif(sum(idx),1.5e+268,9e+297)
    }
    # Generate standard background, effective area etc.:  
    par <- sample.pble(pble,N)
    bg  <- par$B
    L   <- par$L
    E   <- par$E      
    
    # Generate intensity of photon counts:
    lambda.src <- S*E/gamma
    
    # Determine whether observed or missing:
    gp <- g(lambda=lambda.src,bg=bg,E=E,L=L,g.type=g.type,verbose=verbose3)
    
    J <- rbinom(N,1,gp)    #length of gp is N
    n <- sum(J)
    
    ################## ################## ###############
    # Store proportion of observed sources
    pi[i] <- n/N
    
    
    if (verbose) {
      if (use.bp) {
        cat(c("Computed pi(theta,Smin,bp) for Smin =",ppaste(signif(Smin,5)),"; bp =",ppaste(signif(bp,5)),"; theta = ",ppaste(round(theta,5))," (",i,"/",length.theta,") ; Pr(Obs) = ",round(pi[i],5),"\n"))            
      } else if (use.mix){
        cat(c("Computed pi(theta,Smin) for Smin =",ppaste(signif(Smin,5)),"; theta = ",ppaste(round(theta,5))," (",i,"/",length.theta,") ; p.t = ",ppaste(p[i,]), " ; Pr(Obs) = ",round(pi[i],5),"\n"))            
      } else {
        cat(paste("Computed pi(theta,Smin) for Smin =",ppaste(signif(Smin,5)),"; theta = ",round(theta,5)," (",i,"/",length.theta,") ; Pr(Obs) = ",round(pi[i],5),"\n",sep=""))      
      }
    }
    
  }
  
  ################## ################## ###############
  # Evaluate proportion of observed sources
  pi.theta.approx <- list("theta"=theta.grid,"pi.theta"=pi)
  
  return(pi.theta.approx)
  
}
