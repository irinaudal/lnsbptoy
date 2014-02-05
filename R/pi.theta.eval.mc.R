"pi.theta.eval.mc" <- function(nsamples=10000, theta.grid=seq(0.1,2.1,by=0.01), Smin.grid=10^-17, bp.grid=10^-16,
                               m=1, gamma, E=E,
                               g=function(lambda,bg,E,L,g.type){
                                 g.compute(lambda=lambda,bg=bg,E=E,L=L,g.type=g.type)
                               }, verbose=FALSE ){
  
  ####################################################################################################
  # pi.theta.eval.mc   Function evaluates pi(theta,Smin) integral via importance sampling / Monte Carlo integration
  #
  # Input: nsamples   = number of samples to draw for MC
  #        theta.grid = theta values (fixed) for which pi(theta,Smin) is to be evaluated
  #        Smin.grid  = Smin values (fixed) for which pi(theta,Smin) is to be evaluated
  #        bp.grid    = either NULL (no fixed break-points), or, a vector of length p specifying the break-points
  #        m     = dimension of theta, number of Pareto populations
  #        d     = vector if hyper-parameter(s) in dirichlet prior for mixture proportions (p_0,...,p_m)) (broken power-law)
  #        gamma      = constant of transformation, energy per photon
  #        g          = function, probability of observing a source
  #        verbose    = (T/F) display progress of program
  #
  # Output: pi.theta.approx = MC approximation to pi(theta,Smin) integral
  ####################################################################################################
  
  if (verbose) {
    cat("Computing pi(theta,Smin) integral via Monte Carlo.. \n") 
  }
  
  verbose3 <- verbose>2
  if (verbose) {
    if (verbose3) {
      cat("bp.grid:\n"); print(bp.grid) 
    }
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
  
  ################## ################## ###############
  # Create storage vector
  pi <- rep(NA,times=length.theta)
  
  for (i in 1:length.theta) {
    
    theta <- theta.grid[i,]
    Smin  <- Smin.grid[i]
    bp    <- bp.grid[i,]
    
    # sanity check
    if (any(Smin > bp)) {
      pi[i] <- NA
      next
    }
    
    # Generate source flux S:
    ## Generalizing to handle both break-points and untruncated mixtures
    # Generate from broken power-law mixture of Truncated Pareto's:
    S <- rbrokenpareto(n=N,x_min=Smin,k=theta,bp=bp,verbose=verbose3)
    idx <- S==Inf
    # sanity check for poor parameter coeeficients
    if(any(idx)){
      S[idx] <- runif(sum(idx),1.5e+268,9e+297)
    }
    # Generate intensity of photon counts:
    lambda.src <- S*E/gamma
    
    # Determine whether observed or missing:
    gp <- g(lambda=lambda.src)
    
    J <- rbinom(N,1,gp)    #length of gp is N
    n <- sum(J)
    
    ################## ################## ###############
    # Store proportion of observed sources
    pi[i] <- n/N
    
    if (verbose) {
      cat(c("Computed pi(theta,Smin,bp) for Smin =",ppaste(signif(Smin,5)),"; bp =",ppaste(signif(bp,5)),"; theta = ",ppaste(round(theta,5))," (",i,"/",length.theta,") ; Pr(Obs) = ",round(pi[i],5),"\n"))            
    }
    
  }
  
  ################## ################## ###############
  # Evaluate proportion of observed sources
  pi.theta.approx <- list("theta"=theta.grid,"pi.theta"=pi)
  
  return(pi.theta.approx)
  
}
