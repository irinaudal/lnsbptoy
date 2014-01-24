"simulate.mix" <- function(a, b, C=NULL, mu=NULL, am=NULL, bm=NULL, # vectors of length corresponding to the number of Pareto's 
                           alpha, beta, gamma, E,
                           fixed.Smin=FALSE, fixed.bp=NULL, fixed.theta=FALSE, fixed.N=FALSE,
                           g=function(lambda){
                             return(rep(0.8,length(lambda)))
                           }, 
                           verbose=FALSE ){
  
  ####################################################################################################
  # simulate.mix   Function simulates data: observed total(bkg+src) photon counts for LogN-LogS project
  #
  # Input: a     = shape hyper-parameter(s) in gamma prior for theta (single) or (theta_0,...,theta_m) (broken power-law)
  #        b     = rate hyper-parameter(s) in gamma prior for theta (single) or (theta_0,...,theta_m) (broken power-law)
  #        d     = vector if hyper-parameter(s) in dirichlet prior for mixture proportions (p_0,...,p_m)) (broken power-law)
  #        am    = shape hyper-parameter(s) in gamma prior for Smin
  #        bm    = rate hyper-parameter(s) in gamma prior for Smin
  #        alpha = target number of sccessful trials in negbinom prior for N (scalar)
  #        beta  = dispersion parameter in negbinom prior for N (scalar)
  #        mu    = mean parameter in normal proposal for break-points
  #        C     = varianca-covariance matrix of parameters in normal proposal for break-points
  #        gamma = constant of transformation, energy per photon
  #        pble  = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        fixed.N     = Value/FALSE: allow for user specified theta - no mcmc
  #        fixed.theta = Value/FALSE: allow for user specified theta - no mcmc
  #        fixed.Smin  = Value/FALSE: allow for user specified theta - no mcmc
  #        fixed.p     = Value/FALSE, mixture Pareto probabilities: allow for user specified theta - no mcmc
  #        fixed.bp    = either NULL (no fixed break-points), or, a vector of length p specifying the break-points
  #        g.type   = type of g-function {"step","smooth","table"}
  #        g        = function, probability of observing a source
  #        area.conv.func = function, convert from off-axis angle to area [pix]
  #        verbose  = (T/F) display progress of program
  #
  # Output: result  = list object with values: obs = total photon counts, (Ytot, n, E.obs, L.obs)
  #                                            par = vector of parameters (N, theta, S.obs, Smin)
  #                                            mis = list of missing photon counts (N.mis,S.mis,Y.obs.src, 
  #                                                   Y.mis.src,Y.mis.tot,Y.obs.bkg,Y.mis.bkg)
  #                                            mix = list of mixture parameters (I.idx, p)
  #                                            pble = pble object used to generate data
  ####################################################################################################
  
  if (verbose==TRUE) {
    cat("Generating data.. \n") 
  }
  
  m <- length(a)
  if (length(b) != m){
    stop("'a' and 'b' must be the same length")
  }
  
  # Generate the total number of sources:
  if(!is.logical(fixed.N)){
    N <- fixed.N
  } else {
    N <- rnbinom(1,size=alpha,prob=beta/(1+beta))
  }
  if(verbose){
    cat("N =\n"); print(N)
  }
  
  # Generate theta(s):
  idx <- (fixed.theta==TRUE | fixed.theta==FALSE)
  if (any(idx)){     #some/all theta values are not fixed and need to be generated
    theta <- rgamma(n=m,shape=a,rate=b)
    idx2 <- !idx   # index of fixed values
    theta[idx2] <- fixed.theta[idx2]  # overwrite with fixed values
  } else {   #all values for theta are already fixed/given
    theta <- fixed.theta
  }
  if(verbose){
    cat("theta =\n"); print(theta)
  }
  
  
  # Generate Smin, the minimum source flux, informative prior:
  idx <- is.logical(fixed.Smin)
  if (idx){     #some/all Smin values are not fixed and need to be generated
    Smin <- rgamma(n=1,shape=am,rate=bm)
    idx2 <- !idx   # index of fixed values
    Smin[idx2] <- fixed.Smin[idx2]  # overwrite with fixed values
  } else {   #all values for Smin are already fixed/given
    Smin  <- fixed.Smin
  }
  if(verbose){
    cat("Smin =\n"); print(Smin)
  }
  
  
  # Generate the breakpoints:
  # fixed.bp has 3 options:
  # Option 2 :: fixed.bp == (vector of scalars)  ==> Break-points, fixed to user-specified values
  # Option 3 :: fixed.bp == FALSE or TRUE        ==> Break-points, must be estimated
  
  if (any(is.logical(fixed.bp))){
    #       if (any(!fixed.bp)){  
    # Option 3 :: Break-points, must be estimated
    bp <- cumsum(exp(rnorm(n=m-1, mean=mu, sd=C))) + Smin
  } else {
    # Option 2 :: Break-points, fixed to user-specified values
    # Must be of length (m-1):
    if (length(fixed.bp) != (m-1)){
      stop("When 'fixed.bp' is a vector of scalars, it must be of length 'length(a)-1' and 'length(b)-1'")
    }
    bp <- fixed.bp
  } # END generating bp
  if(verbose){
    cat("bp =\n"); print(bp)
  }

  # Generate source flux S:
  ## Generalizing to handle both break-points and untruncated mixtures
  # Generate from broken power-law mixture of Truncated Pareto's:
  S     <- rbrokenpareto(n=N,x_min=Smin,k=theta,bp=bp,verbose=verbose)
  
  # Generate photon counts:
  Ybkg <- 0
  lambda.src <- S*E/gamma    #CDFS: src_count [ct] = (flux [erg/s/cm^2] * exposure_time [s])/ (eff_mean_src_exp * 1.06e-11 [ergs/cm^2/ct])
  Ysrc <- rpois(N,lambda.src)
  
  # sanity check
  if(any(is.na(Ysrc))) {
    id <- is.na(Ysrc)
    Ysrc[id] <- lambda.src[id]
  }
  # Manual hack for unusual large sources, when model is wrong.
  if(any(Ysrc==Inf)){
    id <- Ysrc==Inf
    Ysrc[id] <- 1e299
    lambda.src[id] <- 1e299
  }
  
  Ytot <- Ysrc + Ybkg
  
  # Determine whether observed or missing:
  gp <- g(lambda=lambda.src)
  
  J <- rbinom(N,1,gp)    #length of gp is N
  n <- sum(J)
  
  if(verbose){
    cat("gp in simulate.mix.R:  \n"); print(gp)
    cat("J:  \n"); print(J)
    cat("S:  \n"); print(S)
    cat("E:  \n"); print(E)
    cat("lambda.src:  \n"); print(lambda.src)
  }  
  if (verbose) {
    cat("Storing result.. \n") 
  }   
  
  par <- list("N"=N,"theta."=theta,"tau."=c(Smin,bp),"S.obs."=S[J==1])

  result <- list("obs"=list("Y.obs.tot"=Ytot[J==1],"n"=n),
                 "par"=par, 
                 "mis"=list("N.mis"=sum(1-J),"S.mis"=S[J==0],"Y.mis.tot"=Ytot[J==0]),
                 "orig"=list("S"=S,"E"=E,"gp"=gp,"lambda.src"=lambda.src,"Ybkg"=Ybkg,"Ysrc"=Ysrc))
  return(result)
}

simulate.mix <- cmpfun(simulate.mix)
