"simulate.mix" <- function(a, b, d=NULL, C=NULL, mu=NULL, am=NULL, bm=NULL, # vectors of length corresponding to the number of Pareto's 
                           alpha, beta, gamma, pble,
                           fixed.bp=NULL, fixed.theta=FALSE, fixed.N=FALSE, fixed.Smin=FALSE, fixed.p=NULL,
                           g.type="step", g=function(lambda,bg,E,L,g.type){
                             g.compute(lambda=lambda,bg=bg,E=E,L=L,g.type=g.type)
                           }, area.conv.func=function(off){10.57508 * exp(0.5232*off)},
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
  
  # Create pble object, if needed
  if (class(pble)=="character") {
    pble <- new.pble.table(effects.file=pble)
  } else if (class(pble)!="pble.basic" && class(pble)!="pble.table"){
    stop("'pble' could not be created.")
  }  
  
  m <- length(a)
  if (length(b) != m){
    stop("'a' and 'b' must be the same length")
  }
  use.bp <- !is.null(fixed.bp)
  use.mix <- !is.null(fixed.p)
  # sanity check
  if(use.bp & use.mix){
    stop("Check parameter specification! Cannot be both bp and mix case.")
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
  # Option 1 :: fixed.bp == NULL                 ==> No break-points
  # Option 2 :: fixed.bp == (vector of scalars)  ==> Break-points, fixed to user-specified values
  # Option 3 :: fixed.bp == FALSE or TRUE        ==> Break-points, must be estimated
  
  if (!use.bp){
    # Option 1 :: No break-points
    bp <- NULL
  } else {
    if (any(is.logical(fixed.bp))){
      #       if (any(!fixed.bp)){  
      # Option 3 :: Break-points, must be estimated
      if (verbose) {
        cat("Generating K from prior.\n")
      }
      if (length(C)==m) {
        # Generate (Smin.bp) from prior, overwrite previous version of Smin
        K <- cumsum(exp(rnorm(n=m, mean=mu, sd=C)))
        Smin <- K[1]
        bp <- K[-1]
      } else if (length(C)==m-1) {
        bp <- cumsum(exp(rnorm(n=m-1, mean=mu, sd=C))) + Smin
        
      }      
      #       } 
      #       else {
      #         stop("When specifying a break-point, you must use 'fixed.bp=value_you_want_to_fix_to'")
      #         # Use the user-specified break-points
      #       }
    } else {
      # Option 2 :: Break-points, fixed to user-specified values
      # Must be of length (m-1):
      if (length(fixed.bp) != (m-1)){
        stop("When 'fixed.bp' is a vector of scalars, it must be of length 'length(a)-1' and 'length(b)-1'")
      }
      bp <- fixed.bp
    }
  } # END generating bp
  if(verbose){
    cat("bp =\n"); print(bp)
  }
  
  
  # Mixture-pareto has 2 options:
  # Option 1 :: d == NULL                        ==> No mixtures
  # Option 2 :: d == (vector of scalars)         ==> Mixtures, user-specified
  
  if (!use.mix){
    # Option 1 :: No mixtures
    p <- 1
  } else {  
    # Option 2 :: Mixtures, total number fixed to user-specified values
    if (length(d) != (m)){
      stop("When specifying mixtures, 'd' is a vector of scalars, it must be of length 'length(a)' and 'length(b)'")
    }
    
    # Generate mixture proportion(s):
    idx <- (fixed.p==TRUE | fixed.p==FALSE) #index of non-fixed-at-value entries
    if (all(idx)){
      p <- rdirichlet(n=1,alpha=d)
    } else { #some/all p values are fixed
      ## Suppose (X1,X2,X3,X4,X5)~Dirichlet(d1,d2,d3,d4,d5), s.t. sum(X_i)=1. 
      ## Then, (X1,X2,X3,X4+X5)~Dirichlet(d1,d2,d3,d4+d5), s.t. sum(X_i)=1.
      ## And, (X1,X2,X3)*1/(1-(X4+X5)) ~ Dirichlet(d1,d2,d3)  indep of  Y=(X4+X5) ~ Beta(d4+d5,d1+d2+d3) 
      ## So, X1,X2,X3 | X4, X5  ~  (1-(X4+X5))*Z, where Z~Dirichlet(d1,d2,d3)  
      p <- fixed.p
      p.tilde <- rdirichlet(n=1,alpha=d[idx])
      p[idx] <- p.tilde*(1-sum(p[!idx]))
    } 
  } # END generating mixture probabilities  
  
  # Generate source flux S:
  ## Generalizing to handle both break-points and untruncated mixtures
  if (use.bp){
    # Generate from broken power-law mixture of Truncated Pareto's:
    S     <- rbrokenpareto(n=N,x_min=Smin,k=theta,bp=bp,verbose=verbose)
    I.idx <- NULL
  } else if (use.mix){
    # Generate from mixture of Pareto's:
    tmp   <- rmixpareto(n=N,p=p,k=theta,x_min=Smin) #this function is already compatible with single Pareto case
    I.idx <- tmp$I.idx
    S     <- tmp$S
  } else {
    # Generate from single Pareto:
    S     <- rpareto(n=N,k=theta,x_min=Smin)
    I.idx <- NULL
  }
  
  # Generate standard background, effective area etc.:  
  par <- sample.pble(pble,N)
  bg  <- par$B
  L   <- par$L
  E   <- par$E  
  pix.area <- area.conv.func(off=L)
  k    <- bg*pix.area   
  if(verbose){
    cat("bg =\n"); print(bg)
    cat("L =\n"); print(L)
    cat("E =\n"); print(E)
    cat("pix.area =\n"); print(pix.area)
    cat("k =\n"); print(k)
  }
  
  
  # Generate photon counts:
  Ybkg <- rpois(N,k)
  lambda.src <- S*E/gamma    #CDFS: src_count [ct] = (flux [erg/s/cm^2] * exposure_time [s])/ (eff_mean_src_exp * 1.06e-11 [ergs/cm^2/ct])
  #lambda.src <- min(lambda.src, 1e299)
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
  gp <- g(lambda=lambda.src,bg=bg,E=E,L=L,g.type=g.type,verbose=verbose)
  
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
  
  if (use.mix){
    par <- list("N"=N,"theta."=theta,"Smin"=Smin,"p."=p,"S.obs."=S[J==1])    
  } else if(use.bp) {
    if (m-1==1) {
      par <- list("N"=N,"theta."=theta,"Smin"=Smin,"bp"=bp,"S.obs."=S[J==1])
    } else {
      par <- list("N"=N,"theta."=theta,"Smin"=Smin,"bp."=bp,"S.obs."=S[J==1])
    }
  } else {
    par <- list("N"=N,"theta"=theta,"Smin"=Smin,"S.obs."=S[J==1])
  }
  
  result <- list("obs"=list("Y.obs.tot"=Ytot[J==1],"n"=n,"E.obs"=E[J==1],"bg.obs"=bg[J==1],"L.obs"=L[J==1],"k.obs"=k[J==1]),
                 "par"=par, 
                 "mis"=list("N.mis"=sum(1-J),"S.mis"=S[J==0],"Y.obs.src"=Ysrc[J==1], 
                            "Y.mis.src"=Ysrc[J==0],"Y.mis.tot"=Ytot[J==0], 
                            "Y.obs.bkg"=Ybkg[J==1],"Y.mis.bkg"=Ybkg[J==0]),
                 "mix"=list("I.idx"=I.idx),
                 "pble"=pble,
                 "orig"=list("S"=S,"k"=k,"bg"=bg,"gp"=gp,"lambda.src"=lambda.src,"Ybkg"=Ybkg,"Ysrc"=Ysrc))
  return(result)
}

simulate.mix <- cmpfun(simulate.mix)
