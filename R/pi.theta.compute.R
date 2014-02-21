"pi.theta.compute" <- function(theta.grid, Smin.grid, gamma, g, E,
                               a,b,C,mu,am,bm, length.theta, length.Smin, fixed.Smin, grid.eps=0.0001,
                               fixed.bp=NULL, length.bp=0, bp.grid=NULL,
                               theta.lo=0.1, theta.hi=2.5,
                               Smin.hi=qgamma(p=1.0-grid.eps/5,shape=am,rate=bm),
                               nsamples=100000, mc.integral=TRUE,
                               truncate=TRUE, verbose=FALSE, do.interpolate=FALSE,
                               expand.grid=FALSE, just.return.grid=FALSE){
  
  ####################################################################################################
  # pi.theta.compute   Pre-compute vector of pi(theta,Smin): 4D integral over S,B,L,E 
  #                    Numerical Integration or Monte Carlo for computing probability of average detection/observing a source
  #
  # Input: theta.grid  = NULL or current iteration of theta (column vectors (theta_0,...,theta_m))
  #        Smin   = NULL or minimum flux the sources can be detected to, hyper-parameter in pareto prior for S
  #        gamma  = constant of transformation, energy per photon
  #        g      = function, probability of observing a source
  #        E      = Exposure map
  #        fixed.bp = vector of break-points, fixed for now
  #        fixed.S  = (override) fix the value of S (used mainly for debugging)
  #        truncate = (T/F) crop probability to be between [0,1]
  #        mc.integral = (T/F) method of evaluation of integral: Monte Carlo or Numerical integration
  #        verbose  = (T/F) display progress of program
  #
  # Output: p = whole list["theta","pi"] of unconditional probability of detection, function of theta only
  ####################################################################################################
  
  ###### Construct initial grid
  m <- length(a)
  theta.list <- vector("list",m)
  if (is.null(theta.grid)) {
    
    #### Set-up grid for theta, UN-vectorized, better: storage to list
    for (j in 1:m) {
      theta.grid <- grid.select(a=a[j], b=b[j], grid.eps=grid.eps, grid.length=length.theta, 
                                lowest=theta.lo[j], highest=theta.hi[j])
      theta.list[[j]] <- theta.grid
    }
    
   } else {
    theta.list <- theta.grid # already a list of m vectors
  }
  names(theta.list) <- paste("theta.",1:m,sep="")
  length.theta <- length(theta.list[[1]])
  
  if (all(is.numeric(fixed.bp))) {
    length.bp <- 1
  }
  
  # Begin constructing grid:
  var.grid <- theta.list
  var.length <- list("Smin"=length.Smin, "theta"=length.theta, "bp"=length.bp)
  
  if (is.null(Smin.grid)) {

    #### Set-up grid for Smin
    if (!fixed.Smin){ # Varying Smin case: need a grid over Smin
      Smin.grid <- grid.select(a=am, b=bm, grid.eps=grid.eps, grid.length=length.Smin,
                               highest=Smin.hi)
    } else {
      Smin.grid <- fixed.Smin
    }
  }
  if (all(is.numeric(fixed.bp))) {
    idx <- which(Smin.grid>min(fixed.bp))
    if (any(idx)) {
      warning("Fixed.bp should be larger than Smin when generating grid. Try to reduce the Smin.hi.")
      Smin.grid <- Smin.grid[-idx]
      length.Smin <- var.length$Smin <- length(Smin.grid)
    }
  }
  var.grid$Smin <- Smin.grid
  Smin.list <- list("Smin"=Smin.grid)
  length.Smin <- length(Smin.grid)
  
  if (m>1) {
    # Generate (bp > Smin)
    if (is.null(bp.grid)) {   
      
      #### Set-up grid for bp
      bp.list <- vector("list",0)
      if (any(!fixed.bp)) {
        for (j in 1:(m-1)) { # There are total of m-1 break-points
          # Update each bp subgrid based on bp > smin[k]
          bp.grid <- NULL
          len.k <- ifelse(j==1, length.Smin, length(bp.list[[j-1]]))
          for (k in 1:len.k) {
            low.bound <- ifelse(j==1, Smin.grid[k], bp.list[[j-1]][k])
            lo.bp  <- low.bound*1.001
            hi.bp  <- lo.bp + exp(qnorm(1-grid.eps/100, mean=mu, sd=C))
            tmp1 <- seq(lo.bp, hi.bp, length.out=length.bp-3)
            tmp2 <- exp(qnorm(c(0.15,0.5,0.85),mean=mu,sd=C))
            if (all(tmp2 >= lo.bp)) {
              tmp    <- sort(c(tmp2, tmp1))
            } else { # disregard tmp2 completely
              tmp <- sort(seq(lo.bp, hi.bp, length.out=length.bp)) #similar to tmp1, but longer
            }
            #shape.bp <- am*(1+.0001*10^j)
            #lo.bp  <- pgamma(low.bound*1.0001, shape=am, rate=bm)
            #hi.bp  <- 1-grid.eps/100
            #tmp    <- qgamma(seq(lo.bp, hi.bp, length.out=length.bp), shape=shape.bp, rate=bm)
            bp.grid <- c(bp.grid, tmp)
          }
          bp.list[[j]] <- bp.grid 
          # Replicate entries according to new bp grid
          Smin.grid <- rep(Smin.grid, each=length.bp)
          length.Smin <- length(Smin.grid)
          if (j>1) {
            for (i in 1:(j-1))
            bp.list[[i]] <- rep(bp.list[[i]], each=length.bp)
          }
        }
        length.bp <- length(bp.list[[1]]) 

        # Sanity check
        if (length.Smin != length.bp) {
          stop("There was a problem generating grid for Smin and bp's s.t. Smin[k] < bp.grid. 
               Results should be of same length.")
        }
        
        # Update the storage
        var.grid$Smin <- Smin.grid
        
      } else {
        for (j in 1:(m-1)){
          bp.list[[j]] <- fixed.bp[j] 
        }
      }
      names(bp.list) <- paste0("bp.",1:(m-1))
      
    } else {
      bp.list <- bp.grid 
    }
    
    # Make a matrix of all bp's, Store to be expanded and added AFTER Smin expansion.
    if (length(bp.list[[1]])==1 && length(bp.list)==1) {
      bp.grid <- cbind("bp"=bp.list[[1]])
    } else {
      bp.grid <- sapply(1:(m-1), function(x){cbind(bp.list[[x]])}) # matrix form
      if (class(bp.grid)!="matrix") {
        bp.grid <- matrix(bp.grid,ncol=2)
      }
    }
        
  } else {
    bp.grid <- bp.list <- NULL
  }
    
  if (!expand.grid) {
    if (just.return.grid) {
      ret <- c(list("grid"=NULL), theta.list, Smin.list, bp.list)
      ret$length <- var.length
      return(ret)
    }
    
    stop("Must set 'expand.grid=TRUE' to proceed to evaluate pi.theta.")
    
  } 
  
  ########## Expand grid on theta and Smin ########
  large.grid <- expand.grid(var.grid)
  
  if (m>1) {
    if (all(is.numeric(fixed.bp))) {
      nrep <- nrow(large.grid)
    } else {
      # add expanding grid of bp (must be same length as expanded Smin)
      nrep <- nrow(large.grid)/length.Smin
    }
    large.grid[,(1:(m-1))+m+1] <- apply(bp.grid, 2, rep, each=nrep)
    colnames(large.grid)[(1:(m-1))+m+1] <- paste0("bp.",1:(m-1))
  }
  
  # Update grid variables in matrix/vector form
  large.grid <- data.matrix( large.grid )
  id.th <- grep("theta",colnames(large.grid))
  id.sm <- grep("Smin",colnames(large.grid))
  id.bp <- grep("bp",colnames(large.grid))
  # Retrieve separate grids
  theta.grid <- large.grid[,id.th]
  Smin.grid  <- large.grid[,id.sm]
  length.theta <- length.Smin <- nrow(large.grid)
  if (m>1) {
    bp.grid   <- large.grid[,id.bp]
    length.bp <- length.theta
  } else {
    bp.grid <- length.bp <- NULL
  }
  
  
  if (just.return.grid) {
    ret <- c(list("grid"=large.grid), theta.list, Smin.list, bp.list)
    ret$length <- var.length
    return(ret)
  }
  
  ############################################################
  #### Precompute pi(theta)
  pi <- matrix(0, nrow=length.theta, ncol=1)
  
  if (mc.integral) {
    
    # Use Monte Carlo integration 
    tmp <- pi.theta.eval.mc(nsamples=nsamples, theta.grid=theta.grid, Smin.grid=Smin.grid, bp.grid=bp.grid, 
                            m=m, gamma=gamma, E=E, g=g, verbose=verbose)    
    } # END of integration by MC
  
  pi <- tmp$pi
  
  if (truncate){
    cat('Truncating pi to be in [0,1]...\n')
    pi <- pmin(pi,1)
    pi <- pmax(pi,0) 
  }
  
  # Re-order grid, to be useful for pi.theta.get
  if (m>1) {
    large.grid <- large.grid[,c(id.sm,id.bp[(m-1):1],id.th[m:1])]
  } else {
    large.grid <- large.grid[,c(id.sm,id.th[m:1])]
  }
  # combine results
  if (m>1){
    pi <- c(list("pi"=pi, "grid"=large.grid), theta.list, Smin.list, bp.list)
  } else {
    pi <- c(list("pi"=pi, "grid"=large.grid), theta.list, Smin.list)
  }  
  # add original lengths
  pi$length <- var.length
  
  if (do.interpolate) {
    # Prepare required data structures to use winterpolate()
    # Note: vars must be in SAME order of expand.grid output! 
    tmp <- data.prep(pi=pi$pi, var1=unique(pi$theta[,1]), var2=unique(pi$theta[,2]))
    #tmp <- data.prep(pi=pi$pi, var1=unique(pi$theta[,1]), var2=unique(pi$theta[,2]), var3=unique(pi$Smin))
    x <- tmp$x
    Z <- tmp$Z
    # new.x is a matrix of new values, each column is a variable, each row is a new coordinate for pi
    new.x <- expand.grid("x"=seq(from=range(x$x)[1], to=range(x$x)[2], length.out=500),
                         "y"=seq(from=range(x$y)[1], to=range(x$y)[2], length.out=500))
    #  "u"=seq(from=range(x$u)[1], to=range(x$u)[2], length.out=100))
    # Interpolate
    wintest <- winterpolate(x=x,Z=Z,new.x=new.x,dist.p=0.1,verbose=verbose)
    # Update pi
    pi$pi <- wintest
    pi$theta <- new.x[,1:2]
  }
  
  return(pi)
}

pi.theta.compute <- cmpfun(pi.theta.compute)
