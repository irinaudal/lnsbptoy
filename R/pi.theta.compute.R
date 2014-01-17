"pi.theta.compute" <- function(theta.grid, Smin.grid, gamma, g.type, g, pble,
                               a,b,C,mu,am,bm, length.theta, length.Smin, fixed.Smin, grid.eps=0.0001,
                               fixed.bp=NULL, length.bp=0, bp.grid=NULL, p.t=NULL, 
                               theta.lo=0.1, theta.hi=2.5,
                               Smin.hi=qgamma(p=1.0-grid.eps/5,shape=am,rate=bm),
                               delta=10^(-3), # can be zero...
                               nsamples=1000000, length.S=1000, fixed.S=NULL, use.bp=FALSE, use.mix=FALSE,
                               debug=FALSE, truncate=TRUE, mc.integral=TRUE, interpolate.pi=FALSE, verbose=FALSE, 
                               just.return.grid=FALSE){
  
  ####################################################################################################
  # pi.theta.compute   Pre-compute vector of pi(theta,Smin): 4D integral over S,B,L,E 
  #                    Numerical Integration or Monte Carlo for computing probability of average detection/observing a source
  #
  # Input: theta.grid  = NULL or current iteration of theta (column vectors (theta_0,...,theta_m))
  #        Smin   = NULL or minimum flux the sources can be detected to, hyper-parameter in pareto prior for S
  #        gamma  = constant of transformation, energy per photon
  #        g.type = type of g-function {"step","smooth","table"}
  #        g      = function, probability of observing a source
  #        pble   = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        p.t    = previous vector of mixture component probabilities
  #        fixed.bp = vector of break-points, fixed for now
  #        delta    = small quantle error for irregular grid estimation 
  #        length.S = number of points used in the numerical integration for S
  #        fixed.S  = (override) fix the value of S (used mainly for debugging)
  #        use.bp   = (T/F) break-point pareto version
  #        use.mix  = (T/F) mixture pareto version
  #        debug    = (T/F) used for debugging purposes
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
                                lowest=theta.lo, highest=theta.hi)
      theta.list[[j]] <- theta.grid
    }
    
    #UPDATE grid of theta for individual
    #theta.tmp <- theta.grid[50,2]
    #theta.grid[,2] <- theta.tmp #update theta2 to be identical everywhere.
    #ii <- sim_num-sim_start
    #theta.grid <- theta.grid[1:2 + 2*(ii-1),]
  } else {
    theta.list <- theta.grid # already a list of m vectors
  }
  names(theta.list) <- paste("theta.",1:m,sep="")
  length.theta <- length(theta.list[[1]])
  
  # Begin comstructing grid:
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
  var.grid$Smin <- Smin.grid
  Smin.list <- list("Smin"=Smin.grid)
  length.Smin <- length(Smin.grid)
  
  if (use.bp) {
    # Generate (bp > Smin)
    if (is.null(bp.grid)) {    
      #### Set-up grid for bp
      if (any(!fixed.bp)) {
        bp.list <- vector("list",0)
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
        names(bp.list) <- paste("bp.",1:(m-1),sep="")
        var.grid$Smin <- Smin.grid
        
      } else {
        stop("Code has not been implemented yet to specify fixed version of bp.grid. Don't know if it is a list or vector or matrix?")
        bp.grid <- fixed.bp
        length.bp <- nrow(bp.grid)
      }
    } else {
      bp.list <- bp.grid 
    }
    
    # DO NOT add bp to the list yet
    #var.grid <- c(var.grid, bp.list)
   
    # Make a matrix of all bp's, Store to be expanded and added AFTER Smin expansion.
    if (length(bp.list[[1]])==1 && length(bp.list)==1) {
      bp.grid <- cbind("bp"=bp.list[[1]])
    } else {
      bp.grid <- sapply(1:(m-1), function(x){cbind(bp.list[[x]])}) # matrix form
    }
        
  } else {
    bp.grid <- bp.list <- NULL
  }
  
  if (use.mix){
    # TODO: check that this part of code works well with LIST version of var.grid and pi
    #### Set-up grid for p.t (the mixture proportion)
    p.t.grid <- matrix(seq(0,1,by=0.1),ncol=1)
    p.t.list <- list("p.t.1"=p.t.grid) #, "p.t.2"=1-p.t.grid)
    var.grid <- c(var.grid, p.t.list)
  } else {
    p.t.grid <- p.t.list <- NULL
  }
    
  ########## Expand grid on theta and Smin (also: bp, or p.t) ########
  large.grid <- expand.grid(var.grid)
    
  if (use.bp) {
    # add expanding grid of bp (must be same length as expanded Smin)
    nrep <- nrow(large.grid)/length.Smin
    large.grid$bp <- apply(bp.grid, 2, rep, each=nrep)
  } else if (use.mix) {
    # add expanding grid of 1-p.t (must be same length as expanded p.t)
    nrep <- nrow(large.grid)/length.p.t
    large.grid$p.t <- apply(1-p.t.grid, 2, rep, each=nrep)
  }
    
  # Update grid variables in matrix/vector form
  large.grid <- data.matrix( large.grid )
  id.th <- grep("theta",colnames(large.grid))
  id.sm <- grep("Smin",colnames(large.grid))
  id.bp <- grep("bp",colnames(large.grid))
  id.p  <- grep("p.",colnames(large.grid))
  # Retrieve separate grids
  theta.grid <- large.grid[,id.th]
  Smin.grid  <- large.grid[,id.sm]
  length.theta <- length.Smin <- nrow(large.grid)
  if (use.bp) {
    bp.grid   <- large.grid[,id.bp]
    length.bp <- length.theta
  } else if (use.mix) {
    p.t.grid   <- large.grid[,id.p]
    length.p.t <- length.theta
  } else {
    p.t.grid <- bp.grid <- length.bp <- length.p.t <- NULL
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
    tmp <- pi.theta.eval.mc(nsamples=nsamples, theta.grid=theta.grid, Smin.grid=Smin.grid, bp.grid=bp.grid, p.t.grid=p.t.grid, 
                            use.bp=use.bp, use.mix=use.mix,
                            m=m, gamma=gamma, pble=pble, g=g, g.type=g.type, verbose=verbose)    
  } else {
    
    # Use numerical integration:  
    tmp <- pi.theta.eval.numint(theta.grid=theta.grid, Smin.grid=Smin.grid, bp.grid=bp.grid, p.t.grid=p.t.grid, 
                                gamma=gamma, pble=pble, g=g, g.type=g.type, delta=delta,
                                m=m, fixed.S=fixed.S, length.S=length.S,
                                use.bp=use.bp, use.mix=use.mix, debug=debug, verbose=verbose)    
  } # END of integration by MC or NI
  
  pi <- tmp$pi
  
  if (truncate){
    cat('Truncating pi to be in [0,1]...\n')
    pi <- pmin(pi,1)
    pi <- pmax(pi,0) 
  }
  
  # Re-order grid, to be useful for pi.theta.get
  if (use.bp) {
    large.grid <- large.grid[,c(id.sm,id.bp[(m-1):1],id.th[m:1])]
  } else if (use.mix) {
    large.grid <- large.grid[,c(id.sm,id.p[(m-1):1],id.th[m:1])]
  } else {
    large.grid <- large.grid[,c(id.sm,id.th[m:1])]
  }
  # combine results
  if (use.bp){
    pi <- c(list("pi"=pi, "grid"=large.grid), theta.list, Smin.list, bp.list)
  } else if (use.mix) {
    pi <- c(list("pi"=pi, "grid"=large.grid), theta.list, Smin.list, p.t.list)
  } else {
    pi <- c(list("pi"=pi, "grid"=large.grid), theta.list, Smin.list)
  }  
  # add original lengths
  pi$length <- var.length
  if (debug) { # only for numerical-integration case
    pi$grid.frame    <- tmp$grid.frame
    pi$db.grid.theta <- tmp$db.grid.theta
  }
  
  if (interpolate.pi) {
    
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
    ##pi$Smin  <- new.x[,1]
    ##pi$theta <- new.x[,2:3]
    
    #     library(akima)
    #     wintest <- interp(x=pi$theta[,1],y=pi$theta[,2], z=Z, xo=new.x[,1], yo=new.x[,2], linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL, ncp = NULL)
    #     pi$pi <- winterp$z
    
    #  stop("HERE")
    
    #     interp2d <- function(old, newx, newy) {
    #       interp.surface.grid(list(x=seq(nrow(old)),y=seq(ncol(old)),z=old),
    #                           list(x=seq(1,nrow(old),length=newx),
    #                                y=seq(1,ncol(old),length=newy)))$z
    #     }
    #     
    #     library(fields)
    #     wintest <- interp2d(old=Z, newx=100, newy=100)
    
  } # END inverse-weight-distance interpolation
  
  return(pi)
}

pi.theta.compute <- cmpfun(pi.theta.compute)
