"pi.theta.get" <- function(theta, Smin, gamma, g.type, g, pble, bp=NULL, p.t=NULL, pi=NULL,
                           delta=10^(-3), # can be zero...
                           length.S=1000, fixed.S=NULL, use.bp=FALSE, use.mix=FALSE, debug=FALSE, nsamples=10000,
                           truncate=TRUE, mc.integral=TRUE, verbose=FALSE, print.time=FALSE, 
                           new.method=2){
  
  ####################################################################################################
  # pi.theta.get   Given a single vector of theta & Smin, perform Monte Carlo or Numerical Integration function for 
  #                computing probability of detection/observing a source with broken power-law
  #                OR
  #                Given object pi[$pi, $theta, $p.t, $Smin], access required probability
  #
  # Input: theta  = current iteration of theta (vector (theta_0,...,theta_m))
  #        Smin   = minimum flux the sources can be detected to, hyper-parameter in pareto prior for S
  #        gamma  = constant of transformation, energy per photon
  #        g.type = type of g-function {"step","smooth","table"}
  #        g      = function, probability of observing a source
  #        pble   = B,L,E parameterss /or joint distribution of B,L,E (from file)
  #        bp     = vector of break-points
  #        p.t    = previous vector of mixture component probabilities
  #        pi     = NULL/R.object: a list("pi","theta") with detection probalities: pi(theta)
  #        delta    = small quantle error for irregular grid estimation 
  #        length.S = number of points used in the numerical integration for S
  #        fixed.S  = (override) fix the value of S (used mainly for debugging)
  #        use.bp   = (T/F) break-point pareto version
  #        use.mix  = (T/F) mixture pareto version
  #        truncate = (T/F) crop probability to be between [0,1]
  #        debug    = (T/F) debug numerical integral
  #        verbose  = (T/F) display progress of program
  #
  # Output: p = marginal probability of detection, function of theta only
  ####################################################################################################
  
  return(0.9*rep(1,length(Smin)) + rnorm(n=length(Smin), mean=0, sd=0.00001))
  
  ### NOTE: This function assumes the following order of the grid in pi:
  # [Smin, bp, p.t, theta]
  
  # In the precomputed-pi case, access the closest value of pi directly
  if (!is.null(pi)) { 
    # Decide if pi is evaluated over grid of theta alone or combinations of Smin,bp,p.t, and theta.
    var.grid       <- pi$grid  # class of grid is matrix
    grid.names     <- colnames(var.grid)
    m              <- length(theta)
    pi.idx         <- 1:nrow(var.grid)
    
    if (new.method==1) {
      
      Smin.grid.idx  <- grep("Smin",grid.names)
      bp.grid.idx    <- grep("bp",grid.names)
      p.t.grid.idx   <- grep("p.t",grid.names)
      theta.grid.idx <- grep("theta",grid.names)
      
      # Add dummy variables for replicates calculations
      if (length(bp.grid.idx)==0) {
        pi$length$bp <- 1
      }
      
      if (length(Smin.grid.idx)>0) {
        # Use binary search to find best entries of Smin
        curr.last.idx <- findInterval(Smin, var.grid[,Smin.grid.idx])
        n.reps        <- pi$length$theta^m * pi$length$bp # formula for number of replications of Smin value in grid
        curr.idx      <- (curr.last.idx - n.reps + 1):curr.last.idx
        # Shorten search grid 
        var.grid      <- var.grid[curr.idx,]
        pi.idx   <- pi.idx[curr.idx]
      }
      if (length(bp.grid.idx)>0) { # TODO: change this to multivariate case for general m, instead of univariate.
        # Use binary search to find best entries of bp
        bp.grid.idx   <- grep("bp",colnames(var.grid))
        curr.last.idx <- findInterval(bp, var.grid[,bp.grid.idx])
        n.reps        <- pi$length$theta^m # formula for number of replications of bp value in grid
        curr.idx      <- (curr.last.idx - n.reps + 1):curr.last.idx
        # Shorten search grid 
        var.grid      <- var.grid[curr.idx,]
        pi.idx   <- pi.idx[curr.idx]
      }
      if (length(p.t.grid.idx)>0) { # TODO: change this to multivariate case for general m, instead of univariate.
        # Use binary search to find best entries of p.t
        p.t.grid.idx  <- grep("p.t",colnames(var.grid))
        curr.last.idx <- findInterval(p.t, var.grid[,p.t.grid.idx])
        n.reps        <- pi$length$theta^m # formula for number of replications of p.t value in grid
        curr.idx      <- (curr.last.idx - n.reps + 1):curr.last.idx
        # Shorten search grid 
        var.grid      <- var.grid[curr.idx,]
        pi.idx   <- pi.idx[curr.idx]
      }
      for (j in m:1) {
        # Use binary search to find best entries of theta
        theta.grid.idx   <- grep("theta",colnames(var.grid))[j]
        curr.last.idx <- findInterval(theta[j], var.grid[,theta.grid.idx])
        n.reps        <- pi$length$theta^(j-1) # formula for number of replications of bp value in grid
        curr.idx      <- (curr.last.idx - n.reps + 1):curr.last.idx
        # Shorten search grid 
        var.grid      <- var.grid[curr.idx,]
        pi.idx   <- pi.idx[curr.idx]
      }
      
      return( pi$pi[ pi.idx ] ) # Found index of pi
      
    } else if (new.method==2) {
       
      "splitNfind" <- function (grid, pi.idx, coord, reps, ct=0, verb=verbose) {
        ################################################
        # Pseudocode for Recursive Function for finding the correct index of all voxel boundaries 
        # around the new.coord p-dim variable. Binary search is used once per execution of function.
        # function split-and-find-index:
        # inputs: pi.idx=vetor of indices to look in=nrows of var.grid.
        #         var.grid=(n x p) matrix of coordinates
        #
        # 1. if length(pi.idx)==1, then return list of pi.idx value and var.grid vector; else
        #
        # 2. find best match coordinate
        #    split var.grid in two
        #    (a) concatinate result from CALL to split-and-find-index(var.grid.1, pi.idx.1)
        #    (b) concatinate result from CALL to split-and-find-index(var.grid.2, pi.idx.2)
        #    return combined result as a list
        #
        # end split-and-find-index
        ################################################
        
        # Check the base case
        if (length(pi.idx)==1) {
          return( list("corner"=grid, "idx"=pi.idx) )
        } 
        # Increment counter variable
        ct <- ct+1 
        # Use binary search to find best match bin-location entries of current variable
        n.rep <- reps[ct]
        curr.end.idx1 <- findInterval(coord[ct], grid[,ct]) #, rightmost.closed=TRUE)
        if (curr.end.idx1==0) { # bin found below minimum value
          curr.idx1 <- 1:n.rep
        } else {
          curr.idx1 <- (curr.end.idx1 - n.rep + 1):curr.end.idx1
        }
        # Save index of upper bin, as well, for voxel. No binary search needed
        curr.end.idx2 <- curr.end.idx1 + 1
        if (curr.end.idx2>nrow(grid)) { # bin found above minimum value
          curr.idx2 <- curr.idx1
        } else {
          curr.idx2 <- curr.end.idx2:(curr.end.idx1 + n.rep)
        }
        if (verb) {
          cat("curr.end.idx1:\n"); print(curr.end.idx1)
          cat("curr.end.idx2:\n"); print(curr.end.idx2)
          cat("curr.idx1:\n"); print(curr.idx1)
          cat("curr.idx2:\n"); print(curr.idx2)
          cat("nrow(grid):\n"); print(nrow(grid))
          cat("coordinate:\n"); print(coord[ct])
        }
        
        # Split and Shorten search grid 
        grid1   <- grid[curr.idx1,]
        pi.idx1 <- pi.idx[curr.idx1]
        grid2   <- grid[curr.idx2,]
        pi.idx2 <- pi.idx[curr.idx2]
        
        if (verb) {
          cat("pi.idx1:\n"); print(pi.idx1)
          cat("pi.idx2:\n"); print(pi.idx2)
        }
        
        # Call split-and-find-index function recursively on reduced grids
        res1 <- splitNfind(grid=grid1, pi.idx=pi.idx1, coord=coord, reps=reps, ct=ct)
        res2 <- splitNfind(grid=grid2, pi.idx=pi.idx2, coord=coord, reps=reps, ct=ct)
        # Combine and return results
        corner.coord <- rbind(res1$corner, res2$corner)
        pi.found.idx <- c(res1$idx, res2$idx)
        return( list("corner"=corner.coord, "idx"=pi.found.idx) )
        
      } # END splitNfind recursive function
      
      if (verbose) {
        cat("~*~*~*~*~*~*~  Inside pi.theta.get()  ~*~*~*~*~*~*~\n")
      }
      
      # Formula for number of replications of [Smin, bp, theta.2, theta.1] value in partially-reduced grid
      if (use.bp) {
        n.reps <- c(pi$length$theta^m * pi$length$bp^((m-1):1), pi$length$theta^(m:1), 1) 
      } else if (use.mix) {
        n.reps <- c(pi$length$theta^m * pi$length$p.t^((m-1):1), pi$length$theta^(m:1), 1) 
      } else {
        n.reps <- c(pi$length$theta^(m:1), 1) 
      }
      
      # Define coordinate to look-up. NOTE: pi$grid MUST be sorted the same way already.
      new.coord <- matrix(c(Smin, bp, p.t, theta),nrow=1)
      nvars     <- ncol(new.coord)
      
      # Find voxel using recursion and binary search
      voxel <- splitNfind(grid=var.grid, pi.idx=pi.idx, coord=new.coord, reps=n.reps)
      
      # Store voxel information for interpolation within a voxel of pi, stored as array
      pi.voxel   <- array(pi$pi[voxel$idx],rep(2,nvars))
      unique.coord <- voxel$corner[c(1,2^nvars),]
      # Check for errors
      zero.dim <- sapply( apply(unique.coord,2,unique), length)==1
      if (any(zero.dim)) {
        return (mean(pi.voxel))
      }
      coord.list <- split(unique.coord, col(unique.coord))
      
      if (verbose) {
        cat("Voxel boundaries and value of pi:\n"); print(cbind(voxel$corners, c(pi.voxel)))
      }
      
      # Perform w-interpolation on found interval/voxel
      pi.res <- winterpolate(x=coord.list, Z=pi.voxel, new.x=new.coord, dist.p=0.1, verbose=verbose)
      
      if (verbose) {
        cat("Selected Interpolated pi:\n"); print(pi.res)
        cat("~*~*~*~*~*~*~  Finished pi.theta.get()  ~*~*~*~*~*~*~\n")
      }
      return(pi.res)
      
    } else if (new.method==FALSE) { 
      # Perform exhaustive search within each variable to select nearby-bin match (unordered variables)
      mm <- ncol(var.grid) # number of all variables (vectors) for grid of pi
      idx <- 1:nrow(var.grid)
      var.look <- c(theta, Smin, bp, p.t)
      for (j in 1:mm){
        temp <- abs(var.grid[idx,j] - var.look[j])
        curr.idx <- which(temp == min(temp))
        idx <- idx[curr.idx]
      }
      
      return( pi$pi[ idx ] )
      
    } 
        
  } # END !is.null(pi)
  
  
  m <- length(theta)
  
  if (mc.integral) {
    
    # Use Monte Carlo integration 
    tmp <- pi.theta.eval.mc(nsamples=nsamples, theta.grid=theta, Smin.grid=Smin, bp.grid=bp, p.t.grid=p.t, 
                            use.bp=use.bp, use.mix=use.mix,
                            m=m, gamma=gamma, pble=pble, g=g, g.type=g.type, verbose=verbose)    
    
  } else {
    
    # Use numerical integration:
    tmp <- pi.theta.eval.numint(theta.grid=theta, Smin.grid=Smin, bp.grid=bp, p.t.grid=p.t,  
                                gamma=gamma, pble=pble, g=g, g.type=g.type, delta=delta,
                                m=m, fixed.S=fixed.S, length.S=length.S,
                                use.bp=use.bp, use.mix=use.mix, debug=debug, verbose=verbose)    
  }
  
  pi <- tmp$pi
  
  if (truncate) {
    if(verbose)
      cat('Truncating pi to be in [0,1]...\n')
    pi <- pmin(pi,1)
    pi <- pmax(pi,0) 
  }
  
  return(pi)
  
}

pi.theta.get <- cmpfun(pi.theta.get)
