"pi.theta.get" <- function(theta, Smin, bp, gamma, g, E,
                           pi=NULL, sigma=0,
                           mc.integral=TRUE, nsamples=10000,
                           truncate=TRUE, do.interpolate=FALSE, verbose=FALSE,vertex.method=1,method=1){
  
  ####################################################################################################
  # pi.theta.get   Given a single vector of theta & Smin, perform Monte Carlo or Numerical Integration function for 
  #                computing probability of detection/observing a source with broken power-law
  #
  # Input: theta  = current iteration of theta (vector (theta_0,...,theta_m))
  #        Smin   = minimum flux the sources can be detected to, hyper-parameter in pareto prior for S
  #        bp     = vector of break-points
  #        gamma  = constant of transformation, energy per photon
  #        g      = function, probability of observing a source
  #        E      = parameter for exposure time
  #        sigma  = sd parameter of fake approximation to MC error
  #        mc.integral = (T/F) do Monte Carlo approximation (T) or rejection sampling (F)
  #        nsamples = Number of samples for MC approx
  #        truncate = (T/F) crop probability to be between [0,1]
  #        verbose  = (T/F) display progress of program
  #        method   = (1/2) 1=original findInterval function, 2=C-interval findINterval function
  #
  # Output: p = marginal probability of detection, function of theta only
  ####################################################################################################
  
  m <- length(theta)
  
  ### NOTE: This function assumes the following order of the grid in pi:
  # [Smin, bp, theta]
  ### In the precomputed-pi case, access the closest value of pi directly
  if (!is.null(pi)) { 
    # Decide if pi is evaluated over grid of theta alone or combinations of Smin,bp, and theta.
    var.grid       <- pi$grid  # class of grid is matrix
    grid.names     <- colnames(var.grid)
    
    if (vertex.method==1) {
      
      if (verbose) {
        cat("~*~*~*~*~*~*~  Inside pi.theta.get()  ~*~*~*~*~*~*~\n")
      }
      
      # Formula for number of replications of [Smin, bp, theta.2, theta.1] value in partially-reduced grid
      if (m>1) {
        n.reps <- c(pi$length$theta^m * pi$length$bp^((m-1):1), pi$length$theta^(m:1), 1) 
      } else {
        n.reps <- c(pi$length$theta^(m:1), 1) 
      }
      
      # Define coordinate to look-up. NOTE: pi$grid MUST be sorted the same way already.
      new.coord <- matrix(c(Smin, bp, theta),nrow=1)
      nvars     <- ncol(new.coord)
      
      # Try to search the shortened grid first, and then access pi by magic-index
      len.bp <- pi$length$bp
      new.idx <- numeric(nvars)
      new.idx[1] <- findInterval2(new.coord[1], pi$Smin, all.inside=TRUE)
      vec.bp <- pi$bp.1[ (new.idx[1]-1)*len.bp + 1:len.bp ] # sorted portion of bp vector
      new.idx[2] <- findInterval2(new.coord[2], vec.bp, all.inside=TRUE)
      new.idx[3] <- findInterval2(new.coord[3], pi$theta.2, all.inside=TRUE)
      new.idx[4] <- findInterval2(new.coord[4], pi$theta.1, all.inside=TRUE)
      
      ### SUPER FAST INDEX SEARCH, just like magic:
      pi.idx <- sum( (new.idx-1)*n.reps ) + 1
      
      if (verbose) {
        # matched coordinate
        mtc <- c(pi$Smin[new.idx[1]],
                 vec.bp[new.idx[2]], 
                 pi$theta.2[new.idx[3]], 
                 pi$theta.1[new.idx[4]])
        
        #check reality
        idx1 <- pi$grid[,1]==mtc[1]
        idx2 <- pi$grid[,2]==mtc[2]
        idx3 <- pi$grid[,3]==mtc[3]
        idx4 <- pi$grid[,4]==mtc[4]
        pi$grid[idx1 & idx2 & idx3 & idx4]
        pi.idx.true <- which(idx1 & idx2 & idx3 & idx4)
        
        cat("Original coordinate is:\n"); print(new.coord)
        cat("Matched coordinate found  :\n"); print(mtc)
        cat("True matched coordinate is:\n"); print(pi$grid[idx1 & idx2 & idx3 & idx4])
        
        cat("pi.index selected:\n"); print(pi.idx)
        cat("True pi.index is :\n"); print(pi.idx.true)
      }
                  
      # Snap-to-grid idea
      pi.res <- pi$pi[pi.idx]
      
    } else if (vertex.method==2) {
      
      pi.idx         <- 1:nrow(var.grid)
      
      # concept:  rewrite this function using ONLY the call to C function of FinaInterval. 
      # In this case it will avoid checking of the list is sorted, which for is ALWAYS TRUE, 
      # and is the grid of pi is contains NAs. It is not necessary to check this with pi.theta.get.
      # Instead, it will be faster to check it once, before running analyze.mis code. :)
      # f5.i <- function(v)
      #    .Internal(findInterval(v, 0 - .Machine$double.neg.eps, FALSE, FALSE))
      #    .Internal(findInterval(coord[ct], grid[,ct], FALSE, FALSE))
      
      ### Use split-and-find method of searching through the grid
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
        #       if(method==1){
        curr.end.idx1 <- findInterval(coord[ct], grid[,ct]) #, rightmost.closed=TRUE)
        #       } else {
        #         curr.end.idx1 <- .Internal(findInterval(coord[ct], grid[,ct], FALSE, FALSE))
        #       }
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
      if (m>1) {
        n.reps <- c(pi$length$theta^m * pi$length$bp^((m-1):1), pi$length$theta^(m:1), 1) 
      } else {
        n.reps <- c(pi$length$theta^(m:1), 1) 
      }
      
      # Define coordinate to look-up. NOTE: pi$grid MUST be sorted the same way already.
      new.coord <- matrix(c(Smin, bp, theta),nrow=1)
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
      
      if (do.interpolate) {
        # Perform w-interpolation on found interval/voxel
        pi.res <- winterpolate(x=coord.list, Z=pi.voxel, new.x=new.coord, dist.p=0.1, verbose=verbose)
      } else {
        # Snap-to-grid idea
        pi.res <- pi.voxel[1]
      }
      
    } else if (vertex.method==0) {
      
      pi.idx         <- 1:nrow(var.grid)
      
      Smin.grid.idx  <- grep("Smin",grid.names)
      bp.grid.idx    <- grep("bp",grid.names)
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
      for (j in 1:m) {
        # Use binary search to find best entries of theta
        theta.grid.idx   <- grep("theta",colnames(var.grid))[j]
        curr.last.idx <- findInterval(theta[j], var.grid[,theta.grid.idx])
        n.reps        <- pi$length$theta^(m-j) # formula for number of replications of bp value in grid
        if (curr.last.idx==0) { # bin found below minimum value
          curr.idx <- 1:n.reps
        } else {
          curr.idx <- (curr.last.idx - n.reps + 1):curr.last.idx
        }
        # Shorten search grid 
        var.grid      <- var.grid[curr.idx,]
        pi.idx   <- pi.idx[curr.idx]
      }
      
      pi.res <- pi$pi[pi.idx] # Found index of pi
      
    } else {
      # Perform exhaustive search within each variable to select nearby-bin match (unordered variables)
      mm <- ncol(var.grid) # number of all variables (vectors) for grid of pi
      idx <- 1:nrow(var.grid)
      var.look <- c(theta, Smin, bp)
      for (j in 1:mm){
        temp <- abs(var.grid[idx,j] - var.look[j])
        curr.idx <- which(temp == min(temp))
        idx <- idx[curr.idx]
      }
      pi.res <- pi$pi[idx]
      
    }# END vertex method
    
    if (verbose) {
      cat("Selected Interpolated pi:\n"); print(pi.res)
      cat("~*~*~*~*~*~*~  Finished pi.theta.get()  ~*~*~*~*~*~*~\n")
    }
    return(pi.res)
    
  } # END !is.null(pi)
  
  
  # Evaluate constant pi with error
  #   pi <- g(lambda=Smin*E/gamma) + rnorm(1,mean=0,sd=sigma)  # NOTE: this g==CONST
  
  # Use Monte Carlo integration 
  tmp <- pi.theta.eval.mc(nsamples=nsamples, theta.grid=theta, Smin.grid=Smin, bp.grid=bp, 
                          m=m, gamma=gamma, E=E, g=g, verbose=verbose)
  
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
