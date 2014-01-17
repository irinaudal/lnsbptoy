"pi.theta.eval.numint" <- function(theta.grid=matrix(seq(0.1,2.1,by=0.01),ncol=1), Smin.grid=10^-17, bp.grid=10^-16, p.t.grid=1,
                                   m=1, gamma, pble, g.type="smooth",
                                   fixed.S=NULL, length.S=1000, delta=1*10^-3,
                                   g=function(lambda,bg,E,L,g.type){
                                     g.compute(lambda=lambda,bg=bg,E=E,L=L,g.type=g.type)
                                   }, use.bp=FALSE, use.mix=FALSE, debug=FALSE, verbose=FALSE ){
  
  ####################################################################################################
  # pi.theta.eval.numint   Function evaluates pi(theta,Smin) integral via numerical integration
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
  
  if (verbose) {
    cat("fixed.S:"); print(fixed.S)
  }
  
  # Get break-points
  if (!use.bp){
    bp <- NULL
  } else {
    bp <- bp.grid
  }
  
  # Form S grid based on user-specified min, max and number of points:
  if (is.null(fixed.S)){
    
    if (length.S < 2)
      stop("'length.S' must be greater than or equal to 2")
    
  } else {
    # Override with 1-point grid:
    S.fixed  <- TRUE
    S.grid   <- fixed.S
    length.S <- 1
    dS <- 1.0
  }
  length.S.old <- length.S
  
  ############### version PBLE with [ -accessor, without @ -accessor ###########
  pars <- pble["pars"]
  grid <- pble["grid"]
  # Form E grid based on user-specified min, max and number of points:
  fixed.E  <- pars$fixed.E
  length.E <- pars$length.E
  E.min    <- pars$Emin
  E.max    <- pars$Emax
  dE       <- pars$dE
  E.grid   <- grid$E.grid
  
  # Form L grid based on user-specified min, max and number of points:
  fixed.L  <- pars$fixed.L
  length.L <- pars$length.L
  L.min    <- pars$Lmin
  L.max    <- pars$Lmax
  dL       <- pars$dL
  L.grid   <- grid$L.grid
  
  # Form B grid based on user-specified min, max and number of points:
  fixed.B  <- pars$fixed.B
  length.B <- pars$length.B
  B.min    <- pars$Bmin
  B.max    <- pars$Bmax
  dB       <- pars$dB
  B.grid   <- grid$B.grid
  
  if (verbose){
    cat(paste("Creating (",length.E," x ",length.L," x ",length.B,") grid for numerical integration...\n"),sep="")
    cat(paste("E varies from ",E.min," to ",E.max," at increments of ",dE,"\n",sep=""))
    cat(paste("L varies from ",L.min," to ",L.max," at increments of ",dL,"\n",sep=""))
    cat(paste("B varies from ",B.min," to ",B.max," at increments of ",dB,"\n",sep=""))
  }
  
  BLE.grid <- pble["BLE.grid"]
  bg <- BLE.grid$B
  L  <- BLE.grid$L
  E  <- BLE.grid$E
  if (verbose>1){
    cat("Grid for BLE:\n")
    print(BLE.grid)
  }
  
  if (debug){
    if (verbose) {
      cat("Making large data.frame to store debugging results...\n")
    }
    db <- matrix(NA,nrow=length.S*length.theta,ncol=3+length(BLE.grid))
    colnames(db) <- c("theta","S","fq",ppaste("BLE.",1:length(BLE.grid)))
    db[,"theta"] <- rep(theta,each=length.S)
    if (verbose) {
      cat("dim(db):\n")
      print(dim(db))
      cat("partial colnames(db):\n")
      print(head(colnames(db))) 
    }
  }
  
  if (verbose){
    cat("done. Now looping over Theta...\n")
  }
  
  
  length.theta <- nrow(theta.grid)
  
  ### TODO: apply this case to variable Smin ? May be too much computationally infeasible.
  Smin <- Smin.grid[1]
  p.t <- p.t.grid   # same length as nrow(theta)
  
  for (j in 1:length.theta) {
    
    theta <- as.numeric(theta.grid[j,])
    length.S <- length.S.old
    
    # Irregular grid:
    q.seq <- seq(0,1.0-delta,length.out=(length.S))
    if (use.bp) {
      if(verbose>1){
        #NOTE: Bug found in here: dbrokenpareto at q==0, returns a value slightly <Smin - an impossible value, but is a numeric error: Real error is in qtruncpareto()
        S.grid.q <- qbrokenpareto(q.bp=q.seq,x_min=Smin,k=theta,bp=bp,verbose=verbose)
      } else {
        S.grid.q <- qbrokenpareto(q.bp=q.seq,x_min=Smin,k=theta,bp=bp,verbose=FALSE)
      }
      
    } else if (use.mix) { 
      S.grid.q <- matrix(numeric(1),ncol=m,nrow=length.S)
      for (i in 1:m) {
        S.grid.q[,i] <- qpareto(q=q.seq,x_min=Smin[i],k=theta[i])
      }
      S.grid.q <- sort(S.grid.q)
      length.S <- length(S.grid.q)
      
    } else {
      S.grid.q <- qpareto(q=q.seq,x_min=Smin,k=theta)
    }
    dS.q <- diff(S.grid.q)
    
    if (verbose) {
      cat("k (i.e. theta):\n") ; print(theta)
      if (use.bp && verbose>1) {
        check.sum <- sum(dbrokenpareto(x=S.grid.q[2:length(S.grid.q)],x_min=Smin,k=theta,bp=bp,verbose=verbose)*dS.q)
        cat(paste("sum(qbrokenpareto(S.grid.q))*dS.q = ",check.sum,"\n",sep=""))
      }    
      cat(paste("Irregular S grid:\n")) ; print(S.grid.q)
      cat(paste("dS.q:\n")) ; print(dS.q)
      cat(paste("range(dS.q):\n")) ; print(range(dS.q))
      cat("Computing g at each S over BLE grid...\n")
    }
    
    ### Computing Inner Integral 
    # Irregular grid: 
    g.tilde.q <- rep(NA,times=length.S)
    for (i in 1:length.S) {
      
      if (debug){
        cat(ppaste("Computing g-tilde for S[",i,"] = ",S.grid.q[i],":\n"))
        cat("lambda:\n") ; print(S.grid.q[i]*E/gamma)
        cat("bg:\n") ; print(bg)  
        cat("g.type:\n") ; print(g.type)
      }
      
      #Note: p.BLE is evaluated at exactly N=(length.E-1) points:  E_1,...,E_N.
      g.star <- g(lambda=S.grid.q[i]*E/gamma, bg=bg,L=L,E=E, g.type=g.type) 
      if (class(pble)=="pble.basic") {
        g.tilde.q[i] <- sum(g.star)*pble["prob"]*dB*dL*dE
      } else if (class(pble)=="pble.table") {
        # old code:  
        tmp1 <- sum(g.star*pble["table"]$prob) #*dB*dL*dE
        
        # Trapezoid Rule:
        #pi[j,2] <- (sum(f)-0.5*(f[1]-f[length.S])) *dS*dL*dE / (E.range*L.range) 
        
        # Trapezoid Rule in 3D:
        ix.B.start <- bg == B.min
        ix.B.end   <- bg == B.max
        ix.L.start <- L == L.min
        ix.L.end   <- L == L.max
        ix.E.start <- E == E.min
        ix.E.end   <- E == E.max
        
        f <- g.star*pble["table"]$prob    # inner product
        tmp2 <- sum(f) - 0.5*(sum(f[ix.B.start]) + sum(f[ix.B.end]) + 
                                sum(f[ix.L.start]) + sum(f[ix.L.end]) +
                                sum(f[ix.E.start]) + sum(f[ix.E.end])) +
          0.25*(sum(f[ix.B.start & ix.L.start]) + sum(f[ix.B.start & ix.E.start]) + sum(f[ix.L.start & ix.E.start]) +
                  sum(f[ix.B.start & ix.L.end])   + sum(f[ix.B.start & ix.E.end])   + sum(f[ix.L.start & ix.E.end])   +
                  sum(f[ix.B.end   & ix.L.start]) + sum(f[ix.B.end   & ix.E.start]) + sum(f[ix.L.end   & ix.E.start]) +
                  sum(f[ix.B.end   & ix.L.end])   + sum(f[ix.B.end   & ix.E.end])   + sum(f[ix.L.end   & ix.E.end]))   +
          -0.125*(f[ix.B.start & ix.L.start & ix.E.start] + f[ix.B.end   & ix.L.end   & ix.E.end] +
                    f[ix.B.start & ix.L.start & ix.E.end]   + f[ix.B.start & ix.L.end   & ix.E.end] +
                    f[ix.B.start & ix.L.end   & ix.E.start] + f[ix.B.end   & ix.L.start & ix.E.end] +
                    f[ix.B.end   & ix.L.start & ix.E.start] + f[ix.B.end   & ix.L.end   & ix.E.start])
        g.tilde.q[i] <- tmp1
        
        if(verbose==3){
          if(i==1){
            cat("[Riemann Sum,  3D Trapezoid Rule] of: g.tilde = integral{ g*p(ble) db dl el }:\n")
          }
          print(c(tmp1,tmp2))
        }
      }
      
      if (debug) {
        cat(ppaste("S.grid[i]:",i,":\n")); print(S.grid.q[i])
        cat(ppaste("g.star:",i,":\n")); print(g.star)
        cat("g.tilde.q[i]:\n") ; print(g.tilde.q[i])  
        cat("pble['table']$prob:\n") ; print(pble["table"]$prob)
      }
    } # END for-loop over S.grid
    if (verbose)
      cat("Finished computing g* and g-tilde...\n")
    
    # Sanity Check: if all g.tilde=1, this means all g=1   =>  pi=1
    if (all(g.tilde.q==1)){
      if(verbose){
        cat('All g functions are 1 => all sources are observed here... so pi=1. \n')
      }
      pi <- matrix(1, nrow=length.theta, ncol=3)
      ret <- list("pi"=pi,"theta"=theta)  
      return(ret)
    }
    
    # Evaluate functions on irregular grid: (varying theta part)
    if (use.bp) {
      if(verbose>1){        
        pdf.S.grid.q <- dbrokenpareto(x=S.grid.q, x_min=Smin, k=theta, bp=bp, verbose=verbose)
      } else {
        pdf.S.grid.q <- dbrokenpareto(x=S.grid.q, x_min=Smin, k=theta, bp=bp, verbose=FALSE)
      }
      
    } else if (use.mix) {  
      if(verbose>1){        
        pdf.S.grid.q <- dmixpareto(x=S.grid.q,p=as.numeric(p.t[j,]),x_min=Smin,k=theta,verbose=verbose)
      } else {
        pdf.S.grid.q <- dmixpareto(x=S.grid.q,p=as.numeric(p.t[j,]),x_min=Smin,k=theta,verbose=FALSE)
      }
      #TODO: #Question: should we try to use log-scale to improve numerical error?
      
    } else {
      pdf.S.grid.q <- dpareto(x=S.grid.q, x_min=Smin, k=theta)
    }
    fq <- g.tilde.q*pdf.S.grid.q    
    
    if (debug){
      # If debugging, store the integral over B, L, and E, but not S:
      tmp <- (db[,"theta"]==theta)
      db[tmp,"S"] <- S.grid.q
      db[tmp,"fq"] <- fq 
    } 
    
    if (verbose){
      cat("Finished computing dbrokenpareto/dmixpareto/dpareto...\n")
      cat("g.tilde.q:\n") ; print(g.tilde.q)
      if(verbose>1){
        cat("pdf.S.grid.q:\n") ; print(pdf.S.grid.q)
        cat("fq:\n") ; print(fq)
        cat("dS.q:\n") ; print(dS.q)
        cat("Riemann rectangles: [Low Height, Upper Height, dS]:\n") ; print(cbind(fq[1:(length.S-1)],fq[2:length.S],dS.q))
        cat("Areas of rectangles:\n") ; print(apply(cbind(fq[1:(length.S-1)],fq[2:length.S],dS.q),1,function(x){x[3]*0.5*(x[1]+x[2])}))
      }
      cat("Sum of areas of rectangles (Trapezoid):\n") ; print(sum(apply(cbind(fq[1:(length.S-1)],fq[2:length.S],dS.q),1,function(x){x[3]*0.5*(x[1]+x[2])})))
      if(verbose & class(pble)=="pble.basic"){
        cat("dB * dL * dE * p.BLE:\n") ; print(dB*dL*dE * pble["prob"])
      }
      if(verbose>1){
        cat("Correction term :\n") ; print(sum(g.tilde.q[c(1,length(g.tilde.q))]*delta/nrow(BLE.grid)))
      }
      cat("Length.S:"); print(length.S)
      if(verbose>1){
        cat("Length.E:"); print(length.E)
        cat("Length.B:"); print(length.B)
        cat("Length.L:"); print(length.L)
        cat("d.E:"); print(dE)
        cat("d.L:"); print(dL)
        cat("d.B:"); print(dB)
      }
    }
    
    ### Integrate 
    # Irregular grid integral:
    iv <- sum(apply(cbind(fq[1:(length.S-1)],fq[2:length.S],dS.q),1,function(x){x[3]*0.5*(x[1]+x[2])})) 
    # Correction for the bits chopped off at upper limit, still accounting for g
    pi[j] <- iv + sum(g.tilde.q[c(1,length(g.tilde.q))]*delta) #/nrow(BLE.grid)) 
    
    if (verbose){
      cat("Finished trapezoidal integration over irregular grid...\n")
      cat('pi:\n'); print(pi[j])
    }     
    
    if (use.bp) {
      cat(c("Computed pi(theta,Smin) for Smin = ",ppaste(Smin),"; theta = ",ppaste(round(theta,5))," (",j,"/",length.theta,") ; Pr(Obs) = ",round(pi[j],5),"\n"))            
    } else if (use.mix) {
      cat(c("Computed pi(theta,Smin) for Smin = ",ppaste(Smin),"; theta = ",ppaste(round(theta,5))," (",j,"/",length.theta,") ; p.t = ",ppaste(p.t[j,]), " ; Pr(Obs) = ",round(pi[j],5),"\n"))            
    } else {   
      cat(paste("Computed pi(theta,Smin) for Smin = ",ppaste(Smin),"; theta = ",theta," (",j,"/",length.theta,") ; Pr(Obs) = ",round(pi[j],5),"\n",sep=""))      
    }
    
  } # END loop over Theta
  
  if(verbose) {
    cat('pi:\n'); print(pi)
  }
  
  ret <- list("pi"=pi)
  
  # If debugging, append the component grid:
  if (debug){
    ret$db.grid.frame <- db
    ret$db.grid.theta <- theta
  }
  
  ################## ################## ###############
  return(ret)
  
}
