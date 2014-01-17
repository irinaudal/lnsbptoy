"pareto_lns" <- function(theta, x_min, p=NULL, length.x=1000,delta=0.01, verbose=FALSE, 
                         outdir=FALSE, smin.color="red", smin.lwd=1.8, smin.cex=1)  
{
  ##############
  # pareto_lns     Make logN-logS plot for any pareto configuration
  #
  # Input: theta    = vector of power-law slopes
  #        x_min    = vector of ALL flux boundary conditions: S_min, break-points, or mixture S_min's
  #        p        = optional vector of mixture proportions
  #        length.x = size of grid for Pareto samples of fluxes
  #        delta    = upper bundary Pareto quantile for grid of fluxes
  #        verbose  = (T/F) display progress of program
  #        outdir      = FALSE/filename : plot to screen or to file
  #        smin.color  = plotting color for missing fluxes
  #        smin.lwd    = plotting line width for missing sources
  #        smin.cex    = plotting point size for missing sources
  ##############
  
  use.mix=FALSE; use.bp=FALSE
  x_min <- sort(x_min)
  len <- length(x_min)
  m <- length(theta)
  q.seq <- seq(0,1.0-delta,length.out=length.x)
  
  if (len!=m){
    stop("Length of x_min vector and theta vector must be of the same size.\n")
  }
  
  if (len==1){
    # Regular Pareto:
    x <- qpareto(q=q.seq, x_min=x_min, k=theta)
    y <- dpareto(x=x, x_min=x_min, k=theta)
    w <- log(1-ppareto(x=x, x_min=x_min, k=theta), base=10)
    model <- "Pareto"
    model.short <- "pareto"
    
  } else if (len>1){
    # Broken- or Mixture-Pareto:
    if (!is.null(p)){
      # Mixture-Pareto:
      use.mix <- TRUE
      
      x <- matrix(numeric(1),ncol=m,nrow=length.x)
      y <- x
      for (i in 1:(m-1)) {
        x[,i] <- qpareto(q=q.seq, x_min=x_min[i], k=theta[i])
      }
      i <- m
      x[,i] <- qpareto(q=q.seq[1:floor(length.x/2)], x_min=x_min[i], k=theta[i])
      
      x <- sort(x)
      y <- dmixpareto(x=x, p=p, x_min=x_min, k=theta, verbose=verbose)
      w <- log(1-pmixpareto(x=x, p=p, x_min=x_min, k=theta, verbose=verbose), base=10)
      model <- "Mixture-Pareto"
      model.short <- "mix"
      
    } else {
      # Broken-Pareto:
      use.bp <- TRUE
      bp <- x_min[-which.min(x_min)]
      Smin <- min(x_min)
      
      x <- qbrokenpareto(q.bp=q.seq, x_min=Smin, k=theta, bp=bp, verbose=verbose)
      x2 <- qpareto(q=q.seq[1:floor(length.x/2)], x_min=x_min[len], k=theta[len])   # BAD fix
      x <- sort(c(x,x2))
      y <- dbrokenpareto(x=x, x_min=Smin, k=theta, bp=bp, verbose=verbose)
      w <- log(1-pbrokenpareto(x=x, x_min=Smin, k=theta, bp=bp, verbose=verbose), base=10)
      model <- "Broken-Pareto"
      model.short <- "bp"
    }    
  } else {
    stop("Length of x_min must be a positive integer.\n")
  }
  w <- w+abs(w[length(w)]) # adjust all values above 0
  log.x <- log(x, base=10)
  log.x_min <- log(x_min, base=10)
  
  if (verbose){
    cat("pareto_lns() gives:\n")
    cat("x=\n"); print(x)
    cat("y=\n"); print(y)
    cat("w=\n"); print(w)
    cat("model=\n"); print(model)
  }
  
  # Plot density curve
  
  if (is.character(outdir)) {
    pdf(ppaste(outdir,"view_density_",model.short,".pdf"))
  }
  
  plot(x,y,type="l",lwd=smin.lwd, xlab="S", ylab="P(S>Smin)", main=ppaste(model," density"))
  abline(v=x_min, col=smin.color)
  
  if (is.character(outdir)) {
    dev.off()
  }
  
  # Plot LNS curve
  
  if (is.character(outdir)) {
    pdf(ppaste(outdir,"view_lns_",model.short,".pdf"))
  }
  
  colvec <- rep("black",length.x)
  colvec[1] <- paste(rgb( t( col2rgb("black")/255 ) ),"00",sep='') #set the start and end horizontal line to transparent
  
  plot(x=NULL,type="n",xlim=range(range(log.x)),ylim=range(w),cex=0.6,
       main=ppaste(model,": log(N>s) vs. log(s)"),xlab="log(s)",ylab="log(N>s)")
  s.lns <- stepfun(x=log.x, y=c(w[1],w))
  lines(s.lns,col=colvec,lwd=smin.lwd,cex=smin.cex)
  abline(v=log.x_min, col=smin.color)
  
  if (is.character(outdir)) {
    dev.off()
  }
  
  return(NULL)
}