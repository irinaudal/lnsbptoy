"plot.lns" <- function(resmc,logS.truth=NULL,main=NULL,max.lines=1000,cex=0.6,
                       col.obs="darkgrey",col.mis="red",col.truth=rgb(.13,.29,.99),transp.col=NULL)
{
  ####
  # plot.lns    Plot LogN-LogS curve of the posterior samples of fluxes.
  #
  # Input: resmc      = large object of (1) "draws": matrix of posterior draws (rows: MCMC iterations, cols: variables including S.obs)
  #                                     (2) "draws.S.mis": list of posterior draws of missing sources (optional)
  #        logS.truth = true vector of log(flux)
  #        main       = title for the plot
  #        max.lines  = maximum number of posterior iterations to plot
  #        cex        = size of plot points
  #        col.obs    = color of posterior samples of observed fluxes
  #        col.mis    = color of posterior samples of missing fluxes
  #        col.truth  = color of true flux.  rgb(.13,.29,.99) = lighter blue color
  #        transp.col = NULL/"obs"/"mis" make a subset of fluxes transparent (i.e. do not plot)
  ####
  
  
  # Set color 'col' to hex notation
  col.obs.blank <- paste(rgb( t( col2rgb(col.obs)/255 ) ),"00",sep='')
  col.mis.blank <- paste(rgb( t( col2rgb(col.mis)/255 ) ),"00",sep='')
  
  logSobs <- mcmc(log(resmc$draws[,grep("S.",varnames(resmc$draws)),drop=FALSE],base=10))
  logSmis <- lapply(resmc$draws.S.mis,log,base=10)
  
  # Plot loads of faint lines for each MCMC iteration...
  
  # First plot the graph...    
  if (is.null(main)){
    main <- "log(N>s) vs. log(s): Posterior Draws"
  }
  n <- ncol(logSobs)
  N.mis <- sapply(logSmis,length)
  N.mis.max <- max(N.mis)
  N.com <- n+N.mis 
  N.com.max <- max(N.com)
  
  logNgS.max <- log(N.com.max*(1-seq(0,(N.com.max-1)/N.com.max,length.out=N.com.max+1)),base=10)
  
  S.obs.range <- range(logSobs)
  S.mis.range <- range(unlist(logSmis))
  
  # Don't add too many lines...
  niter <- nrow(logSobs)
  if (niter > max.lines){
    # Only plot max.lines iterations:
    iters.todo <- floor(seq(from=1,to=niter,length.out=max.lines))
  } else { 
    # Can plot all of them:
    iters.todo <- c(1:niter)
  }
  
  # First plot the graph...
  if(!all(N.mis==0)){
    plot(x=NULL,type="n",xlim=range(c(S.obs.range,S.mis.range)),ylim=range(logNgS.max),cex=cex,
         main=main,xlab="log(s)",ylab="log(N>s)")
  } else {
    plot(x=NULL,type="n",xlim=range(c(S.obs.range)),ylim=range(logNgS.max),cex=cex,
         main=main,xlab="log(s)",ylab="log(N>s)")
  }  
  
  # Now add a line for each iteration...
  for (i in iters.todo){
    logNgS <- log(N.com[i]*(1-seq(0,(N.com[i]-1)/N.com[i],length.out=N.com[i]+1)),base=10)
    
    logScom <- c(logSobs[i,],logSmis[[i]])
    I <- rep(c(1,0),times=c(n,N.mis[i]))
    
    slog <- sort(logScom,index.return=TRUE)
    logS.tmp <- slog$x
    s.lns <- stepfun(x=logS.tmp,y=logNgS)
    
    if (!is.null(I)){
      # Color the missing data points in red...
      if (length(I)!=N.com[i]){
        stop("'I' must be same length as 'logS'")
      }
      I.star <- I[slog$ix] 
      if (is.null(transp.col)){
        colvec <- ifelse(I.star,col.obs,col.mis)
      }else if(transp.col=="obs"){ 
        colvec <- ifelse(I.star,col.obs.blank,col.mis)
      }else if(transp.col=="mis"){
        colvec <- ifelse(I.star,col.obs,col.mis.blank)
      }
      # Make the first node the transparent color
      colvec[1] <- col.mis.blank  #set the start and end horizontal line to transparent
    } else {
      colvec <- rep(col.obs,n)
    }
    lines(s.lns,col=colvec,lwd=0.1,cex=0.1)
  } # End iteration loop	
  
  if (!is.null(logS.truth)){
    #ADD BLUE LINE of truth
    # Just a vector of S's, plot a single logN-logS
    slog <- sort(logS.truth,index.return=TRUE)
    logS <- slog$x
    n <- length(logS)
    logNgS <- log(n*(1-seq(0,(n-1)/n,length.out=n+1)),base=10)
    s.lns <- stepfun(x=logS,y=logNgS)
    
    lines(s.lns,yaxt="n",col=col.truth,xlab="log(s)",ylab="log(N>s)",cex=0.75,lwd=1.8)
  }
  
  return(NULL)
}
