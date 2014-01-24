"trace.mcmc.plot" <- function(resmc, res, true.par, max.mcmc.plots="all", remote.copy=FALSE) {

  if (all(!resmc$inputs$fixed.theta)) {
    # Make sure lengths of theta match
    id.theta <- grep("theta",colnames(resmc$draws))
    nt.mcmc <- length(id.theta)
    nt.true <- length(res$par$theta)
    if (nt.mcmc>nt.true) {
      # Add NA to the true.par vector of theta
      tmp <- res$par
      tmp$theta[(nt.true+1):nt.mcmc] <- NA
      true.par <- unlist(tmp)
    }    
  }
  
  # Grab variable index to be plotted
  id.logp <- grep("log.Poster",colnames(resmc$draws))
  id.S.obs <- grep("S.obs",colnames(resmc$draws))
  if (max.mcmc.plots=="all") {
    if (length(id.logp>0)) {
      index <- 1:ncol(resmc$draws[,-id.logp])
    } else {
      index <- 1:ncol(resmc$draws)
    }
  } else if (is.numeric(max.mcmc.plots)) {
    index <- 1:min(max.mcmc.plots,length(id.S.obs+2))
  } else {
    index <- 1:(min(id.S.obs)-1)
  }
  
  # Make trace plots
  for( i in index){
    tmp.parname <- varnames(resmc$draws)[i]
    tmp.truth <- true.par[tmp.parname]
    traceplot(resmc$draws[,tmp.parname],main=paste("Posterior draws: ",tmp.parname,sep=""),ylab=tmp.parname)
    densplot(resmc$draws[,tmp.parname],main=paste("Posterior density: ",tmp.parname,sep=""),xlim=c(min(c(tmp.truth,resmc$draws[,tmp.parname])),max(c(tmp.truth,resmc$draws[,tmp.parname]))))
    abline(v=tmp.truth,col="blue",lwd=1.8)
  }
  
  # Copy file to another server and delete the local copy:
  if (remote.copy){
    ## TODO: Try to patch this better!!!
    scp.copy(from.filename=output.mcmcplot.hack,to.filename=remote.mcmcplot,from.delete=from.delete,recursive.dir.create=recursive.dir.create,verbose=TRUE)
  }
  
} # END trace.mcmc.plot()


"posterior.theta" <- function(resmc, res, ci.vals=c(0.025,0.975)) {
  
  n.theta <- length(resmc$inputs$a) #number of theta parameters
  model   <- resmc$inputs$model
  
  for(j in 1:n.theta){
    td <- mcmc(resmc$draws[,ppaste("theta.",j)])
    sorted.td <- sort(td)
    td.ci <- sorted.td[length(td)*ci.vals]    
    densplot(td,main=ppaste("Posterior distribution of theta",j))
    abline(v=res$par$theta[j],col="blue",lwd=1.5)
    abline(v=td.ci,col="red",lwd=1.5)
    text(x=res$par$theta[j]+max(res$par$theta[j])*0.02,y=0.25,labels=expression(theta),col="blue")
    text(x=td.ci+max(td.ci)*0.02,y=0.25,labels=paste(as.character(ci.vals),"%",sep=""),col="red",cex=0.75)
    cat("Theta Quantiles :\n"); print(quantile(td,probs=c(.5,.9)))
  }
  # Clean-up:
  if (exists("td")){
    rm(td)
  }
  if (exists("td.ci")){
    rm(td.ci)
  }
  
} # END posterior.theta()


"posterior.scatterplot" <- function(resmc, res) {
  
  n.theta    <- length(resmc$inputs$a) #number of theta parameters
  fixed.Smin <- resmc$inputs$fixed.Smin # indicator if Smin is fixed
  fixed.bp   <- resmc$inputs$fixed.bp # indicator if bp is fixed
  model      <- resmc$inputs$model
  
  # High Density Scatterplots with Color Transparency 
  plot(cbind(resmc$draws[,2],resmc$draws[,"N"]),xlab="Theta",ylab="N", main="Scatterplot of Posterior Draws", col=rgb(0,0,0,50,maxColorValue=255), pch=16)  
  if(!fixed.Smin){
    plot(cbind(resmc$draws[,"tau.1"],resmc$draws[,"N"]),xlab="Tau.1",ylab="N", main="Scatterplot of Posterior Draws", col=rgb(0,0,0,50,maxColorValue=255), pch=16)
    if (model=="bp" && all(!fixed.bp)){
      plot(cbind(resmc$draws[,"tau.1"],resmc$draws[,"tau.2"]),xlab="Tau.1",ylab="Tau.2", main="Scatterplot of Posterior Draws", col=rgb(0,0,0,50,maxColorValue=255), pch=16)
      plot(cbind(resmc$draws[,"tau.2"],resmc$draws[,"theta.1"]),xlab="Tau.2",ylab="Theta.1", main="Scatterplot of Posterior Draws", col=rgb(0,0,0,50,maxColorValue=255), pch=16)
    }
    plot(cbind(resmc$draws[,"tau.1"],resmc$draws[,2]),xlab="tau.1",ylab="Theta.1", main="Scatterplot of Posterior Draws", col=rgb(0,0,0,50,maxColorValue=255), pch=16)
  } else {
    if(model=="bp" && all(!fixed.bp)){
      plot(cbind(resmc$draws[,"tau.2"],resmc$draws[,"theta.1"]),xlab="Tau.2",ylab="Theta.1", main="Scatterplot of Posterior Draws", col=rgb(0,0,0,50,maxColorValue=255), pch=16)
    }
  }
  # plot bivariate plot of theta, if dim(theta)==2
  if (n.theta==2){
    plot(cbind(resmc$draws[,"theta.1"],resmc$draws[,"theta.2"]),xlab="Theta1",ylab="Theta2", main="Scatterplot of Posterior Draws of Theta", col=rgb(0,0,0,50,maxColorValue=255), pch=16)
  }
  
} # END posterior.scatterplot()


"prior.plots" <- function (resmc, true.par) {

  model      <- resmc$inputs$model
  
  #### prior of Smin
  if (model=="bp") {
    mu    <- resmc$inputs$mu
    C     <- resmc$inputs$C
    n.K   <- length(mu)
    
    cat("BP priors plot function is not written yet.\n")
    
  } else {
    am      <- resmc$inputs$am
    bm      <- resmc$inputs$bm
    
    if(!is.null(am) && !is.null(bm)){
      #Smin.prior.mean <- 0.9*10^(-17) ;  Smin.prior.var  <- (8.0*10^(-17))^2; am  <- Smin.prior.mean^2 / Smin.prior.var  ;  bm  <- am / Smin.prior.mean
      #Smin.prior.mean <- 1.5*10^(-13) ;  Smin.prior.var  <- (9.2*10^(-14))^2 ; am  <- Smin.prior.mean^2 / Smin.prior.var  ;  bm  <- am / Smin.prior.mean    
      #Smin.prior.mean <- 5.0*10^(-13) ;  Smin.prior.var  <- (4.2*10^(-13))^2 ; am  <- Smin.prior.mean^2 / Smin.prior.var  ;  bm  <- am / Smin.prior.mean    
      
      Smin.plot.seq <- qgamma(p=seq(0.001,0.99,length.out=10000),shape=am,rate=bm)
      pSmin.plot.seq <- dgamma(x=Smin.plot.seq,shape=am,rate=bm)
      plot(x=Smin.plot.seq,y=pSmin.plot.seq,main="Prior for Smin",xlab="tau.1",ylab="p(tau.1)",type="l")
      abline(v=true.par["tau.1"],col=2)
      legend("topright",legend=paste("shape=",am,", rate=",bm,sep=""))
      vvv <- c( Smin.plot.seq[1], Smin.plot.seq[92], Smin.plot.seq[10000])
      print(vvv)
      abline(v=vvv,col=3)
    }
  }
  
  #### prior of N
  alpha     <- resmc$inputs$alpha
  bbeta     <- resmc$inputs$beta
  
  N.plot.lo <- qnbinom(p=0.001,size=alpha,prob=beta/(1+beta))
  N.plot.hi <- qnbinom(p=0.999,size=alpha,prob=beta/(1+beta))
  N.plot.seq <- seq(N.plot.lo,N.plot.hi,by=1)
  pN.plot.seq <- dnbinom(x=N.plot.seq,size=alpha,prob=beta/(1+beta))
  
  plot(y=pN.plot.seq,x=N.plot.seq,main="Prior for N",xlab="N",ylab="p(N)")
  legend("topright",legend=paste("size=",alpha,", prob=",beta/(1+beta),sep=""))
  abline(v=true.par["N"],col="red")
  
  #### prior of theta
  a       <- resmc$inputs$a
  b       <- resmc$inputs$b
  n.theta <- length(a)
  
  theta.plot.lo <- qgamma(p=0.001,shape=a,rate=b)
  theta.plot.hi <- qgamma(p=0.999,shape=a,rate=b)
  
  for(j in 1:n.theta){
    if(length(a)>1){
      plot(function(x){dgamma(x=x,shape=a[j],rate=b[j])},theta.plot.lo[j],theta.plot.hi[j],main=ppaste("Prior for theta.",j),xlab=ppaste("theta.",j),ylab=ppaste("p(theta.",j,")"))
    } else {
      plot(function(x){dgamma(x=x,shape=a[j],rate=b[j])},theta.plot.lo,theta.plot.hi,main="Prior for theta",xlab="theta",ylab="p(theta)")
    }
    abline(v=true.par[1+j],col="red")
    legend("right",legend=paste("shape=",a[j],", rate=",b[j],sep=""))
  }
} # END prior.plots()


"posterior.prior.plots" <- function (resmc, true.par) {
  
  model      <- resmc$inputs$model
  
  #### prior of Smin
  if (model=="bp") {
    mu    <- resmc$inputs$mu
    C     <- resmc$inputs$C
    n.K   <- length(mu)
    
    cat("BP posterior+prior plot function is not written yet.\n")
    
  } else {
    am      <- resmc$inputs$am
    bm      <- resmc$inputs$bm
    
    if(!is.null(am) && !is.null(bm)){
      #### plot prior+posterior for Smin:
      Smin.plot.lo <- qgamma(p=0.001,shape=am,rate=bm)
      Smin.plot.hi <- qgamma(p=0.999,shape=am,rate=bm)
      Smin.plot.seq <- seq(Smin.plot.lo,Smin.plot.hi,length.out=1000)
      pSmin.plot.seq <- dgamma(x=Smin.plot.seq,shape=am,rate=bm)
      Smin <- resmc$draws[,"tau.1"]
      Smin_density <- density(Smin)
      xlims <- c(min(Smin.plot.lo, Smin, true.par["tau.1"]),max(Smin.plot.hi, Smin, true.par["tau.1"]))
      ylims <- c(0,max(Smin_density$y))
      plot(y=pSmin.plot.seq,x=Smin.plot.seq,main="Prior for tau.1",xlab="tau.1",ylab="p(tau.1)",type="n",xlim=xlims,ylim=ylims)
      ### ADD shading of 95% Credible Interval of posterior of Smin
      area <- shade.area.under.curve(Smin,q_lo=0.025,q_hi=0.975,nbins=4000,col="lightgrey")
      ### ADD the posterior of Smin as a kernel density estimate
      lines(Smin_density)
      ### ADD the prior of Smin
      lines(y=pSmin.plot.seq,x=Smin.plot.seq,lty="dashed",col="blue")
      ### ADD true N parameter
      abline(v=true.par["tau.1"],col="red")
      legend("topright",legend=c("posterior","prior"),lty=c(1,2),col=c(1,"blue"))
    }
  }
  
  ##### plot prior+posterior for N:
  alpha     <- resmc$inputs$alpha
  beta     <- resmc$inputs$beta
  
  N.plot.lo <- qnbinom(p=0.001,size=alpha,prob=beta/(1+beta))
  N.plot.hi <- qnbinom(p=0.999,size=alpha,prob=beta/(1+beta))
  N.plot.seq <- seq(N.plot.lo,N.plot.hi,by=1)
  pN.plot.seq <- dnbinom(x=N.plot.seq,size=alpha,prob=beta/(1+beta))
  
  N <- resmc$draws[,"N"]
  N_density <- density(N)
  xlims <- c(min(N.plot.lo, N, true.par["N"]),max(N.plot.hi, N, true.par["N"]))
  ylims <- c(0,max(N_density$y))
  plot(y=pN.plot.seq,x=N.plot.seq,main="Prior for N",xlab="N",ylab="p(N)",type="n",xlim=xlims,ylim=ylims)
  ### ADD shading of 95% Credible Interval of posterior of N
  area <- shade.area.under.curve(N,q_lo=0.025,q_hi=0.975,nbins=4000,col="lightgrey")
  ### ADD the posterior of N as a kernel density estimate
  lines(N_density)
  ### ADD the prior of N
  lines(y=pN.plot.seq,x=N.plot.seq,lty="dashed",col="blue")
  ### ADD true N parameter
  abline(v=true.par["N"],col="red")    
  legend("topleft",legend=c("posterior","prior"),lty=c(1,2),col=c(1,"blue"))
  
  ### Prior+posterior plot for theta
  a       <- resmc$inputs$a
  b       <- resmc$inputs$b
  n.theta <- length(a)
  
  theta.plot.lo <- qgamma(p=0.001,shape=a,rate=b)
  theta.plot.hi <- qgamma(p=0.999,shape=a,rate=b)
  
  n.theta <- length(a)  
  for(j in 1:n.theta){     
    theta <- resmc$draws[,1+j]
    th_density <- density(theta)
    xlims <- c(min(theta.plot.lo[j], theta, true.par[1+j]),max(theta.plot.hi[j], theta, true.par[1+j]))
    ylims <- c(0,max(th_density$y))
    if(length(a)>1){
      plot(function(x){dgamma(x=x,shape=a[j],rate=b[j])},theta.plot.lo[j],theta.plot.hi[j],main=ppaste("Prior for theta.",j),xlab=ppaste("theta.",j),ylab=ppaste("p(theta.",j,")"),type="n",xlim=xlims,ylim=ylims)
    } else {
      plot(function(x){dgamma(x=x,shape=a[j],rate=b[j])},theta.plot.lo,theta.plot.hi,main="Density for theta",xlab="theta",ylab="p(theta)",type="n",xlim=xlims,ylim=ylims)
    }
    ### ADD shading of 95% Credible Interval of posterior of theta
    area <- shade.area.under.curve(theta,q_lo=0.025,q_hi=0.975,nbins=4000,col="lightgrey")
    ### ADD the posterior of theta as a kernel density estimate
    lines(th_density)
    ### ADD the prior of theta
    curve <- function(x){dgamma(x=x,shape=a[j],rate=b[j])}
    bc <- seq(theta.plot.lo[j],theta.plot.hi[j],length.out=500)
    yv <- curve(bc)
    lines(bc,yv,lty="dashed",col="blue")
    ### ADD true theta parameter
    abline(v=true.par[1+j],col="red")
    legend("topleft",legend=c("posterior","prior"),lty=c(1,2),col=c(1,"blue"))
  }
  

} # END posterior.prior.plots()
