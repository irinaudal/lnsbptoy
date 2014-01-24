# 
# "rpareto" <- function(n, x_min, k)
# {
#    .Call("rpareto",
#    as.integer(n),as.numeric(x_min),as.numeric(k),
#    package="logNlogS") 
# }
# 
# "dpareto" <- function(x, x_min, k, log=FALSE)
# {
#    .Call("dpareto",
# 	 as.numeric(x),as.numeric(x_min),as.numeric(k),as.logical(log),
# 	 package="logNlogS") 
# }
# 
# "ppareto" <- function(x, x_min, k, check.args=TRUE, log=FALSE)
# {
#    .Call("ppareto",
# 	 as.numeric(x),as.numeric(x_min),as.numeric(k),as.logical(check.args),as.logical(log),
# 	 package="logNlogS") 
# }
# 
# "qpareto" <- function(q, x_min, k)
# {
#    .Call("qpareto",
# 	 as.numeric(q),as.numeric(x_min),as.numeric(k),
# 	 package="logNlogS") 
# }
# 
# # Not in NAMESPACE -- not to be called separately
# "truncpareto.check" <- function(x_min, x_max, k, fn, verbose=FALSE)
# {
#    .Call("truncpareto_check",
# 	 as.numeric(x_min),as.numeric(x_max),as.numeric(k),as.character(fn),as.logical(verbose),
# 	 package="logNlogS")
# }
# 
# "rtruncpareto" <- function(n, x_min, x_max, k)
# {
#    .Call("rtruncpareto",
# 	 as.integer(n),as.numeric(x_min),as.numeric(x_max),as.numeric(k),
# 	 package="logNlogS")
# }
# 
# "dtruncpareto" <- function(x, x_min, x_max, k, log=FALSE, verbose=FALSE)
# {
#    .Call("dtruncpareto",
# 	 as.numeric(x),as.numeric(x_min),as.numeric(x_max),as.numeric(k),as.integer(log),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# "ptruncpareto" <- function(x, x_min, x_max, k, log=FALSE)
# {
#    .Call("ptruncpareto",
# 	 as.numeric(x),as.numeric(x_min),as.numeric(x_max),as.numeric(k),as.integer(log),
# 	 package="logNlogS")
# }
# 
# "qtruncpareto" <- function(q, x_min, x_max, k)
# {
#    .Call("qtruncpareto",
# 	 as.numeric(q),as.numeric(x_min),as.numeric(x_max),as.numeric(k),
# 	 package="logNlogS")
# }
# 
# 
# ############################################################
# ##
# ## *brokenpareto functions:
# ##
# ## Functions for broken power laws with M break-points.
# ##
# ## Density, RNG, quantiles and CDF for a
# ## mixture of M truncated-Pareto's and one
# ## untruncated Pareto. Note that there are M+1 mixture
# ## components when there are M break-points:
# ##
# ## Y = p_0 X_0 + ... + p_{M-1} X_{M-1} + p_{M} X_{M}
# ##
# ## where:
# ##
# ## \sum_{j=0}^{M} p_{j} = 1.0
# ##
# ## p_{0} = 1-(bp_{1}/S_{min})^{-\theta_{0}}
# ## p_{j} = (1-(bp_{j+1}/bp_{j})^{-\theta_{j}})*(1-\sum_{i=0}^{j-1}p_{i})
# ## 
# ## and:
# ##
# ## X_{j} ~ Trunc-Pareto( \theta_{j} , [bp_{j} , bp_{j+1}] )
# ##
# ## where:
# ##
# ## for j=1,...,M-1, and:
# ##
# ## X_{M} ~ Pareto( \theta_{M} , bp_{M} ) .
# ##
# ##
# ############################################################
# 
# "rbrokenpareto" <- function(n, x_min, k, bp, verbose=FALSE)
# {
#    .Call("rbrokenpareto",
# 	 as.integer(n),as.numeric(x_min),as.numeric(k),as.numeric(bp),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# "dbrokenpareto" <- function(x, x_min, k, bp, log=FALSE, verbose=FALSE)
# {
#    .Call("dbrokenpareto",
# 	 as.numeric(x),as.numeric(x_min),as.numeric(k),as.numeric(bp),as.integer(log),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# 
# "pbrokenpareto" <- function(x, x_min, k, bp, log=FALSE, verbose=FALSE)
# {
#    .Call("pbrokenpareto",
# 	 as.numeric(x),as.numeric(x_min),as.numeric(k),as.numeric(bp),as.integer(log),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# "qbrokenpareto" <- function(q, x_min, k, bp, verbose=FALSE)
# {
#    .Call("qbrokenpareto",
# 	 as.numeric(q),as.numeric(x_min),as.numeric(k),as.numeric(bp),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# # Not in NAMESPACE -- C-version does not need to interface with R
# 
# "brokenpareto.posterior" <- function(a,b,S.obs,S.mis,Smin,bp,verbose=FALSE)
# {
#    .Call("brokenpareto_posterior",
# 	 as.numeric(a),as.numeric(b),as.numeric(S.obs),as.numeric(S.mis),as.numeric(Smin),as.numeric(bp),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# ############################################################
# ##
# ## *mixpareto functions:
# ##
# ## Functions for mixture of power laws with M different Smin-points.
# ##
# ## Density, RNG, quantiles and CDF for a
# ## mixture of M Pareto's. 
# ##
# ## Y = p_1 X_1 + ... + p_{M-1} X_{M-1} + p_{M} X_{M}
# ##
# ## where:
# ##
# ## \sum_{j=1}^{M} p_{j} = 1.0
# ##
# ## X_{j} ~ Pareto( \theta_{j} , Smin_{j} ) for j=1,...,M.
# ##
# ############################################################
# 
# 
# "rmixpareto" <- function(n, p, k, x_min, fixed.I.idx=NULL)
# {
#   # Note: fixed.I.idx becomes integer(0) when passed to C-code, but this is fine,
#   # as the C-code checks length(fixed.I.idx) to determine NULL-ness
#   .Call("rmixpareto",
# 	 as.integer(n),as.numeric(p),as.numeric(k),as.numeric(x_min),as.integer(fixed.I.idx),
# 	 package="logNlogS")
# }
# 
# "dmixpareto" <- function(x, p, x_min, k, log=FALSE, verbose=FALSE)
# {
#    .Call("dmixpareto",
# 	 as.numeric(x),as.numeric(p),as.numeric(x_min),as.numeric(k),as.integer(log),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# "pmixpareto" <- function(x, p, x_min, k, log=FALSE, verbose=FALSE)
# {
#    .Call("pmixpareto",
# 	 as.numeric(x),as.numeric(p),as.numeric(x_min),as.numeric(k),as.integer(log),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# "qmixpareto" <- function(q, p, x_min, k, tol=1e-5, n.tries.max=2000, use.C=TRUE, verbose=FALSE)
# {
#     # TODO: Remove use_C argument when we fully switch to C code:
#    .Call("qmixpareto",
# 	 as.numeric(q),as.numeric(p),as.numeric(x_min),as.numeric(k),as.numeric(tol),as.integer(n.tries.max),as.integer(verbose),
# 	 package="logNlogS")
# }
# 
# # Not in NAMESPACE -- C-version does not need to interface with R
# 
# "mixpareto.posterior" <- function(a,b,S.obs,S.mis,Smin,I.idx,verbose=FALSE)
# {
#    .Call("mixpareto_posterior",
# 	 as.numeric(a),as.numeric(b),as.numeric(S.obs),as.numeric(S.mis),as.numeric(Smin),as.integer(I.idx),as.integer(verbose),
# 	 package="logNlogS")
# }

######################## OLD R CODE ########################


"rpareto" <- function(n, x_min, k, verbose=FALSE)
{
  if(verbose){
    cat("rpareto()...\n n:\n"); print(n)
    cat("s_min:\n"); print(x_min)
    cat("slope:\n"); print(k)
  }
  
  if (any(k<0)){
    k[k<0] <- NA
    warning("'k' must be positive.\n") 
  }
  if (any(x_min<0)){
    x_min[x_min<0] <- NA
    warning("'x_min' must be positive.\n")
  }
  
  return(x_min/(runif(n)^(1/k)))
}
rpareto <- cmpfun(rpareto)

"dpareto" <- function(x, x_min, k, log=FALSE)
{
  if (length(log)>1)
    stop("Log must be TRUE/FALSE")
  if (any(k<0)){
    k[k<0] <- NA
    warning("'k' must be positive.\n") 
  }
  if (any(x_min<0)){
    x_min[x_min<0] <- NA
    warning("'x_min' must be positive.\n") 
  }
  
  # Allow result to be of length x, of k, or of x_min. Vectors recycled.  
  if(log){
    if (any(x<0)){ #Supress warning message of evaluating log(neg.value)
      x[x<0] <- NA 
    }
    ret <- log(k/x_min)-(k+1)*log(x/x_min)
    
    ret[x<x_min] <- -Inf
    
  } else {
    ret <- (k/x_min)*((x/x_min)^(-(k+1)))  
    
    ret[x<x_min] <- 0
  }
  
  return(ret)
}
dpareto <- cmpfun(dpareto)


"ppareto" <- function(x, x_min, k, check.args=TRUE, log=FALSE)
{
  if (check.args){
    if (length(log)>1)
      stop("Log must be TRUE/FALSE")
    if (any(k<0)){
      k[k<0] <- NA
      warning("'k' must be positive.\n") 
    }
    if (any(x_min<0)){
      x_min[x_min<0] <- NA
      warning("'x_min' must be positive.\n") 
    }
  }
  
  # Allow result to be of length x, of k, or of x_min. Vectors recycled.  
  if (log){
    ret <- log(1-(x_min/x)^k)
    ret[x<x_min | x<0] <- -Inf
  } else {
    ret <- 1-(x_min/x)^k
    ret[x<x_min | x<0] <- 0
  }
  
  return(ret)
}
ppareto <- cmpfun(ppareto)

"qpareto" <- function(q, x_min, k)
{
  ret <- x_min / ((1-q)^(1/k))
  ret[(q<0) | (q>1) | x_min<0 | k<0] <- NA
  if (any(x_min<0)) { warning("'x_min' must be positive.\n") }
  if (any(k<0))     { warning("'k' must be positive.\n") }
  return(ret)  
}
qpareto <- cmpfun(qpareto)

####

"truncpareto.check" <- function(x_min, x_max, k, fn, verbose=FALSE)
{
  if (length(x_min) != length(x_max)){
    stop(ppaste("'x_min' and 'x_max' must be the same length in '",fn,"'"))
  }
  
  if (verbose){
    cat(ppaste("Debugging inside truncpareto.check called by '",fn,"'\n"))
    cat("x_min:\n") ; print(x_min)
    cat("x_max:\n") ; print(x_max)
    cat("k:\n")     ; print(k)
  }
  
  if (any(is.na(k)) | any(k<0)){
    k.NA.indices <- is.na(k)
    k.0.indices <- (k<0)
    if (any(k.NA.indices)){
      stop(ppaste("'k' cannot be NA in '",fn,"'")) 
    }
    if (!any(k.NA.indices) & any(k.0.indices)){
      warning(ppaste("'k' must be positive in '",fn,"'")) 
    }
    k[k.NA.indices] <- NA
    k[k.0.indices] <- NA
  }
  
  if (any(x_min<0) || any(x_min>x_max)){ # PDB 4/10/12: Changed to any(x_min>x_max)
    bad.x_min.indices <- NULL
    bad.x_max.indices <- NULL
    if (any(x_min<0)){
      bad.x_min.indices <- (x_min<0)
      warning(ppaste("'x_min' must be positive in '",fn,"'"))
    } else {
      bad.x_max.indices <- (x_min>x_max)
      if (!any(is.na(bad.x_max.indices)) & any(bad.x_max.indices)){
        cat("(x_min,x_max):\n") ; print(data.frame("x_min"=x_min,"x_max"=x_max,"condition"=x_min>x_max)) 
        stop(ppaste("'x_max' must be larger than x_min in '",fn,"'")) 
      }
    }
    x_min[bad.x_min.indices] <- NA
    x_max[bad.x_max.indices] <- NA
  }
  
  return(list("k"=k,"x_min"=x_min,"x_max"=x_max))
}
truncpareto.check <- cmpfun(truncpareto.check)


"rtruncpareto" <- function(n, x_min, x_max, k, check=FALSE)
{
  if(check){
    # Perform default error checks:
    chk <- truncpareto.check(x_min=x_min, x_max=x_max, k=k, fn="rtruncpareto")
    k <- chk$k
    x_min <- chk$x_min
    x_max <- chk$x_max
  }
  # Note: works with x_max==Inf too :)
  return(x_min/((1+runif(n)*(((x_min/x_max)^k)-1))^(1/k)))
}
rtruncpareto <- cmpfun(rtruncpareto)

"dtruncpareto" <- function(x, x_min, x_max, k, log=FALSE,verbose=FALSE)
{
  if (length(log)>1)
    stop("Log must be TRUE/FALSE")
  
  # Perform default error checks:
  chk <- truncpareto.check(x_min=x_min, x_max=x_max, k=k, fn="dtruncpareto",verbose=verbose)
  k <- chk$k
  x_min <- chk$x_min
  x_max <- chk$x_max
  
  if(log){
    ret <- log(k)-(k+1)*log(x)-log((x_min^(-k))-(x_max^(-k)))
    # Override NA's with 0 or -Inf if parameter values were valid, but x is out of range:
    ret[(x<x_min | x>x_max) & !is.na(k) & !is.na(x_min) & !is.na(x_max)] <- -Inf
    
  } else {
    ret <- k/((x^(k+1))*((x_min^(-k))-(x_max^(-k))))
    
    ret[(x<x_min | x>x_max) & !is.na(k) & !is.na(x_min) & !is.na(x_max)] <- 0
  }
  
  return(ret)
}
dtruncpareto <- cmpfun(dtruncpareto)


"ptruncpareto" <- function(x, x_min, x_max, k, log=FALSE)
{
  if (length(log)>1)
    stop("Log must be TRUE/FALSE")
  
  # Perform default error checks:
  chk <- truncpareto.check(x_min=x_min, x_max=x_max, k=k, fn="ptruncpareto")
  k <- chk$k
  x_min <- chk$x_min
  x_max <- chk$x_max
  
  # Allow result to be of length x, of k, or of x_min. Vectors recycled.  
  if(log){
    ret <- log((x_min^(-k)-x^(-k))) - log((x_min^(-k)-x_max^(-k)))
    
    ret[ x<x_min & !is.na(k) & !is.na(x_min) & !is.na(x_max) ] <- -Inf
    ret[ x>x_max & !is.na(k) & !is.na(x_min) & !is.na(x_max) ] <- 0
  } else {
    ret <- (x_min^(-k)-x^(-k))/(x_min^(-k)-x_max^(-k))
    
    ret[ x<x_min & !is.na(k) & !is.na(x_min) & !is.na(x_max) ] <- 0
    ret[ x>x_max & !is.na(k) & !is.na(x_min) & !is.na(x_max) ] <- 1
  }
  
  return(ret) 
}
ptruncpareto <- cmpfun(ptruncpareto)

"qtruncpareto" <- function(q, x_min, x_max, k)
{
  # TODO: Note -- this requires x_min to be a vector of same length as q. We want it to work for both vector and scalar arguments!
  
  # Perform default error checks:
  chk <- truncpareto.check(x_min=x_min, x_max=x_max, k=k, fn="qtruncpareto")
  k <- chk$k
  x_min <- chk$x_min
  x_max <- chk$x_max
  
  ret <- (x_min^(-k) - q*(x_min^(-k)-x_max^(-k)))^(-1/k) #TODO: NOTE: Found ERROR here: value returned can be < x_min!!!
  
  # Any NA's in the parameter values will already produce an NA, so just handle q:
  ret[(q<0) | (q>1)] <- NA
  
  # Error fix: if q=0, then return strictly x_min
  ret[(q==0)] <- x_min[(q==0)]
  
  return(ret)
}
qtruncpareto <- cmpfun(qtruncpareto)


############################################################
##
## *brokenpareto functions:
##
## Functions for broken power laws with M break-points.
##
## Density, RNG, quantiles and CDF for a
## mixture of M truncated-Pareto's and one
## untruncated Pareto. Note that there are M+1 mixture
## components when there are M break-points:
##
## Y = p_0 X_0 + ... + p_{M-1} X_{M-1} + p_{M} X_{M}
##
## where:
##
## \sum_{j=0}^{M} p_{j} = 1.0
##
## p_{0} = 1-(bp_{1}/S_{min})^{-\theta_{0}}
## p_{j} = (1-(bp_{j+1}/bp_{j})^{-\theta_{j}})*(1-\sum_{i=0}^{j-1}p_{i})
## 
## and:
##
## X_{j} ~ Trunc-Pareto( \theta_{j} , [bp_{j} , bp_{j+1}] )
##
## where:
##
## for j=1,...,M-1, and:
##
## X_{M} ~ Pareto( \theta_{M} , bp_{M} ) .
##
##
############################################################

".brokenpareto.compute.p" <- function(m,x_min,k,bp,verbose=FALSE)
{
  ####
  ## Note: m is the number of mixture components -- NOT the number of break-points!!!!
  ## k is theta, and should be of length m
  ## bp is the vector of breakpoints, and should be length m-1
  ## x_min is the minimum value and must be a scalar
  ####
  p <- rep(NA,m)
  p[1] <- 1-((bp[1]/x_min)^(-k[1]))
  
  if (m>2){
    for (j in 2:(m-1)){
      p[j] <- (1-(bp[j]/bp[j-1])^(-k[j])) * (1-sum(p[1:(j-1)]))
    }
  }
  p[m] <- 1.0-sum(p[1:(m-1)])
  
  if (verbose){
    cat(paste("m = ",m,"\n",sep=""))
    cat(paste("x_min = ",x_min,"\n",sep=""))
    cat("bp:\n") ; print(bp)
    cat("k:\n") ; print(k)
    cat("p:\n") ; print(p)
  }
  
  if (any(p<0.0 | p>1.0))
    stop("Function error: computed 'p' out of range in '.brokenpareto.compute.p'")
  
  return(p)
}
.brokenpareto.compute.p <- cmpfun(.brokenpareto.compute.p)

# Not in NAMESPACE -- C-version does not need to interface with R

"rbrokenpareto" <- function(n, x_min, k, bp, verbose=FALSE)
{
  ####
  ## NOTE: rbrokenpareto cannot be vectorized!!!
  ####
  ## k is theta, and should be of length m
  ## bp is the vector of breakpoints, and should be length m-1
  ## x_min is the minimum value and must be a scalar
  ####
  
  if (is.null(bp)){
    # Divert back to regular Pareto:
    return(rpareto(n=n,x_min=x_min,k=k))
  }
  
  m <- length(k) # number of mixture components
  p <- rep(NA,m)
  if ((length(bp)+1) != m)
    stop("length of 'bp' must be length(k)-1")
  if (length(x_min)!=1)
    stop("'x_min' must be a scalar: 'rbrokenpareto' is not vectorizable")
  if (any(x_min>bp))
    stop("'x_min' must be less than all breakpoints")
  if (any(sort(bp)!=bp))
    stop("'bp' must be in increasing order")
  
  # Compute p's recursively:
  p <- .brokenpareto.compute.p(m=m,x_min=x_min,k=k,bp=bp,verbose=verbose)
  p.cumsum <- c(0,cumsum(p))
  
  # Decide which mixture component to use:
  u.idx <- findInterval(runif(n),p.cumsum)
  k.u  <- k[u.idx]
  bp <- c(x_min,bp,Inf)
  bpl.u <- bp[u.idx] 
  bpu.u <- bp[u.idx+1]
  
  # Generate from the corresponding mixture component:
  s <- rtruncpareto(n=n,x_min=bpl.u,x_max=bpu.u,k=k.u)
  
  if (verbose){
    cat("k:\n") ; print(k)
    cat("break-point(s):\n") ; print(bp)
    cat("p:\n") ; print(p)
    cat("mixture components chosen:\n") ; print(u.idx)
    cat("lower breakpoints:\n") ; print(bpl.u)
    cat("Upper breakpoints:\n") ; print(bpu.u)
    cat("Samples:\n") ; print(s)
    cat("p.cumsum:\n") ; print(p.cumsum)
  }
  
  return(s)
}
rbrokenpareto <- cmpfun(rbrokenpareto)

"dbrokenpareto" <- function(x, x_min, k, bp, log=FALSE, verbose=FALSE)
{
  if (is.null(bp)){
    # Divert back to regular Pareto:
    return(dpareto(x=x,x_min=x_min,k=k,log=log))
  }
  
  ## NOTE: dbrokenpareto cannot be vectorized in the parameters, only in x!!!
  m <- length(k) # number of mixture components
  if ((length(bp)+1) != m)
    stop("length of 'bp' must be length(k)-1")
  if (any(x_min<0)){
    warning("'x_min' must be a positive scalar in 'dbrokenpareto'.\n") 
    return(rep(NA,length(x)))
  }
  if (length(x_min)!=1)
    stop("'x_min' must be a scalar: 'dbrokenpareto' is not vectorizable")
  if (any(x_min>bp))
    stop("'x_min' must be less than all breakpoints")
  if (any(sort(bp)!=bp))
    stop("'bp' must be in increasing order")
  if (any(is.infinite(x) & x<0)) # x==-Inf
    stop("'x' cannot be -Inf in in 'drokenpareto'")
  
  # Separate x variable into two groups if some (x<x_min). x_min=scalar is required:
  bad.idx <- (x<x_min)
  if(verbose){
    if(any(bad.idx)){
      cat("bad.idx is actually BAD here!!!")
    }
  }
  n.old <- length(x)
  x.bad <- x[bad.idx]
  x     <- x[!bad.idx]   #NOTE: some values of x may be removed if some(x<x_min). 
  
  # Compute p's recursively:
  p <- .brokenpareto.compute.p(m=m,x_min=x_min,k=k,bp=bp,verbose=verbose)
  p.cumsum <- c(0,cumsum(p))
  
  # Decide which mixture component each x belongs to:
  n <- length(x)
  bp <- c(x_min,bp,Inf)
  
  u.idx <- findInterval(x,bp)
  
  k.u   <- k[u.idx]
  bpl.u <- bp[u.idx] 
  bpu.u <- bp[u.idx+1]
  p.u   <- p[u.idx]
  
  # Compute density from the corresponding mixture component:
  ret <- dtruncpareto(x=x,x_min=bpl.u,x_max=bpu.u,k=k.u,log=log)
  if (verbose){
    cat("step0: parameters sent to dtruncpareto:\n")
    print(data.frame(x=x,x_min=bpl.u,x_max=bpu.u,k=k.u,log=log))
    cat("step1: densities:\n") ; print(ret)
  }
  if (log){
    ret <- ret + log(p.u)
  } else {
    ret <- ret*p.u
  }
  
  # Merge back results from two groups if some (x<x_min), i.e. fill-in NA's to 'ret': 
  ret.final <- rep(NA,n.old)  
  ret.final[!bad.idx] <- ret 
  
  if (verbose){
    cat(ppaste("Computing density of brokenpareto on ",ifelse(log,"log","unlogged")," scale...\n"))
    cat("x:\n") ; print(x)
    cat("x.bad:\n") ; print(x.bad)
    cat("bad.idx:\n") ; print(bad.idx)
    cat("p:\n") ; print(p)
    cat("mixture component memberships:\n") ; print(u.idx)
    cat("mixture component probabilities:\n") ; print(p.u)
    cat("mixture-specific thetas:\n") ; print(k.u)
    cat("lower breakpoints:\n") ; print(bpl.u)
    cat("Upper breakpoints:\n") ; print(bpu.u)
    cat("p.cumsum:\n") ; print(p.cumsum)
    cat("densities:\n") ; print(ret)
    cat("final densities:\n") ; print(ret.final)
    cat("Finished computing 'dbrokenpareto'\n")
  }
  
  return(ret.final)
}
dbrokenpareto <- cmpfun(dbrokenpareto)


"pbrokenpareto" <- function(x, x_min, k, bp, log=FALSE, verbose=FALSE)
{
  if (is.null(bp)){
    # Divert back to regular Pareto:
    return(ppareto(x=x,x_min=x_min,k=k,log=log))
  }
  
  m <- length(k) # number of mixture components
  p <- rep(NA,m)
  if ((length(bp)+1) != m)
    stop("length of 'bp' must be length(k)-1")
  if (any(x_min<0)){
    warning("'x_min' must be a positive scalar in 'pbrokenpareto'.\n") 
    return(rep(NA,length(x)))
  }
  if (length(x_min)!=1)
    stop("'x_min' must be a scalar: 'pbrokenpareto' is not vectorizable in the parameters (only in x)")
  if (any(x_min>bp))
    stop("'x_min' must be less than all breakpoints")
  if (any(sort(bp)!=bp))
    stop("'bp' must be in increasing order")
  
  # Separate variables into two groups if some (x<x_min). x_min=scalar is required:
  bad.idx <- (x<x_min)
  n.old <- length(x)
  x.bad <- x[bad.idx]
  x     <- x[!bad.idx]   #NOTE: some values of x may be removed if (x<x_min). 
  
  # Compute p's recursively:
  p <- .brokenpareto.compute.p(m=m,x_min=x_min,k=k,bp=bp,verbose=verbose)
  p.cumsum <- c(0,cumsum(p))
  
  # Decide which mixture component each x belongs to:
  n <- length(x)
  bp <- c(x_min,bp,Inf)
  
  u.idx <- findInterval(x,bp)
  
  k.u  <- k[u.idx]
  bpl.u <- bp[u.idx] 
  bpu.u <- bp[u.idx+1]
  p.u <- p[u.idx]
  
  # Integrate over all previous mixture +
  # Integrate over the corresponding mixture component:
  ret <- p.cumsum[u.idx] + p.u*ptruncpareto(x=x,x_min=bpl.u,x_max=bpu.u,k=k.u)
  
  # Merge back results from two groups if some (x<x_min), i.e. fill-in NA's to 'ret': 
  ret.final <- rep(NA,n.old)  
  ret.final[!bad.idx] <- ret  
  
  if (log){
    ret.final <- log(ret.final)
  }
  
  if (verbose){
    cat("x:\n") ; print(x)
    cat("x.bad:\n") ; print(x.bad)
    cat("bad.idx:\n") ; print(bad.idx)
    cat("p:\n") ; print(p)
    cat("mixture component memberships:\n") ; print(u.idx)
    cat("mixture-specific thetas:\n") ; print(k.u)
    cat("lower breakpoints:\n") ; print(bpl.u)
    cat("Upper breakpoints:\n") ; print(bpu.u)
    cat("p.u:\n") ; print(p.u)
    cat("p.cumsum:\n") ; print(p.cumsum)
    cat("cdfs:\n") ; print(ret)
    cat("final cdf's (or on log-scale if log=TRUE):\n") ; print(ret.final)
  }
  
  return(ret.final)
}
pbrokenpareto <- cmpfun(pbrokenpareto)


"qbrokenpareto" <- function(q.bp, x_min, k, bp, verbose=FALSE)
{
  if (is.null(bp)){
    # Divert back to regular Pareto:
    return(qpareto(q=q.bp,x_min=x_min,k=k))
  }
  
  n <- length(q.bp)
  m <- length(k) # number of mixture components
  p <- rep(NA,m)
  if ((length(bp)+1) != m)
    stop("length of 'bp' must be length(k)-1")
  if (length(x_min)!=1)
    stop("'x_min' must be a scalar: 'rbrokenpareto' is not vectorizable")
  if (any(x_min>bp))
    stop("'x_min' must be less than all breakpoints")
  if (any(sort(bp)!=bp))
    stop("'bp' must be in increasing order")
  
  # Compute p's recursively:
  p <- .brokenpareto.compute.p(m=m,x_min=x_min,k=k,bp=bp,verbose=verbose)
  p.cumsum <- c(0,cumsum(p))
  
  # TODO: If q==1 return Inf
  # Make sure that 
  # Decide which mixture component to use:
  q.idx <- rep(NA,n)
  for (i in 1:n){
    q.idx[i] <- min(m,max(c(1:(m+1))[q.bp[i] >= p.cumsum])) #NOTE: min only needed to handle q==1 case
  }
  k.q  <- k[q.idx]
  bp <- c(x_min,bp,Inf)
  bpl.q <- bp[q.idx] 
  bpu.q <- bp[q.idx+1]
  p.idx <- p[q.idx]
  plo.idx <- p.cumsum[q.idx]
  
  # Look up scaled quantile from resulting truncated pareto:
  ret <- qtruncpareto(q=(q.bp-plo.idx)/p.idx, k=k.q, x_min=bpl.q, x_max=bpu.q)
  
  if (verbose){
    cat("q:\n") ; print(q.bp)
    cat("p:\n") ; print(p)
    cat("mixture-specific thetas:\n") ; print(k.q)
    cat("mixture component memberships:\n") ; print(q.idx)
    cat("lower breakpoints:\n") ; print(bpl.q)
    cat("Upper breakpoints:\n") ; print(bpu.q)
    cat("mixture-specific rescaled quantiles:\n") ; print((q.bp-plo.idx)/p.idx)
    cat("plo.idx:\n") ; print(plo.idx)
    cat("p.cumsum:\n") ; print(p.cumsum)
    cat("quantiles:\n") ; print(ret)
  }
  
  return(ret)
}
qbrokenpareto <- cmpfun(qbrokenpareto)


"brokenpareto.posterior" <- function(a,b,S.obs,S.mis,Smin,bp,verbose=FALSE)
{
  S <- c(S.obs,S.mis)
  m <- length(bp)+1  # number of mixture components
  bp <- c(Smin,bp,Inf)
  
  # Decide which mixture component each x belongs to:
  N <- length(S)
  S.idx <- sapply(1:N,function(x){
    min(m,max(c(1:(m+1))[S[x]>=bp])) #NOTE: min only useful if S == Inf (which shouldn't happen!)
  })
  
  # Pick out the lower and upper breakpoints for each data points:
  bpl.S <- bp[S.idx] 
  bpu.S <- bp[S.idx+1]
  
  # Count the number of sources for each mixture component:
  n <- rep(NA,m)
  
  # Compute sum_{i\in\mathcal{I}(j)} log(s_{i}/bp_{j}) for j=1,...,m:
  s3.log.ratio <- rep(NA,m)
  
  for (i in 1:m){
    tmp.idx <- (S.idx==i)
    n[i] <- sum(tmp.idx)
    s3.log.ratio[i] <- sum(log(S[tmp.idx]/bpl.S[tmp.idx]))
  }
  
  # Compute I{j~=m}*log(bp_{j+1}/bp_{j}) for j=1,...,m:
  bp.log.ratio <- c(log(bp[2:m]/bp[1:(m-1)]), 0)
  
  # Compute f(j) = I{j~=m}*\sum_{l=j+1}^{m} n(l) for j=1,...,m:
  n.cumsum <- c(sort(cumsum(n[m:2]),decreasing = TRUE), 0) # vector of partial sums
  
  # Compute partial posterior parameters of gamma distr of theta:
  a.post <- a + n
  b.post <- b + n.cumsum*bp.log.ratio + s3.log.ratio
  
  if (verbose>1){
    cat("a:\n") ; print(a)
    cat("b:\n") ; print(b)
    cat("S.obs:\n") ; print(S.obs)
    cat("S.mis:\n") ; print(S.mis)
    cat(paste("Smin = ",Smin,"\n",sep=""))
    cat("bp:\n") ; print(bp)
    cat("mixture component memberships:\n") ; print(S.idx)
    cat("lower breakpoints:\n") ; print(bpl.S)
    cat("Upper breakpoints:\n") ; print(bpu.S)
    cat("Mixture component sizes:\n") ; print(n)
    cat("Cumulative mixture component sums:\n") ; print(n.cumsum)
    cat("Breakpoint log-ratios:\n") ; print(bp.log.ratio)
    cat("log(s/bp) log-ratios:\n") ; print(s3.log.ratio)
    cat("a.post:\n") ; print(a.post)
    cat("b.post:\n") ; print(b.post)
  }
  
  return(list("a"=a.post,"b"=b.post))
}
brokenpareto.posterior <- cmpfun(brokenpareto.posterior)


############################################################
##
## *mixpareto functions:
##
## Functions for mixture of power laws with M different Smin-points.
##
## Density, RNG, quantiles and CDF for a
## mixture of M Pareto's. 
##
## Y = p_1 X_1 + ... + p_{M-1} X_{M-1} + p_{M} X_{M}
##
## where:
##
## \sum_{j=1}^{M} p_{j} = 1.0
##
## X_{j} ~ Pareto( \theta_{j} , Smin_{j} ) for j=1,...,M.
##
############################################################

# 
"rmixpareto" <- function(n, p, k, x_min, fixed.I.idx=NULL,verbose=FALSE)
{
  if(verbose){
    cat("rmixpareto()...\n n:\n"); print(n)
    cat("s_min:\n"); print(x_min)
    cat("slope:\n"); print(k)
    cat("p:\n"); print(p)
    cat("fixed.I.idx:\n"); print(fixed.I.idx)
    
  }
  
  if (n<=0){
    return(numeric(0))
  }
  
  if (any(k<0)){
    k[k<0] <- NA
    warning("'k' must be positive.\n") 
  }
  if (any(x_min<0)){
    x_min[x_min<0] <- NA
    warning("'x_min' must be positive.\n") 
  }
  
  if (!is.null(p)){
    
    # Generate mixture indicators:
    if (is.null(fixed.I.idx)){
      # Generate all new elements of I.idx
      I <- rmultinom(n=n,size=1,prob=p)
      I.idx <- apply(I==1,2,which)  
      
    } else {    
      I.idx <- fixed.I.idx
      
      #Sometimes, if N.t > N.t-1, a number of last elements of fixed.I.idx are NA. They must be generated.
      nai <- is.na(fixed.I.idx) 
      if (any(nai)){
        I <- rmultinom(n=sum(nai),size=1,prob=p)
        I.idx[nai] <- apply(I==1,2,which)  
      }
    }
    ###Problem is here:   n is sometimes different from length(I.idx) ....
    
    # Generate fluxes from mixture of Pareto's:
    S <- rpareto(n=n,x_min=x_min[I.idx],k=k[I.idx],verbose=verbose)
    
  } else {
    
    # Compatibility code for regular Pareto distribution:
    I.idx <- NULL  
    
    # Generate fluxes from Pareto's:
    S <- rpareto(n=n,x_min=x_min,k=k,verbose=verbose)
    
  }
  
  return(list("S"=S,"I.idx"=I.idx))
}
rmixpareto <- cmpfun(rmixpareto)

"dmixpareto" <- function(x, p, x_min, k, log=FALSE, verbose=FALSE)
{
  # 'dmixpareto' cannot be vectorized. length(p)=length(k)=length(x_min)
  
  m <- length(p)
  n <- length(x)
  
  if (length(log)>1)
    stop("Log must be TRUE/FALSE")
  if (any(k<0)){
    k[k<0] <- NA
    warning("'k' must be positive.\n") 
  }
  if (any(x_min<0)){
    x_min[x_min<0] <- NA
    warning("'x_min' must be positive.\n") 
  }
  
  if ( (length(x_min) != m) | (length(k) != m) )
    stop("length of parameters must be equal in 'dmixpareto': length(p)=length(k)=length(x_min).")
  
  # Loop over vector of values in x:
  ret <- numeric(n)
  for (j in 1:n) {
    # Compute (weighted) probability of membership for each mixture component, and add them up:
    #    ret[j] <- sum( sapply(1:m, function(i){p[i]*dpareto(x[j], x_min[i], k[i], log=FALSE)}) )
    tmp <- numeric(m)
    for(i in 1:m){  
      tmp[i] = p[i]*dpareto(x[j], x_min[i], k[i], log=FALSE) 
    }
    ret[j] <- sum(tmp)  
  }
  
  if (log){
    ret <- log(ret)
  }
  
  if (verbose){
    cat("density:\n") ; print(ret)
  }
  
  return(ret) 
}
dmixpareto <- cmpfun(dmixpareto)


"pmixpareto" <- function(x, p, x_min, k, log=FALSE, verbose=FALSE)
{
  # 'pmixpareto' cannot be vectorized. length(p)=length(k)=length(x_min)
  
  m <- length(p)
  n <- length(x)
  
  if ( (length(x_min) != m) | (length(k) != m) )
    stop("length of parameters must be equal in 'pmixpareto': length(p)=length(k)=length(x_min).")
  if (length(log)>1)
    stop("Log must be TRUE/FALSE")
  if (any(k<0)){
    k[k<0] <- NA
    warning("'k' must be positive.\n") 
  }
  if (any(x_min<0)){
    x_min[x_min<0] <- NA
    warning("'x_min' must be positive.\n") 
  }
  
  # Loop over vector of values in x:
  ret <- numeric(n)
  for (j in 1:n) {
    # Compute (weighted) probability of membership for each mixture component, and add them up:
    #ret[j] <- sum( sapply(1:m, function(i){p[i]*ppareto(x[j], x_min[i], k[i], log=FALSE)}) )
    ret[j] <- sum( unlist(lapply(1:m, function(i){p[i]*ppareto(x[j], x_min[i], k[i], check.args=FALSE, log=FALSE)})) )
  }
  
  if (log){
    ret <- log(ret)
  }
  if (verbose){
    cat("cdf of mix-pareto:\n") ; print(ret)
  }
  
  return(ret)
}  
pmixpareto <- cmpfun(pmixpareto)
####################################################################################################
# VERY SLOW FUNCTION!!!!
# Write C Code to replace this..
# 
"qmixpareto" <- function(q, p, x_min, k, tol=1e-5, n.tries.max=2000, use.C=TRUE, verbose=FALSE)
{ 
  # Invalid arguments should return an error before this (new C code does this):
  if (any(x_min<0)){
    stop("'x_min' must be positive in 'qmixpareto'")
  }
  if (any(k<0)){
    stop("'k' must be positive in 'qmixpareto'")
  }
  
  if (use.C){
    input_args <- list("q"=as.double(q),
                       "p"=as.double(p),
                       "x_min"=as.double(x_min),
                       "k"=as.double(k),
                       "tol"=as.double(tol),
                       "n.tries.max"=as.integer(n.tries.max),
                       "verbose.level"=as.integer(verbose*1))
    ret <- .Call("qmixpareto_C", input_args, package = "logNlogS")    #tries to loop up R package
    #    ret <- .Call("qmixpareto_C", input_args, PACKAGE = "logNlogS")   #tries to look up dynamix library
    
  } else {
    
    # Use R code... 
    
    use.bisection <- TRUE
    #  use.newton <- FALSE #TRUE
    
    tol <- tol#*min(x_min)   #resize the tolerance to be on the same scale as x_min
    n <- length(q)
    m <- length(p) # number of mixture components
    if ( (length(x_min) != m) | (length(k) != m) )
      stop("length of parameters must be equal in 'pmixpareto': length(p)=length(k)=length(x_min).") 
    
    if (any(x_min<0)) { warning("'x_min' must be positive.\n") }
    if (any(k<0))     { warning("'k' must be positive.\n") }
    
    ret <- numeric(n)
    
    if(use.bisection){
      # Use Bisection Algorithm to find quantile quickly
      if (verbose){
        cat("Bisection method to get quantiles of mix-pareto:\n")
        cat("q: \n") ; print(q)
        cat("p: \n") ; print(p)
        cat("x_min: \n") ; print(x_min)
        cat("k: \n") ; print(k)
        cat("n.tries.max = ",n.tries.max,"\n")
      }
      
      for (j in 1:n){ #for each element of q...
        n.tries     <- 0
        start.q     <- q[j]
        a <- min( x_min / ((1-start.q)^(1/k)) )
        b <- max( x_min / ((1-start.q)^(1/k)) ) 
        
        if (verbose){
          cat("Bisection method: Lower and Upper boundaries has been selected:\n")
          cat("a = \n"); print(a)    
          cat("b = \n"); print(b)
          cat("q[j] = \n"); print(q[j])
        }
        
        # Check upper boundary  
        if (is.infinite(b) && b>0){
          ret[j] <- Inf
          next
        }
        
        c <- a
        G.a <- 1-sum(p*(x_min/a)^(k)) - q[j]
        G.c <- G.a
        
        while (n.tries<=n.tries.max){ # limit iterations to prevent infinite loop
          if (verbose){
            cat(c(" try ",n.tries,": G(",c,") = ",G.c,"\n"))
          }
          if (abs(G.c)<tol){# | (b-a)/2 < tol)  # solution found
            sol <- c
            break
          }
          
          # Bisection Algorithm:
          n.tries <- n.tries + 1  # increment step counter
          c   <- (a+b)/2  # new midpoint
          G.c  <- 1-sum(p*(x_min/c)^(k)) - q[j]
          
          if ( sign(G.c)==sign(G.a) ){
            a <- c  # lower boundary moves up
            G.a <- G.c
          } else { 
            b <- c  # upper boundary moves down
          }
          
        } #END besection while-loop
        
        if (verbose){
          cat("bisection n.tries = ") ; cat(n.tries); cat("\n")
          cat("cdf of mix-pareto F(sol):\n") ; print(G.c+q[j])
          cat("current solution is sol = \n"); print(sol)
        }
        
        if(n.tries>=n.tries.max){
          stop("Bisection method failed in qmixpareto: max number of steps exceeded")
        }
        
        ret[j] <- sol # update solution of bisection
        
      } # END for-loop looping over each element of q
      
      #   } else if(use.newton){
      #         # Use Newton Method to find quantile quickly
      # ###PRINT
      # #    if(TRUE)
      #     if (verbose){
      #       cat("Newton method to get quantiles of mix-pareto:\n")
      #       cat("q: \n") ; print(q)
      #       cat("p: \n") ; print(p)
      #       cat("x_min: \n") ; print(x_min)
      #       cat("k: \n") ; print(k)
      #       cat("n: \n") ; print(n)
      #       cat("m: \n") ; print(m)
      #       cat("n.tries.max = ",n.tries.max,"\n")
      #     }
      #     
      #     for (j in 1:n){ #for each element of q...
      #       n.tries     <- 0
      #       a <- min(x_min)
      #       b <- max(qpareto(q=0.99999,x_min=x_min,k=k))
      #       c <- a + q[j]*(b-a) #good starting guess
      #       if(verbose){
      #         cat("a: \n") ; print(a)
      #         cat("b: \n") ; print(b)
      #       }
      #       
      #       while (n.tries<n.tries.max){ # limit iterations to prevent infinite loop
      #         # Newton Algorithm:
      #         G.c  <- 1-sum(p*(x_min/c)^(k)) - q[j]
      #         dG.c <- sum(p*(x_min^k)*k*(c^(-k-1)))
      #         
      #         if(verbose){
      #           cat("c: \n") ; print(c)
      #           cat("G.c: \n") ; print(G.c)
      #           cat("dG.c: \n") ; print(dG.c)
      #         }
      #         
      #         if (abs(G.c)<tol){# | abs(G.c/dG.c) < tol)  # solution found
      #           sol <- c
      #           break
      #         }
      #         n.tries <- n.tries + 1  # increment step counter
      #         ###NOTE:  A problem exists here if c is negative!!!
      #         c <- c - G.c/dG.c  # new interval
      #         
      #         if(c<0){
      #           stop("c<0 ! Newton's method fails. :(' ")
      #         }
      #       } #END Newton while-loop
      # 
      #       if (verbose){
      #         cat("newton n.tries = ") ; cat(n.tries); cat("\n")
      #         cat("cdf of mix-pareto F(sol):\n") ; print(G.c+q[j])
      #         cat("current solution is sol = \n"); print(sol)
      #       }
      # 
      #       if(n.tries>=n.tries.max){
      #         stop("Newton method failed in qmixpareto: max number of steps exceeded")
      #       }
      # 
      #       ret[j] <- sol # update solution of bisection
      #       
      #     } # END for-loop looping over each element of q
      
    } else {  # Use cdf approach, evaluated on the grid
      # Evaluate cdf everywhere / for all range of x
      N <- 1000
      max.id <- which.max(x_min)
      x <- seq(from=min(x_min), to=qpareto(q=0.99999,x_min=x_min[max.id],k=k[max.id]), length.out=N)
      cdf <- pmixpareto(x=x, p=p, x_min=x_min, k=k, verbose=verbose)
      
      # Numerically draw from cdf
      for (j in 1:n){
        ret[j] <- cdf[ sum(cdf <= q[j]) ]
      }
    }
  } # END if (use.C){ ... } else { ... }
  
  # If individual q's are out of bounds give NA's (new C code does this):
  ret[(q<0) | (q>1)] <- NA
  
  if (verbose){
    cat("quantiles:\n") ; print(ret) 
  }
  
  return(ret)  
}
qmixpareto <- cmpfun(qmixpareto)


####################################################################################################

"mixpareto.posterior" <- function(a,b,S.obs,S.mis,Smin,I.idx,verbose=FALSE)
{
  
  n <- length(S.obs)
  N <- length(c(S.obs,S.mis))
  m <- length(a)
  
  # retreive 0's and 1's indices of mixture components
  I <- matrix(numeric(0), nrow=m, ncol=N)
  n.j <- numeric(m)
  
  for(j in 1:m){
    tfvec   <- I.idx == j
    I[j,]   <- tfvec
    n.j[j] <- sum(tfvec)
  }
  #  I <- t(sapply(1:max(I.idx),function(i){I.idx==i})) 
  #  n.j <- table(I.idx)
  
  if(n<N){
    a.post <- a + n.j
    b.post <- b + as.numeric(I[,1:n,drop=FALSE] %*% (log(S.obs/Smin[I.idx[1:n]])) + I[,(n+1):N,drop=FALSE] %*% (log(S.mis/Smin[I.idx[(n+1):N]])))
    
  } else { #n==N
    a.post <- a + n.j
    b.post <- b + as.numeric(I[,1:n,drop=FALSE] %*% (log(S.obs/Smin[I.idx[1:n]]))) 
  }
  
  if (verbose){
    cat("a:\n") ; print(a)
    cat("b:\n") ; print(b)
    cat("a.post:\n") ; print(a.post)
    cat("b.post:\n") ; print(b.post)
  }
  
  return(list("a"=a.post,"b"=b.post))
}

mixpareto.posterior <- cmpfun(mixpareto.posterior)