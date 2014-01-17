"g.compute" <- function(lambda,bg,E,L,g.type, min.photon.limit=1){
  
  ###################################################################################
  # g.compute   function computes the probability of observing source, 
  #             conditional on flux of source through lambda(S)=S*E/gamma.
  #             It is similar to g.func described in the test functions
  #
  ##### NOTE: g MUST BE VECTORIZED!!!!!!
  #
  # Input: lambda  = intensity of observed number of source photons
  #        bg      = background ct/pix of observed sources
  #        E       = exposure for observed sources
  #        L       = location for observed sources
  #        g.type  = type of g-function {"step","smooth","table"}
  #        min.photon.limit = assumption: fixed value of the minimum photon limit allowed to be observed
  #
  # Output: g  = (function of lambda(S),B,L,E) probability of observing source
  #################################################################################
  
  obs.intensity <- lambda
  
  if (is.character(g.type)){
    
    if (g.type == "step"){
      ret <- ifelse(obs.intensity > min.photon.limit, 1,0)
    } else if (g.type == "smooth"){
      
      use.C <- TRUE
      
      #Put this into C code
      if(use.C){
        
        input_args <- as.numeric(obs.intensity)
        ret <- .Call("g_smooth_C", input_args, package = "logNlogS")    #tries to loop up R package
        
      } else {
        
        #Incompleteness curve based on LOW background in: http://iopscience.iop.org/0004-637X/661/1/135/fulltext    (Table 1)
        a0 <- 11.12
        a1 <- -0.83
        a2 <- -0.43
        ret <- 1.0 - a0*(obs.intensity^a1)*exp(a2*obs.intensity)
        ret[ret<=0] <- 0  #Truncate negative values!
        
      } #Use C-code
      
    } else {
      stop("g.table has not been properly read.")
    }
    
  } else { 
    # Case: g.type == g.table    
    # Use binary search to quickly find index, and interpolate g
    # Find closest values of grid to B,L,E using binary search
    i   <- findInterval(L,g.type$pars$L)
    j   <- findInterval(bg,g.type$pars$bg)
    m   <- findInterval(obs.intensity,g.type$pars$lambda)
    if (any(c(i,j,m)<=0,na.rm=TRUE)) {
      stop("Index in findInterval from g.type$table is outside of L, k, or lambda plausible values.")
    }
    length.lambda <- g.type$pars$length.lambda
    length.bg     <- g.type$pars$length.bg
    idx <- 1+(m-1) + (j-1)*length.lambda + (i-1)*length.bg*length.lambda    
    ret <- g.type$table$g[idx] 
    # Linear interpolation along lambda (src ct) direction only (Note: we do not perform trilinear interpolation)
    ret <- ret + (obs.intensity-g.type$pars$lambda[m]) * (g.type$table$g[idx+1]-g.type$table$g[idx]) / (g.type$pars$dlambda)
    # Truncate to 1
    ret[ret>1] <- 1.0
  }
  return(ret)
}

g.compute <- cmpfun(g.compute)
