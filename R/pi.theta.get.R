"pi.theta.get" <- function(theta, Smin, bp, gamma, g, E,
                           sigma=0,
                           mc.integral=TRUE, nsamples=10000,
                           truncate=TRUE, verbose=FALSE){
  
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
  #
  # Output: p = marginal probability of detection, function of theta only
  ####################################################################################################
    
  m <- length(theta)
 
  # evaluate pi with error
  pi <- g(lambda=Smin*E/gamma) + rnorm(1,mean=0,sd=sigma)  # NOTE: this g==CONST
  
#     # Use Monte Carlo integration 
#     tmp <- pi.theta.eval.mc(nsamples=nsamples, theta.grid=theta, Smin.grid=Smin, bp.grid=bp, p.t.grid=p.t, 
#                             use.bp=TRUE, use.mix=FALSE,
#                             m=m, gamma=gamma, pble=pble, g=g, g.type=g.type, verbose=verbose)    
#   
#   pi <- tmp$pi
  
  if (truncate) {
    if(verbose)
      cat('Truncating pi to be in [0,1]...\n')
    pi <- pmin(pi,1)
    pi <- pmax(pi,0) 
  }
  
  return(pi)
  
}

pi.theta.get <- cmpfun(pi.theta.get)
