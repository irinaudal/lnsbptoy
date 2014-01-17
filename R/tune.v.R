"tune.v" <- function(v.S, v.theta, v.bp, v.sm, prop.ct.S.obs, prop.ct.theta, prop.ct.bp, prop.ct.Smin, prop.state, 
                     tune.iter, N.t, n, tune.const=5, verbose=FALSE){
  
  ####################################################################################################
  # tune.v   Reduce or increase tuning SD parameters according for acceptance of proposals for S.obs, theta, bp, Smin to be within 20-60%
  #
  # Input: v.S     = value SD, tuning parameter for observed fluxes S.obs to accept 20-60% of proposals
  #        v.theta = value SD, tuning parameter for slope theta to accept 20-60% of proposals
  #        v.bp    = value SD, tuning parameter for bp to accept 20-60% of proposals
  #        v.sm    = value SD, tuning parameter for Smin to accept 20-60% of proposals
  #        N.t     = total number of sources, N.com
  #        n       = number of observed sources
  #        prop.ct.S.obs = number of accepted proposals for S.obs
  #        prop.ct.theta = number of accepted proposals for theta
  #        prop.ct.bp    = number of accepted proposals for bp
  #        prop.ct.Smin  = number of accepted proposals for Smin
  #        prop.state    = default state of number of accepted proposals for S.com
  #        tune.iter     = chain iteration at which to tune
  #        verbose       = (T/F) display progress of program
  #
  # Output: v.S     = SD vector, tuning parameter to accept 20-60% of proposals of S.obs
  #         v.theta = SD vector, tuning parameter to accept 20-60% of proposals of theta
  #         v.bp    = SD vector, tuning parameter to accept 20-60% of proposals of bp
  #         v.sm    = SD value, tuning parameter to accept 20-60% of proposals of Smin
  #         prop.state = storage count of accepted proposals at last tuning state
  ####################################################################################################
  #verbose <- TRUE
  
  if (verbose) {
    cat("_____Tune v parameters_____\n") 
  } 
  
  m <- length(v.theta)
  prop.state.s <- prop.ct.S.obs - prop.state$prop.ct.S.obs     #retreive proposal state of last 200 iterations
  prop.state.t <- prop.ct.theta - prop.state$prop.ct.theta
  prop.state.b <- prop.ct.bp    - prop.state$prop.ct.bp
  prop.state.m <- prop.ct.Smin  - prop.state$prop.ct.Smin
  spa <- round(prop.state.s/tune.iter,3)*100       #proportion of accepted proposals
  tpa <- round(prop.state.t/tune.iter,3)*100      
  bpa <- round(prop.state.b/tune.iter,3)*100      
  mpa <- round(prop.state.m/tune.iter,3)*100      
  v.S[spa<20] <- v.S[spa<20]/tune.const        #reject too many => reduce SD
  v.S[spa>60] <- v.S[spa>60]*tune.const        #accept too many => increase SD     
  v.theta[tpa<20] <- v.theta[tpa<20]/tune.const
  v.theta[tpa>60] <- v.theta[tpa>60]*tune.const
  v.bp[bpa<20] <- v.bp[bpa<20]/tune.const
  v.bp[bpa>60] <- v.bp[bpa>60]*tune.const
  v.sm[mpa<20] <- v.sm[mpa<20]/tune.const
  v.sm[mpa>60] <- v.sm[mpa>60]*tune.const
  prop.state <- list("prop.ct.S.obs"=prop.ct.S.obs,"prop.ct.theta"=prop.ct.theta,
                     "prop.ct.Smin"=prop.ct.Smin,  "prop.ct.bp"=prop.ct.bp) #update current proposal state
  
  if(verbose) {
    cat("prop.state.s\n"); print(prop.state.s)
    cat("prop.state.t\n"); print(prop.state.t)
    cat("prop.state.b\n"); print(prop.state.b)
    cat("prop.state.m\n"); print(prop.state.m)
    cat("prop.state\n"); print(prop.state)
    cat("n\n"); print(n)
    cat("spa\n"); print(spa)
    cat("tpa\n"); print(tpa)
    cat("bpa\n"); print(bpa)
    cat("mpa\n"); print(mpa)
    cat("v.s\n"); print(v.S)
    cat("v.theta\n"); print(v.theta)
    cat("v.bp\n"); print(v.bp)
    cat("v.sm\n"); print(v.sm)
  }
  return(list("v.S"=v.S,"v.theta"=v.theta,"v.bp"=v.bp,"v.sm"=v.sm,"prop.state"=prop.state))
}

tune.v <- cmpfun(tune.v)
