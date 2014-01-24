"tune.v" <- function(v.th, v.sm, v.bp, v.so, prop.ct.theta, prop.ct.Smin, prop.ct.bp, prop.ct.S.obs, prop.state, 
                     tune.iter, amt=sqrt(5), verbose=FALSE){
  
  ####################################################################################################
  # tune.v   Reduce or increase tuning SD parameters according for acceptance of proposals to be within 20-60%
  #
  # Input: v.th           = value SD, tuning parameter for slope theta to accept 20-60% of proposals
  #        v.sm           = value SD, tuning parameter for Smin to accept 20-60% of proposals
  #        v.bp           = value SD, tuning parameter for eta/bp to accept 20-60% of proposals
  #        v.so           = value SD, tuning parameter for observed fluxes S.obs to accept 20-60% of proposals
  #        prop.ct.theta = number of accepted proposals for theta
  #        prop.ct.Smin  = number of accepted proposals for Smin
  #        prop.ct.bp    = number of accepted proposals for bp
  #        prop.ct.S.obs = number of accepted proposals for S.obs
  #        prop.state    = default state of number of accepted proposals for S.com
  #        tune.iter     = chain iteration at which to tune
  #        amt           = amount factor by which to tune
  #        verbose       = (T/F) display progress of program
  #
  # Output: v.th        = SD vector, tuning parameter to accept 20-60% of proposals of theta
  #         v.sm        = SD value, tuning parameter to accept 20-60% of proposals of Smin
  #         v.bp        = SD vector, tuning parameter to accept 20-60% of proposals of bp
  #         v.so        = SD vector, tuning parameter to accept 20-60% of proposals of S.obs
  #         prop.state = storage count of accepted proposals at last tuning state
  ####################################################################################################
  #verbose <- TRUE
  
  if (verbose) {
    cat("_____Tune v parameters_____\n")
    cat("Before tuning:\n")
    cat("v.so\n"); print(v.so)
    cat("v.th :theta\n"); print(v.th)
    cat("v.sm :smin\n"); print(v.sm)
    cat("v.bp :bp\n"); print(v.bp)
    cat("prop.state\n"); print(prop.state)
  }
  
  # Retreive proposal state of last <tune.iter> iterations
  prop.state.t <- prop.ct.theta - prop.state$prop.ct.theta
  prop.state.m <- prop.ct.Smin  - prop.state$prop.ct.Smin
  prop.state.e <- prop.ct.bp  - prop.state$prop.ct.bp
  prop.state.s <- prop.ct.S.obs  - prop.state$prop.ct.S.obs
  # Evaluate proportion of accepted proposals
  tpa <- round(prop.state.t/tune.iter,3)*100
  mpa <- round(prop.state.m/tune.iter,3)*100
  epa <- round(prop.state.e/tune.iter,3)*100
  spa <- round(prop.state.s/tune.iter,3)*100
  # Tune 
  v.th[tpa<20] <- v.th[tpa<20]/amt      #reject too many => reduce SD
  v.th[tpa>60] <- v.th[tpa>60]*amt      #accept too many => increase SD
  v.sm[mpa<20] <- v.sm[mpa<20]/amt
  v.sm[mpa>60] <- v.sm[mpa>60]*amt
  v.bp[epa<20] <- v.bp[epa<20]/amt
  v.bp[epa>60] <- v.bp[epa>60]*amt
  v.so[spa<20] <- v.so[spa<20]/amt
  v.so[spa>60] <- v.so[spa>60]*amt
  # Update current proposal state
  prop.state <- list("prop.ct.theta"=prop.ct.theta, "prop.ct.Smin"=prop.ct.Smin, 
                     "prop.ct.bp"=prop.ct.bp, "prop.ct.S.obs"=prop.ct.S.obs) #update current proposal state
  
  if (verbose) {
    cat("After tuning:\n")
    cat("spa\n"); print(spa)
    cat("tpa\n"); print(tpa)
    cat("mpa\n"); print(mpa)
    cat("epa\n"); print(epa)
    cat("v.so\n"); print(v.so)
    cat("v.th:theta\n"); print(v.th)
    cat("v.sm:smin\n"); print(v.sm)
    cat("v.bp:bp\n"); print(v.bp)
    cat("prop.state, after tuning:\n"); print(prop.state)
  }
  
  return(list("v.th"=v.th,"v.sm"=v.sm,"v.bp"=v.bp,"v.so"=v.so,"prop.state"=prop.state))
}

tune.v <- cmpfun(tune.v)
