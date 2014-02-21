
############################################################################################################
##
## CREATE pi(theta,Smin) by parts, then interpolate all components.
##
##
############################################################################################################


# Model: 
# Y_i ~ Poisson(lambda = S_i*E/gamma)
# S_i ~ Pareto(theta,Smin,bp)
# theta ~ gamma(a,b)
# N ~ negbinom(alpha,beta)
# Smin ~ gamma(am,bm)
# bp ~ fcn of eta|Smin
# eta | Smin ~ N(mu,C), where eta = log(bp-Smin)
# I_i ~ Bin(1,g(lambda_i))
# g(lambda_i = S_i*E/gamma) = P(observe i-th source) = const
# pi(theta,Smin,bp) = \int g(lambda_i) * p(S_i | theta,Smin,bp)  dS_i 
#    \approx  g+e,   where e ~ N(0, sigma)
# OR
#    \approx via MC integration <- TODO: implement this later after the code validates
# sigma = 0.001, say is error that is incurred while approximating pi


# Start from clean workspace:
rm(list=ls())

print(Sys.info()["sysname"])

# Setup correct paths, R history etc...
if (Sys.info()["sysname"]=="Darwin"){
  # Note: main.dir should NOT have a trailing /
  main.dir   <- "~/Dropbox/Log_N_Log_S/R_Code"
  setwd(main.dir)
  libsrc.dir <- "."
  library(lnsbptoy)
  
} else if (Sys.info()["sysname"]=="Linux"){
  libdir <- "/home/isudal/R/x86_64-pc-linux-gnu-library/3.0/"
  library(lnsbptoy,lib.loc=libdir)
  main.dir   <- "/home/isudal/bplns" #lns_data/"
  setwd(main.dir)
  libsrc.dir <- "."
}

## Handle batch job arguments:
args <- commandArgs(TRUE)  # 1-indexed version is used now.
#args <- 4
cat(paste("Command-line arguments:\n"))
print(args)

###################
sim_start <- 0 #1000  #HERE
if (length(args)==0){
  sinkit <- TRUE #FALSE
  sim_num <- sim_start + 1
} else {
  sinkit <- TRUE
  sim_num <- sim_start + as.numeric(args[1])
}
set.seed(762*sim_num + 1330931 + 10231*92) 
###################

# Set simulation parameters

################################
# Simulation specs
###################
verbose             <- FALSE #TRUE
precompute.pi.theta <- TRUE #FALSE
theta.grid.eps      <- 0.001 # How extreme a grid to compute for theta? 
just.pi.theta       <- TRUE#FALSE #TRUE
true.model <- model <- "bp"    ## Use load.res.image=TRUE to use true model
do.pi.plot          <- FALSE #ifelse(sim_num==1,TRUE,FALSE) # plot of pi(theta,Smin)
truncate            <- TRUE
M <- 1  ## length of subsection of grid of Smin
interpolate.pi       <- FALSE #TRUE
do.profile          <- FALSE#TRUE # FALSE
remove.main.grid    <- TRUE
######################################
data_num <- sim_num  #1-100
job_num  <- 1001 #1006
######################################
# # Define mean and sd parameters for g.func in norm-cdf
# # Define mean and sd parameters for g.func in norm-cdf
# mu.vec <- c(7,13,-12,60) #c(-50,seq(-20,15,by=5),70)
# sd.vec <- c(3,3,15,30) #10
# 
# mu <- mu.vec[job_num-1100]
# s <- sd.vec[job_num-1100]

############################################################################
################### Vary parameters stage ####################################
############################################################################  
model <- true.model
fixed.fit.Smin <- FALSE
fixed.Smin <- fixed.fit.Smin
fixed.bp <- FALSE  #either vector of Value/TRUE/FALSE/NULL

#cat("Full model vector:\n")
#print(data.frame("model"=model.vec,"fixed.fit.Smin"=fixed.fit.Smin.vec,"fixed.bp"=fixed.bp.vec))

cat("\n")
cat("====================================\n")
cat("Model to be fit is as follows...\n\n")
cat("model:\n")
print(model)
cat("fixed.fit.Smin:\n")
print(fixed.fit.Smin)
cat("fixed.Smin:\n")
print(fixed.Smin)
cat("fixed.bp:\n")
print(fixed.bp)
cat("====================================\n")
cat("\n")

model.old <- model

##################################

cat(ppaste("Beginning analysis of dataset number ",data_num,"\n"))
cat(ppaste("Beginning analyzing of job number (fitted model) ",job_num,"\n"))


###################
# Simulation specs
###################
set.seed(762*data_num  - 56*job_num + 1330931)
###################

fixed.theta <- FALSE #c(.6,1.2)   # Vector of values / FALSE
fixed.N     <- FALSE #150 #FALSE  # Value / FALSE
fixed.Smin  <- FALSE #1*10^-13
fixed.bp    <- FALSE #4*10^-13 #NULL
fixed.S.obs <- FALSE
fixed.S.mis <- FALSE

outer.dir.extra     <- paste("/bp_lns_toy_example_1bp_gpnorm_prec", sep="")
output.dir.extra    <- paste("/dataset_",sim_num,sep="")
output.dir.outer    <- paste(main.dir,outer.dir.extra,sep="")
output.dir          <- paste(output.dir.outer,output.dir.extra,sep="")
output.dir.wo.main  <- paste(outer.dir.extra,output.dir.extra,sep="")

sink.file           <- file.path(output.dir,"foo.txt")
profile_file        <- paste(output.dir,"/basicProfile.Rprof",sep="")
save.res.file       <- paste(output.dir,"/res_image.RData",sep="")
save.image.file     <- paste(output.dir,"/store_image.RData",sep="")

pi.file.location    <- ppaste(output.dir.outer,"/pi_theta_",job_num,"_",sim_num,".RData")  
output.pi.plot     <- ppaste(output.dir,"/pi.function_",sim_num,".jpg")


####################################################

if (!file.exists(output.dir.outer)){   #set directory - done after simulate()
  dir.create(output.dir.outer, recursive=TRUE)
  #   setwd(output.dir)
}
if (!file.exists(output.dir)){   #set directory - done after simulate()
  dir.create(output.dir, recursive=TRUE)
  #   setwd(output.dir)
}

if (sinkit) {
  cat("Sinking to...\n")
  print(sink.file)
  sink(sink.file)
}
cat(ppaste("Command-line arguments:\n"))
print(sim_num)

# Incompleteness
"g.func" <- function(lambda) {
  #return(rep(0.8,length(lambda)))   # constant incompleteness
  return(pnorm(lambda,mean=-50,sd=100))   #mean=-25,sd=50))  # 1 minus exponential decay
}
nsamples  <- 10000
E <- 400000 #190000
gamma <- 1.6*(10^(-9))

# Precomputing pi parameters
length.theta <- 100
length.Smin  <- 200
length.bp    <- 150

# Parameters for Normal priors for break-point(s):
mu <- -28.8 # c(-30) # dim: eta(bp)=m-1
C  <- 0.2 # c(0.5) 
##############################
# Parameters for Gamma prior(s) for theta(s):
a <- c(100,200) # c(42,80) #100  #c(10,30)   #c(12,10) #c(25,30)
b <- c( 80, 80)  # c(70,80)  #80   #c(20,30)   #c(18,10) #c(40,30)
##############################
# Parameters for Gamma prior(s) for minimum threshold parameter Smin(s):
Smin.prior.mean <- 1.0*10^(-13)    # Mean = am * bm
Smin.prior.var  <- (4.5*10^(-14))^2  # Var  = am * bm^2
am  <- Smin.prior.mean^2 / Smin.prior.var   # shape
bm  <- am / Smin.prior.mean                 # rate
##############################
# parameters for neg binom prior of N
#    Mean = alpha / beta
#    Var  = alpha * (1+beta) / beta^2
#    => beta   = prior.mean / (prior.var - prior.mean)
#    => alpha = beta * prior.mean  
N.prior.mean <- 300 #150 #80 ##300
N.prior.var  <- 50^2  #(1*100)^2
beta  <- N.prior.mean / (N.prior.var - N.prior.mean)
alpha <- beta * N.prior.mean
###################
# Which posterior quantiles should be stored:
q <- seq(0.01,0.99,by=0.01) #c(0,.025,.25,.5,.75,.975,1)
# Which posterior quantiles of Smis should be validated?
Smis.q.check <- q
max.lines           <- 1000
max.mcmc.plots      <- ifelse(sim_num<1005,"all",10)      # "all" or value
###################

################################################################

########################################
###### Compute pi(theta,Smin) ##########
########################################

if(do.profile){
  Rprof(profile_file)
}

cat("Determining whether to pre-compute pi(theta)?...\n")

pi <- NULL
if (precompute.pi.theta){
  #pi.file.location    <- paste(output.dir.outer,"/pi_theta_",job_num,".RData",sep="")
  
  cat("pi.theta file location specified as:\n")
  print(pi.file.location)
  cat("Does the file already exist?\n")
  if (file.exists(pi.file.location)){
    # Case 1/4: precompute==TRUE, previous file exists, Smin==fixed/random.
    cat("Yes. Loading from file...\n")
    precompute.pi.theta <- FALSE
    # Note: the following default objects will be re-loaded with pi.theta: 
    #       pi,E,fixed.Smin,fixed.bp,gamma,g.func
    load(pi.file.location) 
    
  } else {
    cat("No. Need to pre-compute pi(theta).\n")
    cat("Note: pi(theta) will be saved, and can be used for future analyses...\n")
  }
  
  if (precompute.pi.theta){
    cat("Pre-computing pi(theta)...\n")
    
#     # Split computation into segments of Smin
#     Smin.prior.mean <-  1.0*10^(-13)    # Mean = am * bm
#     Smin.prior.var  <- (8.0*10^(-14))^2  # Var  = am * bm^2
#     am2  <- Smin.prior.mean^2 / Smin.prior.var   # shape
#     bm2  <- am2 / Smin.prior.mean                 # rate
    
    if (!fixed.Smin){ 
      #### Varying Smin case: Set-up a grid over Smin
      # Determine grid based on prior probabilities, so that the grid will automatically adapt for different scenarios:
#       Smin.grid <- grid.select(a=am2, b=bm2, grid.eps=theta.grid.eps, grid.length=length.Smin,
#                                highest=qgamma(p=1.0-theta.grid.eps/2,shape=am2,rate=bm2))
      Smin.grid <- grid.select(a=am, b=bm, grid.eps=theta.grid.eps, grid.length=length.Smin,
                               highest=qgamma(p=1.0-theta.grid.eps/10,shape=am,rate=bm))

      ####### Additional step dividing Smin.grid into small intervals to be evaluated on server #######
      Smin.grid <- Smin.grid[1:M+M*(sim_num-sim_start-1)] # retrieve only a part of the grid - to be done in different sessions
      
    } else {
      Smin.grid=fixed.Smin
    } 
    
    pi.time <- system.time({
      pi <- pi.theta.compute(theta.grid=NULL, Smin.grid=Smin.grid, bp.grid=NULL, gamma=gamma, g=g.func, 
                           E=E, a=a,b=b,C=C,mu=mu,am=am,bm=bm, 
                           length.theta=length.theta, length.Smin=length.Smin, length.bp=length.bp,
                           fixed.Smin=fixed.Smin, fixed.bp=fixed.bp,
                           theta.lo=c(0.7,1.7), theta.hi=c(1.89,3.4), 
                           #Smin.hi=4*10^-13,
                           grid.eps=10^-4,
                           mc.integral=TRUE, nsamples=nsamples,
                           verbose=TRUE,
                           expand.grid=TRUE, just.return.grid=FALSE)
    })
    
    if (remove.main.grid==TRUE) {
      pi <- pi[-2]
    }
    
    cat("done.\n")
    
    cat("Elapsed time for pi(theta,Smin,bp):\n"); print(pi.time)
    save(pi,E,fixed.Smin,fixed.bp,gamma,g.func, file=pi.file.location)
  }
} # END precompute.pi.theta


if(do.pi.plot){
  xx <- pi$Smin
  xy <- pi$theta.1
  xz <- matrix(pi$pi,nrow=length(xx),ncol=length(xy))
  jpeg(output.pi.plot, quality=100,width=960,height=960)
  if(length(Smin.grid)>1){
    persp(x=xx,y=xy,z=xz,ticktype = "detailed",theta = -60, phi = 15, col = "lightblue", 
          ltheta = -120, shade = 0.75, xlab="Smin",ylab="theta",zlab="pi")
  } else {
    plot(pi$pi~pi$grid[,grep("theta",colnames(pi$grid))],type="l",main=ppaste("pi(theta), Smin=",fixed.Smin))
  }
  dev.off()
}

if (sim_num*M==length.Smin){
  just.pi.theta <- FALSE
}

if (just.pi.theta) {
  cat("Just pre-computing pi(theta) is requested. Stopping here for now.")
} else {

  #Sys.sleep(60)         #Suspend Execution for a Time Interval (seconds)
  
  cat("Collecting pieces of pi(theta,Smin) and combining them into 1 file.\n")
  # Combine pi.theta.Smin from many pieces
  t.Smin  <- numeric(0)
  t.bp    <- numeric(0)
  t.grid  <- numeric(0)
  t.pi    <- numeric(0)
  
  # for (j in 1:(100)) {
  #   partial.pi.file <- ppaste(output.dir.outer,"/pi_theta_",job_num,"_",j,".RData")
  #   load(partial.pi.file)
  #   pi$grid <- data.matrix(pi$grid)
  #   save(pi,length.S,pble,fixed.Smin,gamma,g.type,g.func, file=partial.pi.file)
  # }
  
  for (j in 1:(length.Smin/M)) {
    load(ppaste(output.dir.outer,"/pi_theta_",job_num,"_",j,".RData"))
    t.Smin  <- c(t.Smin, pi$Smin)
    t.bp    <- c(t.bp, pi$bp.1)
    if (remove.main.grid==FALSE) {
      t.grid  <- rbind(t.grid, pi$grid)
    }
    t.pi    <- c(t.pi, pi$pi)
    cat(paste("Finished adding file number: ",j,"\n",sep=""))
  }
  
  # Update pi
  pi$pi   <- t.pi
  pi$Smin <- t.Smin
  pi$bp.1 <- t.bp
  if (remove.main.grid==FALSE) {
    pi$grid <- t.grid
  }
  
  # Save new pi file
  if (remove.main.grid==FALSE) {
    save(pi,E,fixed.Smin,fixed.bp,gamma,g.func, file=ppaste(output.dir.outer,"/pi_theta_",job_num,"_nointerp.RData"))
  }
  save(pi,E,fixed.Smin,fixed.bp,gamma,g.func, file=ppaste(output.dir.outer,"/pi_theta_",job_num,"_simple.RData"))

  # Interpolate pi within grid
  if (length(a)==1 && interpolate==TRUE) {
    pi_new <- interpolate.pi.theta(pi=pi,a=a,b=b,theta.grid.eps=theta.grid.eps,length.theta=length.theta,
                                   am=am,bm=bm,fixed.Smin=fixed.Smin,length.Smin=length.Smin, 
                                   use.bp=(model=="bp"),use.mix=(model=="mix"),verbose=verbose)
    pi_old <- pi
    pi <- pi_new
    
    Smin.grid <- unique(pi$Smin)
    if (m==1) {   
      theta.grid <- unique(pi$theta)
    } else {
      theta.grid <- apply(theta.grid,2,unique)
    }
    # Save new pi file
    save(pi,E,fixed.Smin,fixed.bp,gamma,g.func, file=ppaste(output.dir.outer,"/pi_theta_",job_num,"_final.RData"))
    
  }
  
  #do.pi.plot <- TRUE
  if(do.pi.plot){
    xx <- pi$Smin
    xy <- pi$theta.1
    xz <- matrix(pi$pi,length(xx),length(xy))
    jpeg(output.pi.plot, quality=100,width=960,height=960)
    if(length(Smin.grid)>1){
      persp(x=xx,y=xy,z=xz,ticktype = "detailed",theta = -60, phi = 15, col = "lightblue", 
            ltheta = -120, shade = 0.75, xlab="Smin",ylab="theta",zlab="pi")
    } else if (length(grep("theta",colnames(pi$grid)))==1) {
      plot(pi$pi~pi$grid[,grep("theta.1",colnames(pi$grid))],type="l",main=ppaste("pi(theta.1), Smin=",Smin.grid))
    }
    dev.off()
  }
} #END adjusting pi

cat("\nAll done :)\n\n")

if (sinkit)
  sink()

if(do.profile){
  Rprof(NULL)
  print(summaryRprof(profile_file))
}

# 
# ####
# theta.grid1 <- sort(c(unique(pi$theta[,1]),seq(0.03,2.8,length.out=300)))
# theta.grid2 <- sort(c(unique(pi$theta[,2]),seq(0.03,2.8,length.out=300)))
# new.grid <- make.surface.grid(list(theta.grid2, theta.grid1))
# d <- data.frame(theta1=pi$theta[,1],theta2=pi$theta[,2],pi=pi$pi)
# dd <- acast(d, theta2~theta1, value.var="pi")
# obj <- list(y=unique(d$theta1),x=unique(d$theta2),z=dd)
# look <- interp.surface(obj,new.grid)
# pi_new <- list(pi=look, theta=cbind(theta1=new.grid[,2], theta2=new.grid[,1]),Smin=pi$Smin)
# pi_old <- pi
# pi <- pi_new
# 
# # Save new pi file
# save(pi,length.S,pble,fixed.Smin,Smin.grid,gamma,g.type,g.func, file=ppaste(output.dir.outer,"/pi_theta_",job_num,"_final.RData"))
