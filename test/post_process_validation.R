##
## Post-process a large-scale validation run:
##

# Start from clean workspace:
rm(list=ls())

print(Sys.info()["sysname"])


# Setup correct paths, R history etc...
if (Sys.info()["sysname"]=="Darwin"){
  
  # Note: main.dir should NOT have a trailing /
  main.dir   <- "~/Dropbox/Log_N_Log_S/R_Code"
  setwd(main.dir)
  libsrc.dir <- "."
#  library(logNlogS)
  library(lnsbptoy)
  
  # } else if (Sys.info()["sysname"]=="Linux"){
  #   
  #   libdir <- "/nfs/home/P/pbaines/.R/library-x86_64/"
  #   library(logNlogS,lib.loc=libdir)
  #   # NoteL main.dir should NOT have a trailing /
  #   main.dir   <- "/nfs/home/P/pbaines/Rlib/lns_ms_simulation"
  #   setwd(main.dir)
  #   libsrc.dir <- "/nfs/home/P/pbaines/Rlib"
  #   
} else if (Sys.info()["sysname"]=="Linux"){
  
  libdir <- "/home/isudal/R/x86_64-pc-linux-gnu-library/3.0/"
#  library(logNlogS,lib.loc=libdir)
  library(lnsbptoy,lib.loc=libdir)
  main.dir   <- "/home/isudal/bplns"
  setwd(main.dir)
  libsrc.dir <- "."
  
}


## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)
args <- 1

cat(ppaste("Command-line arguments:\n"))
print(args)


# Read in all of the quantiles:
total_num_sims     <- 200#4#1100
num_sims_completed <- 200#4#1100
sim_start <- 1000 #1200

###################################
## Application Specific Quantities: 
##
## Specify filepaths:

if (length(args)==0){
  fp.input <- file.path(main.dir,"val_table_regular")
  fp.dataend <- "regular_misspecas_regular_jobnum_13006/store_image.RData"
  fp.output <- file.path(main.dir,"validation/cdfn_sim_regular")
  
} else {
  if(args==1){
    fp.input <- file.path(main.dir,"bp_lns_toy_example_1bp_g_8_prec")
    fp.dataend <- "store_image.RData"
    fp.output <- file.path(main.dir,"validation/bp_lns_toy_example_1bp_g_8_prec_ex1")
    total_num_sims     <- 500#4#1100
    num_sims_completed <- 500#4#1100
  } else if(args==2){
    fp.input <- file.path(main.dir,"bp_lns_toy_example_pareto_g08_sigma005_fixedNtheta_const")
    fp.dataend <- "store_image.RData"
    fp.output <- file.path(main.dir,"validation/bp_lns_toy_example_pareto_g08_sigma005_fixedNtheta_const")
    total_num_sims     <- 500#4#1100
    num_sims_completed <- 500#4#1100
  } else if(args==3){
    fp.input <- file.path(main.dir,"bp_lns_toy_example_pareto_gpnorm_mc10k")
    fp.dataend <- "store_image.RData"
    fp.output <- file.path(main.dir,"validation/bp_lns_toy_example_pareto_mc10k")
    total_num_sims     <- 500#4#1100
    num_sims_completed <- 500#4#1100
  } 
}



cover.plot <- file.path(fp.output,"coverage_line_plot.pdf")
smis.cover.plot <- file.path(fp.output,"Smis_coverage_line_plot.pdf")
interval.plot.theta.90 <- file.path(fp.output,"posterior_interval_plot_theta_90.pdf")
interval.plot.theta.96 <- file.path(fp.output,"posterior_interval_plot_theta_96.pdf")
interval.plot.N.90 <- file.path(fp.output,"posterior_interval_plot_N_90.pdf")
interval.plot.N.96 <- file.path(fp.output,"posterior_interval_plot_N_96.pdf")
interval.plot.Smin.90 <- file.path(fp.output,"posterior_interval_plot_Smin_90.pdf")
interval.plot.Smin.96 <- file.path(fp.output,"posterior_interval_plot_Smin_96.pdf")
interval.plot.bp.90 <- file.path(fp.output,"posterior_interval_plot_bp_90.pdf")
interval.plot.bp.96 <- file.path(fp.output,"posterior_interval_plot_bp_96.pdf")
verbose2 <- TRUE
###################################

## Read in the posterior quantiles and true parameter values:

qlist      <- list(NULL)
qlist.Smis <- list(NULL)
truth      <- list(NULL)
truth.Smis <- list(NULL)

num_so_far <- 0

for (iter in 1:num_sims_completed)
{
  sim_num <- sim_start + iter
  qs <- try(load(file.path(fp.input,ppaste("dataset_",sim_num),fp.dataend)),silent=TRUE)
  if (class(qs)!="try-error")
  {
    num_so_far <- num_so_far+1
    qlist[[num_so_far]]      <- final_product$post.q
    qlist.Smis[[num_so_far]] <- final_product$post.q.S.mis
    truth[[num_so_far]]      <- as.numeric(final_product$true.par)
    names(truth[[num_so_far]]) <- names(final_product$true.par)
    truth.Smis[[num_so_far]] <- as.numeric(final_product$true.S.mis)
    
    if (length(truth.Smis[[num_so_far]])<2) {
      cat(ppaste("This is dataset # 1000 + ",num_so_far,". The vector of truth.Smis = \n"))
      print(truth.Smis[[num_so_far]])
    }
    
    # Remove missing source entries if either no missing sources present
    if (is.null(final_product$post.q.S.mis)) {
      truth.Smis[[num_so_far]] <- NULL
    }
    if (is.null(final_product$true.S.mis)) {
      qlist.Smis[[num_so_far]] <- NULL
    }
    
    # fixed.pars
    if(input.pars$fixed.bp!=FALSE && any(rownames(final_product$post.q)!=names(final_product$true.par))){
      nameid.bp.true <- grep("tau.2", names(final_product$true.par))
      truth[[num_so_far]] <- truth[[num_so_far]][-nameid.bp.true]
    }
    if(input.pars$fixed.Smin!=FALSE && any(rownames(final_product$post.q)!=names(final_product$true.par))){
      nameid.sm.true <- grep("tau.1", names(final_product$true.par))
      truth[[num_so_far]] <- truth[[num_so_far]][-nameid.sm.true]
    }
    if(any(input.pars$fixed.theta!=FALSE)){
      nameid.th.true <- grep("theta", names(truth[[num_so_far]]))
      truth[[num_so_far]] <- truth[[num_so_far]][-nameid.th.true]
      qlist[[num_so_far]] <- qlist[[num_so_far]][-nameid.th.true,]
    }
    if(input.pars$fixed.N!=FALSE){
      nameid.N.true <- grep("N", names(truth[[num_so_far]]))
      truth[[num_so_far]] <- truth[[num_so_far]][-nameid.N.true]
      qlist[[num_so_far]] <- qlist[[num_so_far]][-nameid.N.true,]
    }
    
    # Log-Posterior
    nameid.logPost.true <- grep("log.Poster", names(truth[[num_so_far]]))
    nameid.logPost.post <- grep("log.Poster", rownames(qlist[[num_so_far]]))
    nameid.Sobs <- grep("S.obs", names(truth[[num_so_far]]))
    
    if (length(nameid.logPost.true)>0 && length(nameid.logPost.post)>0){
      # Shift log-posterior values before Sobs
      idx <- c(1:(min(nameid.Sobs)-1),nameid.logPost.true,nameid.Sobs)
      qlist[[num_so_far]] <- qlist[[num_so_far]][idx,]
      truth[[num_so_far]] <- truth[[num_so_far]][idx]
      
    } else if (length(nameid.logPost.post)>0) {
      # Remove log-posterior from qlist
      qlist[[num_so_far]] <- qlist[[num_so_far]][-nameid.logPost.post,]
    }      
    
    if(verbose2){
      cat("head(truth vector): \n")
      print(round(head(truth[[num_so_far]]),3))
    }
    
  } else {
    num_so_far <- num_so_far+1
    qlist[[num_so_far]]      <- NULL
    qlist.Smis[[num_so_far]] <- NULL
    truth[[num_so_far]]      <- NULL
    truth.Smis[[num_so_far]] <- NULL
  }
  
  if ( ((iter%%100)==0) & verbose2) {
    cat(ppaste("Dataset ", iter, " has been read.\n"))
  }
  
} # END for loop over iter: num_sims_completed


cat("Computing Validation statistics of MCMC draws...\n")   

# Print output
# cat("\nSimulation specs:\n")
# cat(paste("nvalid = ",total_num_sims,"\n",sep="")) 
# cat(paste("niter  = ",niter,"\n",sep="")) 
# cat(paste("burnin = ",burnin,"\n",sep="")) 
# cat(paste("z      = ",z,"\n",sep="")) 
# cat(paste("gamma  = ",gamma,"\n",sep="")) 
# cat(paste("Smin   = ",Smin," (Mimimum Flux threshold)\n",sep="")) 
# cat(paste("Emin   = ",Emin,"\n",sep="")) 
# cat(paste("Emax   = ",Emax,"\n",sep="")) 
# cat(paste("alpha  = ",alpha," (NegBinom prior of N)\n",sep="")) 
# cat(paste("beta   = ",beta," (NegBinom prior of N)\n",sep="")) 
# cat(paste("a      = ",a," (Gamma prior of theta)\n",sep="")) 
# cat(paste("b      = ",b," (Gamma prior of theta)\n",sep="")) 
# cat("Marginal probability of observing source, pi(theta):\n")
# print(t(pi.theta))
# print(cbind(n.totals,true.theta))

## Compute coverage proportions:
## "cover_func" <- function(truth,X,Y,type="one-sided",return.type=c('all','proportions','indicators'))

# Find first successful non-NULL qlist:
for (i in 1:(num_sims_completed)){
  if (!is.null(qlist[[i]])){
    all.par.names <- rownames(qlist[[i]])
    break
  }
}

N.idx <- grep("N",all.par.names)
sobs.index.start <- min(grep("S.obs",all.par.names))
fixed.par.names <- all.par.names[1:(sobs.index.start-1)]

coverage_out      <- cover_func(X=qlist,truth=truth,type='one-sided',return.type='all',
                                fixed.dims=c(1:(sobs.index.start-1)), nonfixed.names="S.obs",
                                discrete.dims=N.idx, verbose=verbose2)
coverage_out.Smis <- cover_func(X=qlist.Smis,truth=truth.Smis,type='one-sided',return.type='all',fixed.dims=NULL,verbose=verbose2)

cover_indicators  <- coverage_out$indicators
prop_cover        <- coverage_out$proportions
successful_dataset_indices <- coverage_out$successful_dataset_indices
fixed.par.names   <- rownames(coverage_out$proportions)

cover_indicators.Smis  <- coverage_out.Smis$indicators
prop_cover.Smis        <- coverage_out.Smis$proportions
successful_dataset_indices.Smis <- coverage_out.Smis$successful_dataset_indices

num_sims_successfully_completed <- length(successful_dataset_indices)
num_sims_successfully_completed.Smis <- length(successful_dataset_indices.Smis)

####
## logical(0) for datasets not analyzed
####

cat("\nThere were ",num_sims_successfully_completed,
    " datasets that were successfully analyzed.\n\n",sep="")
cat("\nThere were ",num_sims_successfully_completed.Smis,
    " datasets that had MCMC samples of S.mis.\n\n",sep="")

## Create nice looking table for results:

#fbit <- substr(Sys.time(),start=12,stop=19)
#write.table(file=ppaste("cover_",mu_0,fbit,".txt"),prop_cover,row.names=TRUE,col.names=TRUE,quote=FALSE)

## Sometimes need to average across columns:
p_table <- prop_cover
round(p_table,3)

## Print latex output in convenient form:
#library(xtable)
#xtable(round(p_table,3))


## Retrieve failed indicators
nameid <- 1:(nameid.Sobs[1]-1) # for nameid.N
#names  <- c("N","theta.1","theta.2","tau.1","tau.2")
names  <- names(truth[[1]])[nameid]

conf.level=0.96
ci.lo.int <- (1-conf.level)/2*100
ci.hi.int <- (conf.level)*100 + ci.lo.int
#ix.bad <- matrix(NA,nrow=length(successful_dataset_indices), ncol=num_sims_completed)
true.val <- numeric(num_sims_completed)
q.hi <- true.val
q.lo <- true.val
for(j in 1:length(nameid)){
  cat("Creating vector of true values of parameter...\n")
  true.val <- unlist(lapply(truth,function(x){if (!is.null(x) && length(x)>=2){x[nameid[j]]} else {NA}}))  
  cat("Computing central (1-conf.level)*100% posterior intervals...\n")
  q.lo <- unlist(lapply(qlist,function(x){if (!is.null(x) && ncol(x)>=ci.lo.int){x[names[j],ppaste(ci.lo.int,"%")]} else {NA}}))
  q.hi <- unlist(lapply(qlist,function(x){if (!is.null(x) && ncol(x)>=ci.lo.int){x[names[j],ppaste(ci.hi.int,"%")]} else {NA}}))
  mark.below <- which(q.hi<true.val) 
  mark.above <- which(q.lo>true.val)
  cat(paste("Indices of",num_sims_completed,"credible intervals at",conf.level*100, "% confidence for", names[j], "are below true value:\n"))
  print( mark.below )
  cat(paste("Indices of",num_sims_completed,"credible intervals at",conf.level*100, "% confidence for", names[j], "are above true value:\n"))
  print( mark.above )
}  

## Make some plots of the results:

# actual_coverage must be a matrix, each column corresponding to a different parameter,
# each row corresponding to the coverage rate of the quantile specified in desired_coverage.
#actual_coverage <- matrix(nrow=3,ncol=ncol(prop_cover))
#actual_coverage[1:2,]  <- prop_cover[1:2,]
#actual_coverage[3,] <- apply(prop_cover[3:nrow(prop_cover),,drop=FALSE],2,mean)

actual_coverage      <- prop_cover
actual_coverage.Smis <- prop_cover.Smis
  
rownames(actual_coverage.Smis) <- NULL
desired_coverage <- #c(25,10,5,2.5,1,0.5,0.05, 75,90,95,97.5,99,99.5,99.95)/100
  seq(0.01,0.99,by=0.01)
#c(0.5,1,2.5,5,10,25, 40,50,60, 75,90,95,97.5,99,99.5)/100

## SORT COVERAGE ORDER INSIDE global_coverage_plot

cat("Making global coverage plot for theta, N and Sobs...\n")

pdf(cover.plot)
global_coverage_plot(desired_coverage=desired_coverage,actual_coverage=t(actual_coverage),bands=TRUE,num_sims=num_sims_successfully_completed)
dev.off()

cat("Making global coverage plot for Smis...\n")

pdf(smis.cover.plot)
global_coverage_plot(desired_coverage=desired_coverage,actual_coverage=t(actual_coverage.Smis),bands=TRUE,num_sims=num_sims_successfully_completed.Smis)
dev.off()

# Making posterior for theta(s)
nameid.theta <- grep("theta",names(truth[[1]]))
m <- length(nameid.theta)
names.theta <- paste("theta.",1:m,sep="")

# Making posterior for N
nameid.N <- 1
names.N <- "N"
# nameid.Smin <- grep("Smin",names(final_product$true.par))
# names.Smin <- "Smin"
# nameid.bp <- grep("bp",names(final_product$true.par))
# names.bp <- "bp"
nameid.Smin <- grep("tau.1",names(truth[[1]]))
names.Smin <- "tau.1"
nameid.bp <- grep("tau.2",names(truth[[1]]))
names.bp <- "tau.2"

###############
pdf(interval.plot.theta.90)
generic_posterior_interval_plot(qlist=qlist,truth=truth,nameid=nameid.theta,names=names.theta,
                                successful_dataset_indices=successful_dataset_indices,conf.level=0.90)
dev.off()

pdf(interval.plot.theta.96)
ix.bad <- generic_posterior_interval_plot(qlist=qlist,truth=truth,nameid=nameid.theta,names=names.theta,
                                          successful_dataset_indices=successful_dataset_indices,
                                          conf.level=0.96,ret.bad.ix=TRUE)
dev.off()

# xxt <- which(ix.bad)
# yy  <- sapply(1:length(xxt),function(j){quantile(qlist[[xxt[j]]][nameid.theta,],prob=c(0.02,0.98))})
# tt  <- sapply(1:length(xxt),function(j){(truth[[xxt[j]]][nameid.theta])})
# hi.ci <- tt < yy[1,]
# rr  <- round(cbind(truth=tt,hi.ci,t(yy), diff=yy[2,]-yy[1,]),3)
# srt <- sort(tt,index=TRUE)$ix
# rr <- rr[srt,]
# print(cbind(rr,rr[,5]<0.18))

###############
pdf(interval.plot.N.90)
generic_posterior_interval_plot(qlist=qlist,truth=truth,nameid=nameid.N,names=names.N,successful_dataset_indices=successful_dataset_indices,conf.level=0.90)
dev.off()

pdf(interval.plot.N.96)
in.bad <- generic_posterior_interval_plot(qlist=qlist,truth=truth,nameid=nameid.N,names=names.N,
                                          successful_dataset_indices=successful_dataset_indices,
                                          conf.level=0.96,ret.bad.ix=TRUE)
dev.off()

# xxn <- which(in.bad)
# yy  <- sapply(1:length(xxn),function(j){quantile(qlist[[xxn[j]]][nameid.N,],prob=c(0.02,0.98))})
# nn  <- sapply(1:length(xxn),function(j){(truth[[xxn[j]]][nameid.N])})
# hi.ci <- nn < yy[1,]
# rr  <- round(cbind(truth=nn,hi.ci,t(yy), diff=yy[2,]-yy[1,]),3)
# srt <- sort(nn,index=TRUE)$ix
# rr <- rr[srt,]
# print(cbind(rr,rr[,5]<150))


####################
pdf(interval.plot.Smin.90)
generic_posterior_interval_plot(qlist=qlist,truth=truth,nameid=nameid.Smin,names=names.Smin,
                                successful_dataset_indices=successful_dataset_indices,
                                conf.level=0.90,ret.bad.ix=FALSE)
dev.off()

pdf(interval.plot.Smin.96)
im.bad <- generic_posterior_interval_plot(qlist=qlist,truth=truth,nameid=nameid.Smin,names=names.Smin,
                                          successful_dataset_indices=successful_dataset_indices,
                                          conf.level=0.96,ret.bad.ix=TRUE)
dev.off()

# xxn <- which(im.bad)
# yy  <- sapply(1:length(xxn),function(j){quantile(qlist[[xxn[j]]][nameid.Smin,],prob=c(0.02,0.98))})
# nn  <- sapply(1:length(xxn),function(j){(truth[[xxn[j]]][nameid.Smin])})
# hi.ci <- nn < yy[1,]
# rr  <- round(cbind(truth=nn,hi.ci,t(yy), diff=yy[2,]-yy[1,]),3)
# srt <- sort(nn,index=TRUE)$ix
# rr <- rr[srt,]
# print(cbind(rr,rr[,5]<100))

####################
pdf(interval.plot.bp.90)
generic_posterior_interval_plot(qlist=qlist,truth=truth,nameid=nameid.bp,names=names.bp,
                                successful_dataset_indices=successful_dataset_indices,
                                conf.level=0.90,ret.bad.ix=FALSE)
dev.off()

pdf(interval.plot.bp.96)
im.bad <- generic_posterior_interval_plot(qlist=qlist,truth=truth,nameid=nameid.bp,names=names.bp,
                                          successful_dataset_indices=successful_dataset_indices,
                                          conf.level=0.96,ret.bad.ix=TRUE)
dev.off()
