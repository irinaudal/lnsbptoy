# Testing access of pi object using pi.theta.get() and other routines

# load package

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

sim_num <- 1001
job_num <- 1
set.seed(762*sim_num + 1330931 + 10231*92) 

outer.dir.extra     <- paste("/bp_lns_toy_example_1bp_gpnorm_prec", sep="")
output.dir.extra    <- paste("/dataset_",sim_num,sep="")
output.dir.outer    <- paste(main.dir,outer.dir.extra,sep="")
output.dir          <- paste(output.dir.outer,output.dir.extra,sep="")
output.dir.wo.main  <- paste(outer.dir.extra,output.dir.extra,sep="")

pi.file.location    <- paste(output.dir.outer,"/pi_theta_",job_num,"_combined.RData",sep="")

load(pi.file.location)

### Define grid and coordinate to search within this grid
grid <- pi$grid
new.coord <- c(2.347e-16,  # Smin
               4.876e-13,  # bp.1
               2.098746435,# theta.2
               1.124783642)# theta.1
# E <- 400000 #190000
# gamma <- 1.6*(10^(-9))


### Faster routine than findInterval: avoid is.sorted check

# findInterval written in R as a binary search
findInterval2 <- function(x,v,all.inside=TRUE) {
  n <- length(v)
  if (x<v[1]) {
    if (all.inside){
      return(1)
    } else {
      return(0)
    }
  }
  if (x>=v[n]) {
    return (n)
  }
  i <- 1
  k <- n
  while({j = (k-i) %/% 2 + i; !(v[j] <= x && x < v[j+1])}) {
    if (x < v[j]){
      k <- j
    } else {
      i <- j+1
    }
  }
  return(j)
}

library(compiler)
#findInterval2 compilated with cmpfun
findInterval2 <- cmpfun(findInterval2)


library(microbenchmark)

microbenchmark(pi.theta.get(pi=pi, theta=new.coord[c(3,4)], Smin=new.coord[1], bp=new.coord[2],vertex.method=1),
               pi.theta.get(pi=pi, theta=new.coord[c(3,4)], Smin=new.coord[1], bp=new.coord[2],vertex.method=0),
               pi.theta.get(pi=pi, theta=new.coord[c(3,4)], Smin=new.coord[1], bp=new.coord[2],vertex.method=2),
               pi.theta.get(pi=pi, theta=new.coord[c(3,4)], Smin=new.coord[1], bp=new.coord[2],vertex.method=3), 
               times=50)


#### Code below is a complete mess

library(data.table)

# Data.table examples and testing ;)
d = data.table( pi$grid, pi=pi$pi , key=c("Smin","bp.1","theta.2","theta.1"))
d[1] # extract the first row of the data structure
head(d[,pi]) # extract start of pi column vector


# rounding=c(21,21,4,4)
# for(j in 1:length(rounding)){
#   pi$grid[,j] <- round(pi$grid[,j],rounding[j])
# }
# d = data.table( pi$grid, pi=pi$pi , key=c("Smin","bp.1","theta.2","theta.1"))
# 
# new.coord.round <- round(new.coord,rounding)
# d[new.coord.round] # CANNOT PROVIDE AN EXACT MATCH!!!!!!! does not work.


tt.get <- function(pi, new.coord, verbose=FALSE){
  mtc <- NULL
  # Try to search the shortened grid first, and then access by key
  len.bp <- pi$length$bp
  len.th <- pi$length$theta
  reps.all <- c(len.bp, len.th, len.th)
  nvars <- length(new.coord)
  new.idx <- numeric(nvars)
  new.idx[1] <- findInterval2(new.coord[1], pi$Smin, all.inside=TRUE)
  vec.bp <- pi$bp.1[ (new.idx[1]-1)*len.bp + 1:len.bp ] # sorted portion of bp vector
  new.idx[2] <- findInterval2(new.coord[2], vec.bp, all.inside=TRUE)
  new.idx[3] <- findInterval2(new.coord[3], pi$theta.2, all.inside=TRUE)
  new.idx[4] <- findInterval2(new.coord[4], pi$theta.1, all.inside=TRUE)

  ### SUPER FAST INDEX SEARCH:
  pi.idx <- sum( (new.idx-1)*c(1, cumprod(reps.all))[nvars:1] ) + 1
  
  if (verbose) {
    # matched coordinate
    mtc <- c(pi$Smin[new.idx[1]],
             vec.bp[new.idx[2]], 
             pi$theta.2[new.idx[3]], 
             pi$theta.1[new.idx[4]])
    
    #check reality
    idx1 <- pi$grid[,1]==mtc[1]
    idx2 <- pi$grid[,2]==mtc[2]
    idx3 <- pi$grid[,3]==mtc[3]
    idx4 <- pi$grid[,4]==mtc[4]
    pi$grid[idx1 & idx2 & idx3 & idx4]
    pi.idx.true <- which(idx1 & idx2 & idx3 & idx4)
    
    cat("Original coordinate is:\n"); print(new.coord)
    cat("Matched coordinate found  :\n"); print(mtc)
    cat("True matched coordinate is:\n"); print(pi$grid[idx1 & idx2 & idx3 & idx4])
    
    cat("pi.index selected:\n"); print(pi.idx)
    cat("True pi.index is :\n"); print(pi.idx.true)
  }
  
  return(list(pi.idx=pi.idx, mtc=mtc, pi=pi$pi[pi.idx]))
}

result <- tt.get(pi, new.coord, verbose=TRUE)
pi.idx <- result$pi.idx
mtc <- result$mtc

d[mtc] # nonsence!!

d[J(mtc[1], mtc[2], mtc[3], mtc[4] )] # THIS WORKED before, but not now... !!!!!!


pi$pi[pi.idx]

## TOTALLY DIFFERENT INDEX TO THE SOLUTION!!
# Data.table reorders rows somehow.....
d[pi.idx]
#Smin         bp.1  theta.2   theta.1     pi
#1: 4.256136e-17 1.975009e-13 2.606851 0.7833099 0.6952

# test function using benchmark
new.coord <- c(3.23434e-16, 4.854e-13, 3.119837, 1.799999)

microbenchmark(tt.get(pi=pi1, new.coord), 
               pi.theta.get(pi=pi,theta=new.coord[c(3,4)], Smin=new.coord[1], bp=new.coord[2],vertex.method=1),
               times=50)




