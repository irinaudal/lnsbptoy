"grid.select" <- function(a,b,grid.eps,grid.length,lowest=NULL, highest=NULL, 
                          distr.type="gamma", grid.type="seq", outer.length=ceiling(grid.length*0.1)) 
{  
  # Univariate function to create a grid based on some prior distribution
  
  if (grid.eps>=0.5) {
    stop("Please select a smaller grid.epsilon for grid.select().")
  }
  above <- below <- NULL
  
  if (grid.type=="seq") {
    # Determine grid based on prior probabilities, so that the grid will automatically adapt for different scenarios:
    if (distr.type=="gamma") {
      grid.min <- qgamma(p=grid.eps,shape=a,rate=b)
      grid.max <- qgamma(p=1.0-grid.eps,shape=a,rate=b)
    } else if (distr.type=="normal") {
      grid.min <- qnorm(p=grid.eps,mean=a,sd=b)
      grid.max <- qnorm(p=1.0-grid.eps,mean=a,sd=b)
    } else if (!is.null(lowest) && !is.null(highest)) {
      if (highest<lowest) {
        stop("There was an issue generating boundary values for the grid.select(): highest boundary was smaller than lowest.")
      }
      grid.min <- lowest
      grid.max <- highest
    }
    middle <- seq(grid.min, grid.max, length.out=grid.length-2*outer.length)
    
  } else if (grid.type=="quantile") {
    if (distr.type=="gamma") {
      middle <- qgamma(p=seq(grid.eps, 1.0-grid.eps, length.out=grid.length-2*outer.length), shape=a,rate=b)
    } else if (distr.type=="normal") {
      middle <- qnorm(p=seq(grid.eps, 1.0-grid.eps, length.out=grid.length-2*outer.length), mean=a,sd=b)
    }
  }
  
  # Define additional boundaries for coarse grid for pi stability
  if (outer.length>0) {
    # Add irregular grid for theta above and below:  #10% of total length.theta, 
    if (is.null(lowest)) {
      lowest <- grid.min/100
    } else if (lowest>grid.min) {
      stop("Please specify larger grid.eps OR lowest grid value to be smaller than grid.min in grid.select().")
    }
    if (is.null(highest)) {
      highest <- grid.max*10
    } else if (highest<grid.max) {
      stop("Please specify larger grid.eps OR highest grid value to be higher than grid.max in grid.select().")
    }
    if (highest<lowest) {
      stop("There was an issue generating boundary values for the grid.select(): highest boundary was smaller than lowest.")
    }
    
    below <- seq(from=lowest, to=grid.min, length.out = outer.length+1)
    above <- seq(from=grid.max, to=highest,length.out = outer.length+1)
  }
  
  grid <- unique(c(below,middle,above))
    
  return(grid)
}



# 
# "data.prep" <- function(pi,var1,var2=NULL,var3=NULL,var4=NULL) 
# {
#   if (is.null(var2)) {
#     nvar <- 1
#   } else if (is.null(var3) & !is.null(var2)) {
#     nvar <- 2
#   } else if (is.null(var4) & !is.null(var3) & !is.null(var2)) {
#     nvar <- 3
#   } else if (!is.null(var2) & !is.null(var3) & !is.null(var4)) {
#     nvar <- 4
#   } else {
#     stop("Please use consecutive assignment to 'vars'.")
#   }
#   # x is a list of all variables that pi is a function of.
#   if (nvar==1) {
#     x <- list("x"=var1)
#   } else if (nvar==2) {
#     x <- list("x"=var1,
#               "y"=var2)
#   } else if (nvar==3) {
#     x <- list("x"=var1,
#               "y"=var2,
#               "u"=var3)
#   } else if (nvar==4) {
#     x <- list("x"=var1,
#               "y"=var2,
#               "u"=var3,
#               "v"=var4)
#   }
#   # Z is an array form of pi, evaluates for each combination of x entries 
#   Z <- array(NA,dim=sapply(x,length))
#   ii <- 0
#   if (nvar==1) {
#     for (i in 1:length(x$x)){
#       ii <- ii+1
#       Z[i] <- pi[ii]
#     }
#   } else if (nvar==2) {
#     for (i in 1:length(x$x)){
#       for (j in 1:length(x$y)){
#         ii <- ii+1
#         Z[i,j] <- pi[ii]
#       }
#     }
#   } else if (nvar==3) {
#     for (i in 1:length(x$x)){
#       for (j in 1:length(x$y)){
#         for (k in 1:length(x$u)){
#           ii <- ii+1
#           Z[j,i,k] <- pi[ii]
#         }
#       }
#     }
#   } else if (nvar==4){
#     for (i in 1:length(x$x)){
#       for (j in 1:length(x$y)){
#         for (k in 1:length(x$u)){
#           for (l in 1:length(x$v)){
#             ii <- ii+1
#             Z[i,j,k,l] <- pi[ii]
#           }
#         }
#       }
#     } 
#   }
#   
#   ret <- list("x"=x, "Z"=Z)
#   return(ret)
# }
