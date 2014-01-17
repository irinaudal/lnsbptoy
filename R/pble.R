########################## Define Classes ######################
# First, the class 'pble' is defined. It contains one field : par .
setClass("pble.class",
         representation(pars="list",
                        grid="list",
                        BLE.grid="data.frame",
                        "VIRTUAL"      # cannot use pble class by itself, only its children
         ))
# 'pble.basic' is a subClass of 'pble.class', it contains no more fields:
setClass("pble.basic",
         representation(prob="numeric"
         ),
         contains="pble.class",
         validity=function(object){
           if (length(object@pars)!=15) {
             stop ("[pble.basic: validation] the number of parameters in the list must be 15: each of B,L,E with min,max,length,fixed,d(elta).")
           } else {
             # Checking values of parameters
             checkValidPars(object@pars,object@grid,verbose=FALSE)
           }
           if (object@prob<0){
             stop ("Probability in object[\"prob\"] cannot be < 0.")
           } else if (object@prob>1){
             stop ("Probability in object[\"prob\"] cannot be > 1.")
           } else {}
           return(TRUE)
         })

# 'pble.table' is another subClass of 'pble.pble' it contains 3 more fields: 'table', 'table_orig', 'filename' :
setClass("pble.table",
         representation(table="data.frame",
                        table_orig="data.frame",
                        filename="character"
         ),
         contains="pble.class",
         validity=function(object){
           # Checking values of parameters
           checkValidPars(object@pars,object@grid,verbose=FALSE)
           return(TRUE)
         })

########################## Set Methods ######################
# Let's define a method for easy printing of objects 'pble.basic' and 'pble.table'.
setMethod("show", "pble.basic",
          function(object) {         
            cat("Basic object class of p(BLE)\n")
            cat("  Type :", class(object), "\n")
            idx <- grep("B",names(object@pars))
            cat("  Pars (Background) \n")
            printPars(object,idx)
            idx <- grep("L",names(object@pars))
            cat("       (Location) \n")
            printPars(object,idx)
            idx <- grep("E",names(object@pars))
            cat("       (Exposure Map) \n")
            printPars(object,idx)
            cat("  Prob (Uniform Probability Density) \n")
            cat(paste("             p(B,L,E) = ",object@prob,"\n",sep=""))
            cat("  Grid :\n")
            idx <- grep("B",names(object@grid))
            printGrid(object,idx)
            idx <- grep("L",names(object@grid))
            printGrid(object,idx)
            idx <- grep("E",names(object@grid))
            printGrid(object,idx)
            cat("\n")
            cat(paste("  BLE.grid :   (Long format file partially shown) \n         with dim(BLE.grid)=[",nrow(object@BLE.grid),"x",ncol(object@BLE.grid),"]\n",sep=""))
            print(head(object@BLE.grid))
            cat("\n")
          })

setMethod("show", "pble.table",
          function(object) {
            cat("Table object class of p(BLE)\n")
            cat("  Type :", class(object), "\n")
            idx <- grep("B",names(object@pars))
            cat("  Pars (Background) \n")
            printPars(object,idx)
            idx <- grep("L",names(object@pars))
            cat("       (Location) \n")
            printPars(object,idx)
            idx <- grep("E",names(object@pars))
            cat("       (Exposure Map) \n")
            printPars(object,idx)
            cat("  Grid :\n")
            idx <- grep("B",names(object@grid))
            printGrid(object,idx)
            idx <- grep("L",names(object@grid))
            printGrid(object,idx)
            idx <- grep("E",names(object@grid))
            printGrid(object,idx)
            cat("  Table_orig : \n")
            print(object@table_orig)
            cat("  Table    : (Long format file partially shown)\n")
            print(head(object@table))
            cat("  Filename :",object@filename,"\n")
            cat("\n")
          })

# Define the "[" Accessor as a Method for our class.
# Note: since "[" has already been defined by R, we cannot change its argument names.
setMethod(f="[", signature(x="pble.basic"),
          function(x,i,j,drop){
            if(i=="pars"){return(x@pars)}
            if(i=="prob"){return(x@prob)}
            if(i=="grid"){return(x@grid)}    
            if(i=="BLE.grid"){return(x@BLE.grid)}    
          })

setMethod(f="[",signature(x="pble.table"),
          function(x,i,j,drop){
            if(i=="pars"){return(x@pars)}
            if(i=="table"){return(x@table)}
            if(i=="table_orig"){return(x@table_orig)}
            if(i=="filename"){return(x@filename)}
            if(i=="grid"){return(x@grid)}
            if(i=="BLE.grid"){return(x@BLE.grid)}
          })

# # Define LIMITED Replace Method for replacing "pars" our class.
# setReplaceMethod(f="[",signature(x="pble.basic"),
#                  function(x,i,j,value){
#                    if(i=="pars"){x@pars<-value}else {}
#                    validObject(x)
#                    return (x)
#                  })
# setReplaceMethod(f="[",signature(x="pble.table"),
#                  function(x,i,j,value){
#                    if(i=="pars"){x@pars<-value}else {}
#                    validObject(x)
#                    return (x)
#                  })

# Let's define a generic function 'sample.pble' - allows for overloading of this function.
"sample.pble" <- function(object,n=NULL) {
  stop("Please specify sample size 'n'; or Object of Unknown data type! (Consider passing class 'pble.basic' or 'pble.table'")
}
setGeneric("sample.pble")

# Let's define a 'generic' function 'evaluate.pble' - allows for overloading of this function.
"evaluate.pble" <- function(object,B=NULL,L=NULL,E=NULL) {
  stop("Unknown data type! (Consider passing class 'pble.basic' or 'pble.table'")
}
setGeneric("evaluate.pble")

# Check functionality of this generic function
#showMethods(sample.pble)
#class(sample.pble)
#showMethods(evaluate.pble)
#class(evaluate.pble)

# Then, define S4 methods 'sample.pble' and 'evaluate.pble', dispatched when called for 'pble.basic' and 'pble.table' class.
setMethod("sample.pble",signature(object="pble.basic"),
          function(object,n) {
            E <- runif(n,object@pars$Emin,object@pars$Emax)  
            L <- runif(n,object@pars$Lmin,object@pars$Lmax)
            B <- runif(n,object@pars$Bmin,object@pars$Bmax)
            return(list("B"=B, "L"=L, "E"=E))
          })

setMethod("evaluate.pble",signature(object="pble.basic"),
          function(object) {
            # 3-D uniform distribution density at B,L,E
            return(object@pars$prob)
          })

setMethod("sample.pble",signature(object="pble.table"),
          function(object,n) {          
            # Identify cum-sum-interval of probabilities
            cumct <- cumsum(object@table$prob)
            # Sample from mltinomial distribution of pble 'prob' using Uniform sampling approach
            U     <- runif(n,0,1)
            res   <- object@table[findInterval(U,cumct)+1,]
            # cat('Inside sample.pble for table object - sampling index row of prob.table:\n'); cat('res\n:'); print(res)
            
            # Add uniform errors within each grid inverval of B,L,E
            res$B <- res$B+runif(n,0,object@pars$dB)
            res$L <- res$L+runif(n,0,object@pars$dL)
            res$E <- res$E+runif(n,0,object@pars$dE)           
            return(list("B"=res$B, "L"=res$L, "E"=res$E))
          })

setMethod("evaluate.pble",signature(object="pble.table"),
          function(object,B,L,E) {           
            # Find closest values of grid to B,L,E using binary search
            i   <- findInterval(B,object@grid$B.grid)
            j   <- findInterval(L,object@grid$L.grid)
            k   <- findInterval(E,object@grid$E.grid)
            idx <- 1+(k-1) + (j-1)*object@pars$length.E + (i-1)*object@pars$length.E*object@pars$length.L
            res <- object@table$prob[idx]            
            # # Interpolate within grid
            # res <- interpVoxel(object,B,L,E,i,j,k)
            return(res)
          })

# Check functionality of this generic function
#showMethods(sample.pble)
#showMethods(evaluate.pble)

########################## Set Utility functions ######################
# Utility function for interpolation over 3-d grid
"interpVoxel" <- function(object,B,L,E,i,j,k){
  # The following code is necessary for Interpolation of prob inside Voxel
  # Let Y = V[b,l,e] voxel boundaries
  # Find closest values of grid to B,L,E using binary search
  BLEidx <- expand.grid(E=k:(k+1), L=j:(j+1), B=i:(i+1))               
  id <- 1+(BLEidx$E-1) + (BLEidx$L-1)*object["pars"]$length.E + (BLEidx$B-1)*object["pars"]$length.E*object["pars"]$length.L
  Y <- object["table"]$prob[id]            
  
  #   # Linear interpolation between grid of table
  #   # Requires: library(MASS)
  #   x <- expand.grid(B=b,L=l,E=e)
  #   X <- as.matrix(data.frame(inter=rep(1,nrow(x)),x))
  #   Xnew <- matrix(c(1,B,L,E),nrow=1)
  #   res  <- Xnew %*% ginv(t(X) %*% X) %*% t(X) %*% Y
  
  # Trilinear Interpolation <http://en.wikipedia.org/wiki/Trilinear_interpolation>
  temp <- object@table[id[1],]
  dB <- (B - temp$B)/object@pars$dB
  dL <- (L - temp$L)/object@pars$dL
  dE <- (E - temp$E)/object@pars$dE
  # Interpolate along B-direction
  i1 <- Y[1]*(1-dB) + Y[2]*dB
  i2 <- Y[3]*(1-dB) + Y[4]*dB
  j1 <- Y[5]*(1-dB) + Y[6]*dB
  j2 <- Y[7]*(1-dB) + Y[8]*dB    
  # Interpolate along L-direction
  w1 <- i1*(1-dL) + i2*dL
  w2 <- j1*(1-dL) + j2*dL
  # Interpolate along E-direction
  res <- w1*(1-dE) + w2*dE
  return(res)
}

# Utility function for printing parameters of objects
"printPars" <- function(object,idx) {
  for(j in 1:length(idx)){
    value <- object@pars[[idx[j]]]
    value <- ifelse(is.null(value),"NULL",value)
    cat(paste("           :",format(names(object@pars)[idx[j]],width=8),"=",value,sep=" "),"\n")
  }
}     

# Utility function for printing parameters of objects
"printGrid" <- function(object,idx) {
  value <- object@grid[[idx]]
  if (is.null(value)) {
    value <- "NULL"
  } else if (length(value)>6) {
    value <- paste(paste(round(value[1:3],digits=2),collapse=" "),"...",
                   paste(round(value[(length(value)-2):length(value)],digits=2),collapse=" "))
  }
  cat(paste("           :",format(names(object@grid)[idx],width=8),"=",value,sep=" "),"\n")
}

# Utility function for ckecking parameters of variables 
"checkValidPars" <- function(pars,grid,verbose=FALSE) {
  if (verbose){
    cat("Object Inspector: inspecting 'pars' & 'grid'... ")
  }
  Bmin     <- pars$Bmin
  Bmax     <- pars$Bmax
  fixed.B  <- pars$fixed.B
  length.B <- pars$length.B
  dB       <- pars$dB
  Lmin     <- pars$Lmin
  Lmax     <- pars$Lmax
  fixed.L  <- pars$fixed.L
  length.L <- pars$length.L
  dL       <- pars$dL
  Emin     <- pars$Emin
  Emax     <- pars$Emax
  fixed.E  <- pars$fixed.E
  length.E <- pars$length.E
  dE       <- pars$dE
  B.grid   <- grid$B
  L.grid   <- grid$L
  E.grid   <- grid$E
  # Perform basic error checking
  if (is.null(dB)) {
    stop("'dB' must be a positive VALUE.")
  } else if(dB<=0) {
    stop("'dB' must be > 0.")
  }
  if (is.null(B.grid)) {
    stop("'B.grid' must be an array of VALUEs.")
  }
  if (is.null(fixed.B)) { 
    if (is.null(length.B)) {
      stop("One of 'fixed.B' or 'length.B' must be <not> NULL.")
    }
    if (Bmax < Bmin) {
      stop("'Bmax' must be larger than 'Bmin'.")
    }
    if (Bmax==Bmin){
      if (length.B!=1) {
        stop("Since 'fixed.B'=NULL, Bmin==Bmax, then 'length.B' must be equal to 1.")
      }
      if (fixed.B!=Bmin) {
        stop("Since Bmin==Bmax, 'fixed.B' should be set to 'Bmin'.")
      }
      if (dB!=1.0) {
        stop("Since Bmin==Bmax, 'dB' should be set to 1.")
      }
      if (B.grid!=Bmin) {
        stop("Since Bmin==Bmax, 'B.grid' should be set to 'Bmin'.")
      }
    } else {
      if (length.B < 2) {
        stop("Since 'fixed.B'=NULL, Bmin<Bmax, then 'length.B' must be greater than or equal to 2.")
      }
      if (dB!=(Bmax-Bmin)/(length.B-1)) {
        stop("Since Bmin==Bmax, 'dB' should be set to 1.")
      }
      if (!all(B.grid==seq(Bmin,Bmax,length.out=length.B))) {
        stop("Check definition of 'B.grid'.")
      }
    }   
  } else if (length(fixed.B)==1) {  
    if (is.null(length.B) || (length.B!=1)){
      stop("Since 'fixed.B'=VALUE, then 'length.B' must be equal to 1.")
    } 
    if ((Bmax!=fixed.B) || (Bmin!=fixed.B)) {
      stop("Since 'fixed.B'=VALUE, then 'Bmin' and 'Bmax' must be equal to 'fixed.B'.")
    }
    if (dB!=1.0) {
      stop("Since Bmin==Bmax, 'dB' should be set to 1.")
    }
    if (B.grid!=Bmin) {
      stop("Since Bmin==Bmax, 'B.grid' should be set to 'Bmin'.")
    }
  } else { stop("Please check 'fixed.B' variable: must be a VALUE of length=1 or a NULL.") }
  ######################################################  
  if (is.null(dL)) {
    stop("'dL' must be a positive VALUE.")
  } else if(dL<=0) {
    stop("'dL' must be > 0.")
  }
  if (is.null(L.grid)) {
    stop("'L.grid' must be an array of VALUEs.")
  }
  if (is.null(fixed.L)) { 
    if (is.null(length.L)) {
      stop("One of 'fixed.L' or 'length.L' must be <not> NULL.")
    }
    if (Lmax < Lmin) {
      stop("'Lmax' must be larger than 'Lmin'.")
    }
    if (Lmax==Lmin){
      if (length.L!=1) {
        stop("Since 'fixed.L'=NULL, Lmin==Lmax, then 'length.L' must be equal to 1.")
      }
      if (fixed.L!=Lmin) {
        stop("Since Lmin==Lmax, 'fixed.L' should be set to 'Lmin'.")
      }
      if (dL!=1.0) {
        stop("Since Lmin==Lmax, 'dL' should be set to 1.")
      }
      if (L.grid!=Lmin) {
        stop("Since Lmin==Lmax, 'L.grid' should be set to 'Lmin'.")
      }
    } else {
      if (length.L < 2) {
        stop("Since 'fixed.L'=NULL, Lmin<Lmax, then 'length.L' must be greater than or equal to 2.")
      }
      if (dL!=(Lmax-Lmin)/(length.L-1)) {
        stop("Since Lmin==Lmax, 'dL' should be set to 1.")
      }
      if (!all(L.grid==seq(Lmin,Lmax,length.out=length.L))) {
        stop("Check definition of 'L.grid'.")
      }
    }   
  } else if (length(fixed.L)==1) {    
    if (is.null(length.L) || (length.L!=1)){
      stop("Since 'fixed.L'=VALUE, then 'length.L' must be equal to 1.")
    } 
    if ((Lmax!=fixed.L) || (Lmin!=fixed.L)) {
      stop("Since 'fixed.L'=VALUE, then 'Lmin' and 'Lmax' must be equal to 'fixed.L'.")
    }
    if (dL!=1.0) {
      stop("Since Lmin==Lmax, 'dL' should be set to 1.")
    }
    if (L.grid!=Lmin) {
      stop("Since Lmin==Lmax, 'L.grid' should be set to 'Lmin'.")
    } 
  } else { stop("Please check 'fixed.L' variable: must be a VALUE of length=1 or a NULL.") }
  ######################################################             
  if (is.null(dE)) {
    stop("'dE' must be a positive VALUE.")
  } else if(dE<=0) {
    stop("'dE' must be > 0.")
  }
  if (is.null(E.grid)) {
    stop("'E.grid' must be an array of VALUEs.")
  }
  if (is.null(fixed.E)) { 
    if (is.null(length.E)) {
      stop("One of 'fixed.E' or 'length.E' must be <not> NULL.")
    }
    if (Emax < Emin) {
      stop("'Emax' must be larger than 'Emin'.")
    }
    if (Emax==Emin){
      if (length.E!=1) {
        stop("Since 'fixed.E'=NULL, Emin==Emax, then 'length.E' must be equal to 1.")
      }
      if (fixed.E!=Emin) {
        stop("Since Emin==Emax, 'fixed.E' should be set to 'Emin'.")
      }
      if (dE!=1.0) {
        stop("Since Emin==Emax, 'dE' should be set to 1.")
      }
      if (E.grid!=Emin) {
        stop("Since Emin==Emax, 'E.grid' should be set to 'Emin'.")
      }
    } else {
      if (length.E < 2) {
        stop("Since 'fixed.E'=NULL, Emin<Emax, then 'length.E' must be greater than or equal to 2.")
      }
      if (dE!=(Emax-Emin)/(length.E-1)) {
        stop("Since Emin==Emax, 'dE' should be set to 1.")
      }
      if (!all(E.grid==seq(Emin,Emax,length.out=length.E))) {
        stop("Check definition of 'E.grid'.")
      }
    }   
  } else if (length(fixed.E)==1) {    
    if (is.null(length.E) || (length.E!=1)){
      stop("Since 'fixed.E'=VALUE, then 'length.E' must be equal to 1.")
    } 
    if ((Emax!=fixed.E) || (Emin!=fixed.E)) {
      stop("Since 'fixed.E'=VALUE, then 'Emin' and 'Emax' must be equal to 'fixed.E'.")
    }
    if (dE!=1.0) {
      stop("Since Emin==Emax, 'dE' should be set to 1.")
    }
    if (E.grid!=Emin) {
      stop("Since Emin==Emax, 'E.grid' should be set to 'Emin'.")
    } 
  } else { stop("Please check 'fixed.E' variable: must be a VALUE of length=1 or a NULL.") }
  ######################################################   
  if (verbose){
    cat("Done.\n")
  }
  return(NULL)
}

# Utility function for ckecking+fixing parameters of variables 
"fixValidPars" <- function(Bmin, Bmax, Lmin, Lmax, Emin, Emax, length.B=2, length.L=2, length.E=20, fixed.B=NULL, fixed.L=NULL, fixed.E=NULL) {
  # Perform basic error checking
  if (is.null(fixed.B)) { 
    if (is.null(length.B)) {
      stop("One of 'fixed.B' or 'length.B' must be <not> NULL")
    }
    if (Bmax < Bmin) {
      warning("'Bmax' must be larger than 'Bmin'. Switching...")
      tmp  <- Bmax
      Bmax <- Bmin
      Bmin <- tmp
    }
    if (Bmax==Bmin){
      if (length.B!=1) {
        warning("Since 'fixed.B'=NULL, Bmin==Bmax, then 'length.B' must be equal to 1. Updating: 'length.B'=1.")
        length.B <- 1.0
      }
      # Override with 1-point grid:
      fixed.B  <- Bmin
      B.grid   <- Bmin # == Bmax
      dB       <- 1.0
    } else {
      if (length.B < 2) {
        stop("Since 'fixed.B'=NULL, Bmin<Bmax, then 'length.B' must be greater than or equal to 2")
      }
      B.grid  <- seq(Bmin,Bmax,length.out=length.B)
      dB      <- (Bmax-Bmin)/(length.B-1)
    }   
  } else if (length(fixed.B)==1) {   
    if (is.null(length.B) || (length.B!=1)){
      warning("Since 'fixed.B'=VALUE, then 'length.B' must be equal to 1. Updating: 'length.B'=1.")
      length.B <- 1.0
    } 
    if ((Bmax!=fixed.B) || (Bmin!=fixed.B)) {
      warning("Since 'fixed.B'=VALUE, then 'Bmin' and 'Bmax' must be equal to 'fixed.B'. Updating: 'Bmin','Bmax' to be 'fixed.B'.")
      Bmin <- Bmax <- fixed.B
    }
    # Override with 1-point grid:
    B.grid   <- Bmin # == Bmax
    dB       <- 1.0   
  } else { stop("Please check 'fixed.B' variable: must be a VALUE of length=1 or a NULL.") }
  ######################################################  
  if (is.null(fixed.L)) { 
    if (is.null(length.L)) {
      stop("One of 'fixed.L' or 'length.L' must be <not> NULL")
    }
    if (Lmax < Lmin) {
      warning("'Lmax' must be larger than 'Lmin'. Switching...")
      tmp  <- Lmax
      Lmax <- Lmin
      Lmin <- tmp
    }
    if (Lmax==Lmin){
      if (length.L!=1) {
        warning("Since 'fixed.L'=NULL, Lmin==Lmax, then 'length.L' must be equal to 1. Updating: 'length.L'=1.")
        length.L <- 1.0
      }      # Override with 1-point grid:
      fixed.L  <- Lmin
      L.grid   <- Lmin # == Lmax
      dL       <- 1.0
    } else {
      if (length.L < 2) {
        stop("Since 'fixed.L'=NULL, Lmin<Lmax, then 'length.L' must be greater than or equal to 2")
      }
      L.grid  <- seq(Lmin,Lmax,length.out=length.L)
      dL      <- (Lmax-Lmin)/(length.L-1)
    }   
  } else if (length(fixed.L)==1) {    
    if (is.null(length.L) || (length.L!=1)){
      warning("Since 'fixed.L'=VALUE, then 'length.L' must be equal to 1. Updating: 'length.L'=1.")
      length.L <- 1.0
    } 
    if ((Lmax!=fixed.L) || (Lmin!=fixed.L)) {
      warning("Since 'fixed.L'=VALUE, then 'Lmin' and 'Lmax' must be equal to 'fixed.L'. Updating: 'Lmin','Lmax' to be 'fixed.L'.")
      Lmin <- Lmax <- fixed.L
    }
    # Override with 1-point grid:
    L.grid   <- Lmin # == Lmax
    dL       <- 1.0    
  } else { stop("Please check 'fixed.L' variable: must be a VALUE of length=1 or a NULL.") }
  ######################################################  
  if (is.null(fixed.E)) { 
    if (is.null(length.E)) {
      stop("One of 'fixed.E' or 'length.E' must be <not> NULL")
    }
    if (Emax < Emin) {
      warning("'Emax' must be larger than 'Emin'. Switching...")
      tmp  <- Emax
      Emax <- Emin
      Emin <- tmp
    }
    if (Emax==Emin){
      if (length.E!=1) {
        warning("Since 'fixed.E'=NULL, Emin==Emax, then 'length.E' must be equal to 1. Updating: 'length.E'=1.")
        length.E <- 1.0
      }      
      # Override with 1-point grid:
      fixed.E  <- Emin
      E.grid   <- Emin # == Emax
      dE       <- 1.0
    } else {
      if (length.E < 2) {
        stop("Since 'fixed.E'=NULL, Emin<Emax, then 'length.E' must be greater than or equal to 2")
      }
      E.grid  <- seq(Emin,Emax,length.out=length.E)
      dE      <- (Emax-Emin)/(length.E-1)
    }   
  } else if (length(fixed.E)==1) {    
    if (is.null(length.E) || (length.E!=1)){
      warning("Since 'fixed.E'=VALUE, then 'length.E' must be equal to 1. Updating: 'length.E'=1.")
      length.E <- 1.0
    } 
    if ((Emax!=fixed.E) || (Emin!=fixed.E)) {
      warning("Since 'fixed.E'=VALUE, then 'Emin' and 'Emax' must be equal to 'fixed.E'. Updating: 'Emin','Emax' to be 'fixed.E'.")
      Emin <- Emax <- fixed.E
    }
    # Override with 1-point grid:
    E.grid   <- Emin # == Emax
    dE       <- 1.0    
  } else { stop("Please check 'fixed.E' variable: must be a VALUE of length=1 or a NULL.") }
  ######################################################  
  return("pars"=list("Bmin"=Bmin,"Bmax"=Bmax,"length.B"=length.B,"fixed.B"=fixed.B,"dB"=dB,"B.grid"=B.grid,
                     "Lmin"=Lmin,"Lmax"=Lmax,"length.L"=length.L,"fixed.L"=fixed.L,"dL"=dL,"L.grid"=L.grid,
                     "Emin"=Emin,"Emax"=Emax,"length.E"=length.E,"fixed.E"=fixed.E,"dE"=dE,"E.grid"=E.grid
  ))
}  

# Utility function to remove zero columns from the read data from file
"removeDataColumns" <- function(data, Lmin, Lmax, dL, length.L, L.grid,verbose=FALSE){
  # (Optional) Remove outside zero columns of Location
  idxLcol <- seq(from=(grep("all",names(data)) + 1), to=(max(grep(Lmax-dL,names(data)))),by=1) 
  idxAdj  <- max(idxLcol)-(length.L-1)
  fromTop    <- TRUE
  fromBottom <- TRUE
  idxTop    <- numeric(0)
  idxBottom <- numeric(0)
  for (j in 1:length.L){
    # Check from bottom, if found all zeros, mark for removal
    if (fromBottom) {
      idx <- idxLcol[j]
      if (sum(data[,idx])==0) {
        idxBottom <- c(idxBottom,idx)
      } else {
        fromBottom <- FALSE #Stop searching from this side
      }
    }
    # Check from top, if found all zeros, mark for removal
    if (fromTop) {       
      idx <- idxLcol[length.L-1-j+1]
      if (sum(data[,idx])==0) {
        idxTop <- c(idxTop,idx)
        endFlag <- TRUE
      } else {
        fromTop <- FALSE #Stop searching from this side
      }
    }
  }
  # if found all zeros, then remove
  idxRemove <- union(idxBottom,idxTop)
  if (length(idxRemove)>0) {
    if(verbose){
      cat("All-zero counts found in Location column(s) of CDFN table at: ",names(data)[idxRemove],". Removing the Column(s) of Location.\n")
    }
    L.grid <- L.grid[-length.L] # Temporarily remove top boundary
    L.grid <- L.grid[-(idxRemove-idxAdj)]
    Lmin  <- min(L.grid)
    Lmax  <- max(L.grid) + dL
    L.grid <- c(L.grid,Lmax) # Replace top boundary
    length.L <- length(L.grid)
    data <- data[,-idxRemove]
  }  
  
  return(list("data"=data,"Lmin"=Lmin,"Lmax"=Lmax,"length.L"=length.L,"L.grid"=L.grid))
}

########################## Set Tools functions ######################
# Define a specialized S3 function to create 'pble.table' object
"new.pble.table" <- function(effects.file) {
  
  verbose <- FALSE
  
  # Read histogram file of effects B,L,E
  data <- read.table(effects.file,header=TRUE)
  
  # Read file with scan() to retrieve all other variables
  info <- scan(file=effects.file,sep="#",what="character",nlines=7,quiet=TRUE)
  info <- info[-c(1,3,5,7,9,11)]
  
  Bmin <- abs(eval(parse(text=info[1])))
  if (!is.numeric(Bmin)) {
    stop("Could not read 'Bmin' from effects file.")
  }
  dB   <- abs(eval(parse(text=info[2])))
  if (!is.numeric(dB)) {
    stop("Could not read 'dB' from effects file.")
  }
  Bidx <- unique(data$BGidx)
  B.grid <- Bmin + c(Bidx,(max(Bidx)+1))*dB
  Bmax <- max(B.grid) 
  length.B <- length(B.grid)
  fixed.B  <- NULL  
  # Include actual Background values
  data$B <- Bmin + data$BGidx*dB  
  
  Emin <- abs(eval(parse(text=info[3])))
  if (!is.numeric(Emin)) {
    stop("Could not read 'Emin' from effects file.")
  }
  dE   <- abs(eval(parse(text=info[4])))
  if (!is.numeric(dE)) {
    stop("Could not read 'dE' from effects file.")
  }
  Eidx <- unique(data$EMAPidx)
  E.grid <- Emin + c(Eidx,(max(Eidx)+1))*dE
  Emax <- max(E.grid)
  length.E <- length(E.grid)
  fixed.E  <- NULL
  # Include actual Exposure Map values
  data$E <- Emin + data$EMAPidx*dE  
  
  Lmin <- 0
  Lmax <- 18
  dL   <- 2
  L.grid <- seq(Lmin,Lmax,by=dL)
  length.L <- length(L.grid)
  fixed.L  <- NULL 
  
  # (Optional) Remove outside zero columns of Location
  newdata  <- removeDataColumns(data=data, Lmin=Lmin, Lmax=Lmax, dL=dL, length.L=length.L, L.grid=L.grid, verbose=verbose)
  Lmin     <- newdata$Lmin
  Lmax     <- newdata$Lmax
  L.grid    <- newdata$L.grid
  length.L <- newdata$length.L
  data     <- newdata$data
  
  names.L  <- paste("L",L.grid,sep=".")
  
  # Create cdfn table in long easy-to-access format
  BLE.grid <- expand.grid(E=E.grid,L=L.grid,B=B.grid)
  BLE.grid <- BLE.grid[,c(3,2,1)]
  idxLcol <- seq(from=(grep("all",names(data)) + 1), to=(max(grep(Lmax-dL,names(data)))),by=1)  
  counts  <- matrix(data.matrix(data[,idxLcol]),ncol=1,byrow=TRUE)  
  temp    <- data.frame("L"=rep(L.grid[1:(length.L-1)],each=nrow(data)),
                        "B"=rep(data$B,times=length.L-1),
                        "E"=rep(data$E,times=length.L-1),
                        "freq"=counts)
  cdfn    <- merge(BLE.grid, temp, all=TRUE)
  cdfn$freq[is.na(cdfn$freq)] <- 0
  cdfn$prob <- cdfn$freq/sum(cdfn$freq)   #/(dB*dL*dE)
  
  return(new("pble.table",          
             "pars"=list("Bmin"=Bmin,"Bmax"=Bmax,"length.B"=length.B,"fixed.B"=fixed.B,"dB"=dB,
                         "Lmin"=Lmin,"Lmax"=Lmax,"length.L"=length.L,"fixed.L"=fixed.L,"dL"=dL,
                         "Emin"=Emin,"Emax"=Emax,"length.E"=length.E,"fixed.E"=fixed.E,"dE"=dE
             ),
             "grid"=list("B.grid"=B.grid,"L.grid"=L.grid,"E.grid"=E.grid),
             "BLE.grid"=BLE.grid,
             "table_orig"=data,
             "table"=cdfn,
             "filename"=effects.file))
}


# Define a specialized S3 function to create 'pble.basic' object
"new.pble.basic" <- function(Bmin, Bmax, Lmin, Lmax, Emin, Emax, length.B=10, length.L=2, length.E=20, fixed.B=NULL, fixed.L=NULL, fixed.E=NULL,verbose=FALSE) {
  cat("Creating pble.basic object...\n")
  # Perform basic error checking
  pars <- fixValidPars(Bmin=Bmin, Bmax=Bmax, length.B=length.B, fixed.B=fixed.B,
                       Lmin=Lmin, Lmax=Lmax, length.L=length.L, fixed.L=fixed.L, 
                       Emin=Emin, Emax=Emax, length.E=length.E, fixed.E=fixed.E)
  Bmin     <- pars$Bmin
  Bmax     <- pars$Bmax
  length.B <- pars$length.B
  fixed.B  <- pars$fixed.B
  B.grid   <- pars$B.grid
  dB       <- pars$dB
  Lmin     <- pars$Lmin
  Lmax     <- pars$Lmax
  length.L <- pars$length.L
  fixed.L  <- pars$fixed.L
  L.grid   <- pars$L.grid
  dL       <- pars$dL
  Emin     <- pars$Emin
  Emax     <- pars$Emax
  length.E <- pars$length.E
  fixed.E  <- pars$fixed.E
  E.grid   <- pars$E.grid
  dE       <- pars$dE   
  
  # Create BLE expanded grid
  # For Reimann Sum of BLE inner integral, the last value of the grid must be omitted 
  if (length(fixed.B)>0 && length(fixed.L)>0) { 
    BLE.grid <- expand.grid(E=E.grid[-length.E], L=L.grid, B=B.grid)
  } else if (length(fixed.B)>0) {
    BLE.grid <- expand.grid(E=E.grid[-length.E], L=L.grid[-length.L], B=B.grid)
  } else if (length(fixed.L)>0) {
    BLE.grid <- expand.grid(E=E.grid[-length.E], L=L.grid, B=B.grid[-length.B])
  } else {
    BLE.grid <- expand.grid(E=E.grid[-length.E], L=L.grid[-length.L], B=B.grid[-length.B])
  }
  if (verbose){
    cat("Created grid for BLE:\n")
    print(BLE.grid)
  }
  
  # 3-D uniform distribution density at B,L,E
  res <- 1
  if (is.null(fixed.B)) {
    res <- res/(Bmax-Bmin)
  }
  if (is.null(fixed.E)) {
    res <- res/(Emax-Emin)
  }
  if (is.null(fixed.L)) {
    res <- res/(Lmax-Lmin)
  }
  
  return(new("pble.basic",
             "pars"=list("Bmin"=Bmin,"Bmax"=Bmax,"length.B"=length.B,"fixed.B"=fixed.B,"dB"=dB,
                         "Lmin"=Lmin,"Lmax"=Lmax,"length.L"=length.L,"fixed.L"=fixed.L,"dL"=dL,
                         "Emin"=Emin,"Emax"=Emax,"length.E"=length.E,"fixed.E"=fixed.E,"dE"=dE
             ),
             "grid"=list("B.grid"=B.grid,"L.grid"=L.grid,"E.grid"=E.grid),
             "BLE.grid"=BLE.grid,
             "prob"=res))
}


# # 
# # # User Mode:
# # # Create new objects: basic or table 
# # if (use.effects.file) {
# cdfn <- new.pble.table(effects.file=effects.file)
# # } else {
# basic <- new.pble.basic(Bmin=Bmin, Bmax=Bmax, length.B=length.B, fixed.B=fixed.B,
#                         Lmin=Lmin, Lmax=Lmax, length.L=length.L, fixed.L=fixed.L,
#                         Emin=Emin, Emax=Emax, length.E=length.E, fixed.E=fixed.E)  
# # }
