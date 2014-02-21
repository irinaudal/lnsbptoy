
"winterpolate" <- function(x,Z,new.x,dist.p=2,verbose=FALSE)
{

# x :: list of length p, each component being a vector
# Z :: p-dimensional array
# new.x :: matrix of dimension [n x p]
# Require that: dim(Z) = sapply(x,length)

    if (!is.list(x)){
        stop("'x' must be a list")
    }
    new.x <- as.matrix(new.x)  # TODO: Is this a good idea to do this? I ask because you cannot coerce elements to "double" the object which is a list
    ##> storage.mode(new.x) <- "double"
    ##Error in storage.mode(new.x) <- "double" : 
    ##  (list) object cannot be coerced to type 'double'
    
    if (!is.matrix(new.x)){
        stop("'new.x' must be a matrix")
    }
    if (!is.array(Z)){
        stop("'Z' must be an array")
    }
    p <- length(x)
    dims <- sapply(x,length)
    if (any(dim(Z) != dims)){
        cat("Error:\n")
        cat("dim(Z):\n") ; print(dim(Z))
        cat("Grid dimensions:\n") ; print(dims)
        stop("'Z' must have dimension matching 'x'")
    }
    if (ncol(new.x) != p){
        stop("'new.x' must have length equal to the length of 'x'")
    }
    for (i in 1:p){
        storage.mode(x[[i]]) <- "double" 
    }
    storage.mode(new.x) <- "double" 
    storage.mode(Z) <- "double"

   bin.ix <- as.matrix(expand.grid(rep(list(c(0,1)),p)))
   storage.mode(bin.ix) <- "integer"
   dist.p <- as.double(dist.p)

   winterp <- .Call("winterpolate",
                   PACKAGE="logNlogS",
                list("verbose"= as.integer(verbose),
                     "x"=x,"Z"=Z,"new"=new.x,"dist_p"=dist.p,"bin_ix"=bin.ix))

  return(winterp)

}


