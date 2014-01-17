#"ppaste" <- function(...){paste(...,sep="")}

"scp.copy" <- function(from.filename,to.filename,from.delete=FALSE,recursive.dir.create=FALSE,verbose=TRUE)
{
  ####
  # usage:
  # scp.copy("foo.txt","user@server.com:~/foo.txt")
  #
  # Notes:
  # 
  # This function cannot currently handle situations where the
  # destination filename is different from from.filename AND 
  # recursive directory creation is required.
  # 
  # It is possible to change the filename in the copy,
  # but in this situation the directory must already exist.
  # When retaining the same filename, directories can be
  # recursively computed but to.filename should specify the 
  # directory only -- not the full filepath.
  # 
  ####
  
  if (recursive.dir.create){
    copy.file.cmd <- paste("rsync -ave ssh ",from.filename," ",to.filename,sep="")
  } else {
    copy.file.cmd <- paste("scp ",from.filename," ",to.filename,sep="")
  }
  if (verbose){
    cat("Copying local image file...\n")
    cat(paste("Recursive directory creation = ",recursive.dir.create,"\n",sep=""))
    cat(paste("command: ",copy.file.cmd,"\n",sep=""))
  }
  system(copy.file.cmd)
  if (verbose){
    cat("Copy completed.\n")
  }
  if (from.delete){
    delete.file.cmd <- paste("rm ",from.filename,sep="")
    if (verbose){
      cat("Removing local copy...\n")
      cat(paste("command: ",delete.file.cmd,"\n",sep=""))
    }
    system(delete.file.cmd)
    if (verbose){
      cat("Delete completed.\n")
    }
  }
  return()
}

"truehist.fix" <- function(nbins=40,...)
{
  ##
  # Fix for truehist:
  # If there is only a single value then truehist runs fine 
  # when bins is specified, but breaks when it is not.
  # The following fix should work...
  ##
  take.one <- try(truehist(...),silent=TRUE)
  if (class(take.one)=="try-error")
    truehist(...,nbins=nbins)
  return()
}

"expand.grid.alt" <- function(seq1,seq2,seq3=NULL,seq4=NULL,seq5=NULL,seq6=NULL)
{
  ####
  # Fast way to expand grid of multiple vectors through the use of rep.int() function
  # Usage: expand.grid.alt(x,y,z)
  ####
  
  nargin <- nargs()
  x <- cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
  if (nargin>=3){
    temp <- x
    x <- cbind(Var1=rep.int(temp[,1], length(seq3)), 
               Var2=rep.int(temp[,2], length(seq3)), 
               Var3=rep(seq3, each=nrow(temp)))
  }
  if (nargin>=4){
    temp <- x
    x <- cbind(Var1=rep.int(temp[,1], length(seq4)), 
               Var2=rep.int(temp[,2], length(seq4)), 
               Var3=rep.int(temp[,3], length(seq4)), 
               Var4=rep(seq4, each=nrow(temp)))
  } 
  if (nargin>=5){
    temp <- x
    x <- cbind(Var1=rep.int(temp[,1], length(seq5)), 
               Var2=rep.int(temp[,2], length(seq5)), 
               Var3=rep.int(temp[,3], length(seq5)), 
               Var4=rep.int(temp[,4], length(seq5)), 
               Var5=rep(seq5, each=nrow(temp)))
  } 
  if (nargin>=6){
    temp <- x
    x <- cbind(Var1=rep.int(temp[,1], length(seq6)), 
               Var2=rep.int(temp[,2], length(seq6)), 
               Var3=rep.int(temp[,3], length(seq6)), 
               Var4=rep.int(temp[,4], length(seq6)), 
               Var5=rep.int(temp[,5], length(seq6)), 
               Var6=rep(seq6, each=nrow(temp)))
  }
  return(as.data.frame(x))
}

"showMemoryUse" <- function(sort="size", decreasing=FALSE, limit,pos=1) 
{
  ####
  # Script for improved list of objects: create sorted list of objects, sorted by memory
  # NOTE: use this code as script to be called within some function environment, otherwise, environment of shorMemoryUse() is taken!
  ####
  
  if(pos==1){
    objectList <- ls(parent.frame())
  } else if(pos==-1){
    objectList <- ls(pos=pos)
  }
  
  oneKB <- 1024
  oneMB <- 1048576
  oneGB <- 1073741824
  
  memoryUse <- sapply(objectList, function(x) as.numeric(object.size(eval(parse(text=x)))))
  
  memListing <- sapply(memoryUse, function(size) {
    if (size >= oneGB) return(paste(round(size/oneGB,2), "GB"))
    else if (size >= oneMB) return(paste(round(size/oneMB,2), "MB"))
    else if (size >= oneKB) return(paste(round(size/oneKB,2), "kB"))
    else return(paste(size, "bytes"))
  })
  
  memListing <- data.frame(objectName=names(memListing),memorySize=memListing,row.names=NULL)
  
  if (sort=="alphabetical") memListing <- memListing[order(memListing$objectName,decreasing=decreasing),] 
  else memListing <- memListing[order(memoryUse,decreasing=decreasing),] #will run if sort not specified or "size"
  
  if(!missing(limit)) memListing <- memListing[1:limit,]
  
  print(memListing, row.names=FALSE)
  return(invisible(memListing))
}

"shade.area.under.curve" <- function(draws,q_lo,q_hi,ymin=0,density=0.5,col="lightblue",border=NULL,nbins=2000)
{
  ##
  # Shade in an area under a posterior density curve on an existing plot,
  # also, throw in for free the integral (by trapeziodal approximation)
  ##
  if (nbins<2)
    stop("'nbins' must be at least 2")
  
  dens <- density(draws, n=nbins)
  ci <- summary(draws,quantiles=c(q_lo,q_hi))$quantiles
  ci_y <- approx(dens, xout = ci)$y
  lo <- ci[1]
  hi <- ci[2]
  
  id <- dens$x >= lo & dens$x <= hi
  bc <- c(lo, dens$x[id], hi)
  yv <- c(ci_y[1], dens$y[id], ci_y[2])
  
  for (i in 1:nbins){
    polygon(list("x"=bc[c(i,i,i+1,i+1)],"y"=c(ymin,ymin,yv[i],yv[i+1])),density=density,col=col,border=border)
  }
  
  # Throw in a simple trapezoid approximation to the area as well:
  delta.x <- diff(bc[1:2])
  #area <- sum(yv)*delta.x
  area <- (sum(0.5*2*(yv-ymin))-sum(yv[c(1,length(yv))]-ymin))*delta.x
  return(area)
}