"g.table" <- function(incompleteness.file, inc.factor)
{
  
  ####################################################################################################
  # g.table   Function to read the tabulated survey imcompleteness probabilities, g=P(observing a source) 
  #
  # NOTE: user must also specify in the test files that
  #             use.incompleteness.file <- TRUE
  #             eg. incompleteness.file     <- "Incompleteness_Functions/g_full_table.txt"
  #
  # Input: incompleteness.file = name of incompleteness file which summarizes the tabulated g
  #        inc.factor          = incompleteness convergence factor (convert b/w photons and counts)
  #
  # Output: g.tbl2  = (function of lambda(S),B,L,E) probability of observing source in a form of a list:
  #                     table = tabulated values of g,  
  #                     pars  = list of L, length.L, dL, bg, length.bg, dbg, lambda, length.lambda, dlambda,
  #                     g.type = "table"
  ####################################################################################################
  
  d <- read.table(incompleteness.file,sep="\n",as.is=TRUE,comment.char="%")
  
  # Find specific information on variables
  comment.row <- grep("#",d[,1])
  
  off.row <- grep("#off-axis",d[comment.row,1])
  bg.row <- grep("#background",d[comment.row,1])
  src.row <- grep("#source",d[comment.row,1])
  table.row <- grep("Table format",d[comment.row,1])
  header.row <- comment.row[c(bg.row-1, src.row-1, table.row-1)]
  
  # off-axis angle = location: [degrees] = 60x[arcmin]   #TODO: check about this unit with Andreas
  n.rows   <- comment.row[bg.row]-header.row[1]-1
  off.orig <- read.table(incompleteness.file,skip=header.row[1],header=FALSE,nrows=n.rows)[,1]
  off      <- off.orig # *60
  
  # background: [ph/s/cm^2/pixel] = inc.factor x [ct/pixel] 
  n.rows  <- comment.row[src.row]-header.row[2]-1
  bg.orig <- read.table(incompleteness.file,skip=header.row[2],header=FALSE,nrows=n.rows)[,1]
  
  # source counts: [count / s]   #TODO: check about this unit with Andreas
  n.rows   <- comment.row[table.row]-header.row[3]-1
  src.orig <- read.table(incompleteness.file,skip=header.row[3],header=FALSE,nrows=n.rows)[,1]
  src.ct   <- src.orig
  
  # Check that src cts begin with 0
  flag1  <- FALSE
  if (all(src.ct>0)){
    flag1 <- TRUE
    src.ct <- c(0,src.ct)
  }
  # Check that bg begin with 0
  flag2  <- FALSE
  if (all(bg.orig>0)){
    flag2 <- TRUE
  }
  # Check that off axis angle locations begin with 0
  flag3  <- FALSE
  if (all(off>0)){
    flag3 <- TRUE
    off <- c(0,off)
  }
  
  # Add upper boundary to every measurement - required for interpolation of g later
  flag1.5 <- TRUE
  if (flag1.5){
    src.ct <- c(src.ct,Inf)    
  }
  
  # Read g-tables as is
  header.row <- grep("#off=",d[,1])
  n.rows     <- header.row[2]-header.row[1]-1
  n <- length(header.row)
  data <- list(numeric(n))
  
  for (i in 1:n) {
    temp <- read.table(incompleteness.file,skip=header.row[i],header=FALSE,nrows=n.rows)
    
    if (flag1){ # Since we needed to add zero value to src.cts, add zero prob as well
      if (all(temp==temp[1,1])) { # check fake version of g.table, where all values are the same constant.
        temp <- rbind(temp[1,1], temp)
      } else {
        temp <- rbind(0, temp)
      }
    }
    
    if (flag2){ # Add upper boundary to prob
      if (all(temp==temp[1,1])) { # check fake version of g.table, where all values are the same constant.
        temp <- rbind(temp, temp[1,1])
      } else {
        temp <- rbind(temp,1)
      }
    }
    
    # Force each column to be non-decreasing
    for (j in 1:length(bg.orig)){
      while (is.unsorted(temp[,j])){
        idx <- which(diff(temp[,j])<0)[1]
        # Replace bad value by previous good entry
        temp[idx+1,j] <- temp[idx,j]
      }
    }
    
    # Add lower boundary prob. for lower bg and lower offaxis angle
    if (flag2) {
      temp <- cbind(temp[,1], temp) # copy over the first column of g based on lowest bg
    }
    if (flag3) {
      stop("Incompleteness tables must specify g for OffAxis angle = 0.")
    }
    
    data[[i]] <- temp
  }
  
  if (flag2) { # add bg begin at 0
    bg.orig <- c(0,bg.orig)
  }
  
  # Enter all information into long format table
  # To convert to count/pixel multiply the values below by  0.43E3*10000 = 430 pixels * 10000 sec ?
  #   (this is estimated for power-law with Gamma=1.7 NH=6.0E20 cm^-2);  
  #   10000s is the exposure used in the simulations 
  bg        <- bg.orig*(inc.factor)
  main.grid <- expand.grid(lambda=src.ct, bg=bg, L=off)
  main.grid <- main.grid[,3:1]
  
  g.tbl <- main.grid
  g.tbl$g <- NA
  
  for (j in 1:length(data)){ 
    #  temp <- reshape(data[[j]], direction="long", varying=names(data[[j]]), v.names="g", timevar="bg",new.row.var)
    g.tbl$g[main.grid$L==off[j]] <- unlist(data[[j]])
  }
  g.tbl2 <- list("table"=g.tbl, 
                 "pars"=list("L"=off,"length.L"=length(off),"dL"=off[2]-off[1],
                             "bg"=bg,"length.bg"=length(bg),"dbg"=bg[2]-bg[1],
                             "lambda"=src.ct,"length.lambda"=length(src.ct),"dlambda"=src.ct[3]-src.ct[2]),
                 "g.type"="table")
  
  return(g.tbl2)
}
