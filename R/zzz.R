".onLoad" <- function (libname, pkgname)
{
  ## figure out this year automatically 
  this.year <- substr(as.character(Sys.Date( )), 1, 4)
  
  ## echo output to screen
  packageStartupMessage("##\n## logN-logS Analysis\n")
  packageStartupMessage("## Copyright 2011-", this.year, " Paul Baines, Irina Udaltsova\n\n", sep="")
  
  ## Dependencies now handled by NAMESPACE
  #    
  #    quiet <- TRUE
  #    require(MCMCpack, quietly=quiet)
  #    require(fields, quietly=quiet)
  #    require(compiler, quietly=quiet)
  #    require(ggplot2, quietly=quiet)
  #    require(MASS, quietly=quiet)
  #    ...
  
  # is this really needed, or does NAMESPACE take care of it now?
  # Uncomment if compiled code is added:
  # library.dynam(pkgname, pkgname, lib.loc=libname)
}




