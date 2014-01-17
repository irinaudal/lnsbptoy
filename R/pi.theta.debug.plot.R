"pi.theta.debug.plot" <- function(pi, fixed.Smin=FALSE, debug.pi=FALSE, do.pi.plot=FALSE, output.pi.plot=NULL) {
                               
  m <- length(grep("theta",names(pi)))
  
  if (debug.pi){
    par(mfrow=c(1,2))
    plot(y=pi$pi[,1],x=pi$theta,xlab="theta",ylab="pi(theta)",main="pi(theta)",type="l",lwd=1.6) ; abline(h=1.0,col="red")
    bg <- pi$grid.frame
    sub <- bg[,"theta"]>0.0 ; quilt.plot(bg[sub,"theta"],log(bg[sub,"S"]),log(bg[sub,"f"]))
    rm(bg,sub)
  } # END debug old version of pi


  if(do.pi.plot){
 
    # requires: 'graphics' package for persp(). It is already loaded automatically
    jpeg(output.pi.plot, quality=100,width=960,height=960)
    
    Smin.grid <- pi$Smin
    if(length(Smin.grid)>1){
      xx <- pi$Smin
      xy <- pi$theta.1
      xz <- matrix(pi$pi,length(xx),length(xy))
      persp(x=xx,y=xy,z=xz,ticktype = "detailed",theta = -60, phi = 15, col = "lightblue", 
            ltheta = -120, shade = 0.75, xlab="Smin",ylab="theta",zlab="pi")
   
    } else if(m==1) {
      # Just plot pi(theta) vs. theta
      plot(pi$pi~pi$grid[,grep("theta",colnames(pi$grid))],type="l",main=ppaste("pi(theta), Smin=",fixed.Smin))
      
    } else if(m==2) {
      # make three plots. 1. Perspective plot (3d), 2. pi vs. theta1, 3. pi vs. theta2.
      xx <- pi$theta.1
      xy <- pi$theta.2
      xz <- matrix(pi$pi,length(xx),length(xy))
      persp(x=xx,y=xy,z=xz,ticktype = "detailed",theta = 30, phi = 15, col = "lightblue", 
            ltheta = -120, shade = 0.75, xlab="theta1",ylab="theta2",zlab="pi",
            main=ppaste("pi(theta), Smin=",fixed.Smin))

      #   # add 2-D slices of pi, for different values of Smin
      th.val <- pi$theta.1
      plot(c(min(th.val),max(th.val)),c(0,1),type="n",xlab="theta2",ylab="pi",
           main=ppaste("pi(theta,Smin) vs. theta2 @ Smin=",fixed.Smin,
                       ",every theta1 in [",round(min(th.val),2),",",round(max(th.val),2),"]"))
      for(k in 1:length(th.val)){
        theta.id <- pi$grid[,"theta.1"]==th.val[k]
        lines(pi$grid[theta.id,2],pi$pi[theta.id])
      }
      
      th.val <- pi$theta.2
      plot(c(min(th.val),max(th.val)),c(0,1),type="n",xlab="theta1",ylab="pi",
           main=ppaste("pi(theta,Smin) vs. theta1 @ Smin=",fixed.Smin,",every theta2 in [",round(min(th.val),2),",",round(max(th.val),2),"]"))
      for(k in 1:length(th.val)){
        theta.id <- pi$grid[,"theta.2"]==th.val[k]
        lines(pi$grid[theta.id,1],pi$pi[theta.id])
      }
    } # END plot variations of pi(theta) vs. theta
  
    dev.off()
    
  } # END do.pi.plot
  
  return(NULL)
}