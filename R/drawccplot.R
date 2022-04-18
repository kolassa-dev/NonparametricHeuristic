#' Exhibit construction of continuity correction when approximating a 
#' distributoin supported on a lattice by a continuous distribution.
#' Example is approximation of binomial by normal.
#' @param n binomail samle size
#' @param z ordinate for which P[Z<=z] is required.
#' @export
#' @examples
#' drawccplot()
#' @importFrom graphics barplot polygon
#' @importFrom stats dbinom
drawccplot<-function(n=10,z=4){
   x<-0:n
   y<-dbinom(x,10,.5)
   names(y)<-as.character(x)
   xx<-(0:(10*n))/10
   yy<-dnorm(xx,.5*n,sqrt(.25*n))
   barplot(y,space=0,width=1,ylim=range(yy),density=c(rep(10,z+1),rep(0,n-z)))
   lines(xx+.5,yy)
   mm<-round(xx)<=z
   polygon(c(xx[mm],rev(xx[mm]))+.5, c(yy[mm],rep(0,sum(mm))),
       density=10,angle=-45)
   legend(0,max(yy),density=10,angle=c(45,-45),legend=c("True","Approximate"))
}
