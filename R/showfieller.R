#' Diagram the construction of Fieller's confidence interval, for the ratio of expectations of barw and bary.
#'
#' @param barw Denominator of method of moments estimator.
#' @param bary Numerator of method of moments estimator.
#' @param tau Standard deviation of Denominator
#' @param sigma Standard deviation of Numerator
#' @param rho Correlation of numerator and denominator
#' @param alpha test level
#' @param npts number of points for plot.
#'
#' @examples
#' # Standard case with confidence region as an interval.
#' showfieller(1, 1, 0.70, 0.70)
#' # Problematic case with confidence region as complement of interval
#' showfieller(1, 1, 0.20, 0.20)
#' # Problematic case with confidence region as entire real line.
#' showfieller(1, 1, 0.75, 0.75)
#'
#' @importFrom graphics plot.default lines.default
#' @export
showfieller<-function(barw,bary,tau,sigma,rho=0,alpha=0.05,npts=60){
   xihat<-barw/bary
   se<-sqrt(barw^(-2)*(sigma^2-2*rho*sigma*tau*xihat+tau^2*xihat^2))
   xi<-xihat-3.0*se+(0:npts)*se/(npts/6)
   fr<-(xi*barw-bary)^2/(xi^2*tau^2+sigma^2+2*rho*tau*sigma*xi)
   subt<-paste("bar W=",barw,"bar Y=",bary,"tau=",tau,"sigma=",sigma,
      "rho=",rho,"alpha=",alpha)
   z<-qnorm(1-alpha/2)
   plot.default(range(xi),range(c(fr,z^2)),type="n",
      main="Pivot whose distribution is used for ratio CI",
      xlab="Potential Ratio Value",ylab="Pivot",sub=subt)
   abline(h=z^2)
   lines.default(xi,fr)
}
showfieller(1, 1, 0.70, 0.70)
showfieller(1, 1, 0.20, 0.20)
showfieller(1, 1, 0.75, 0.75)
