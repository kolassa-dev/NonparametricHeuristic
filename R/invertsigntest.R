#' Diagram for the inversion of the sign test to get confidence intervals for a quantile, or, in the two-sample case, the Hodges-Lehman estimator.
#'
#' @param dataset Data vector, or list of two vectors representing two samples.
#' @param alpha 1-confidence value; defaults to 0.05
#' @param maint title for plot.
#' @param tau quantile to be estimated.  Defaults to 0.5.
#' @param returndata Logical flag determining whether the data set is included as part of the response.
#' @param npts number of points at which to evaluate test prior to inversion.
#' @examples
#' invertsigntest(rnorm(12))
#' @export
#' @importFrom stats qwilcox qbinom
#' @importFrom graphics axis
invertsigntest<-function(dataset,alpha=0.05,maint=NULL,tau=0.5,returndata=FALSE,npts=10000){
   hl<-!is.atomic(dataset)
   if(hl){
      mm<-length(dataset[[1]])
      nn<-length(dataset[[2]])
      dataset<-as.vector(outer(dataset[[1]],dataset[[2]],"-"))
#     browser()
   }else{
      nn<-length(dataset)
   }
   rd<-range(dataset)
   mu<-rd[1]-.2*diff(rd)+ (0:npts)/npts*1.4*diff(rd)
   statval<-apply(outer(dataset,mu,"-")>0,2,sum)
   if(hl){
      a<-qwilcox(alpha/2,mm,nn)
      b<-mm*nn+1-qwilcox(alpha/2,mm,nn)
      message("a=",a," b=",b)
   }else{
      a<-qbinom(alpha/2,length(dataset),tau)
      b<-length(dataset)+1-qbinom(alpha/2,length(dataset),1-tau)
   }
   maintt<-if(hl) "Median Difference" else paste("Quantile",tau,sep=" ")
   plot(mu,statval,type="l",
      sub=paste("Confidence",1-alpha,if(hl) "Hodges-Lehman" else "Median",
         "Interval uses order statistics",a,"and",b,sep=" "),
      main=paste("Construction of CI for",maint,maintt,sep=" "))
   abline(h=a,lty=2)
   abline(h=b-1,lty=2)
   axis(4,at=c(a,b-1),labels=c("a","b-1"))
   ci<-sort(dataset)[c(a,b)]
   for(i in 1:2) abline(v=ci[i],lty=3)
   legend(ci[2]+.01*diff(range(dataset)),b-1/2,
      lty=c(2,3),legend=c("Critical Values","CI Endpoints"))
   out<-if(returndata) list(ci=ci,dataset=dataset) else ci
   return(out)
}
