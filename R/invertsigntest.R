#' Diagram for the inversion of rank tests to get confidence intervals for a quantile, or, in the two-sample case, location differences.
#'
#' @param dataset Data vector, or list of two vectors representing two samples.
#' @param alpha 1-confidence value; defaults to 0.05
#' @param maint title for plot.
#' @param tau quantile to be estimated.  Defaults to 0.5.
#' @param returndata Logical flag determining whether the data set is included as part of the response.
#' @param npts number of points at which to evaluate test prior to inversion.
#' @param hl logical flag to determine whether to use sign or Mood test (FALSE) or the signed rank or rank sum test (TRUE)
#' @param rescaley logical flag to determine whether to rescale the vertical scale to [0,1].
#' @param newplot logical flag to start a new plot
#' @examples
#' invertsigntest(rnorm(12))
#' @export
#' @importFrom stats qwilcox qbinom
#' @importFrom graphics axis
invertsigntest<-function(dataset,alpha=0.05,maint=NULL,tau=0.5,returndata=FALSE,npts=10000,hl=NULL,rescaley=FALSE,newplot=TRUE){
   twosamp<-!is.atomic(dataset)
   if(is.null(hl)) hl<-twosamp
   if(twosamp){
      mm<-length(dataset[[1]])
      nn<-length(dataset[[2]])
      dataset<-as.vector(outer(dataset[[2]],dataset[[1]],"-"))
#     browser()
   }else{
      nn<-length(dataset)
      if(hl){
         temp<-outer(dataset,dataset,"+")/2
         dataset<-temp[!upper.tri(temp)]
      }
   }
   dataset<-sort(dataset)
   rd<-range(dataset)
   mu<-rd[1]-.2*diff(rd)+ (0:npts)/npts*1.4*diff(rd)
   statval<-apply(outer(dataset,mu,"-")<0,2,sum)
   if(rescaley){
      mx<-max(statval)
      statval<-statval/mx
   }else{
      mx<-1
   }
   if(twosamp){
      if(hl){
         a<-qwilcox(alpha/2,mm,nn)
         b<-mm*nn+1-qwilcox(alpha/2,mm,nn)
         message("a=",a," b=",b)
      }else{
         message("Confidence interval construction from the Mood Median Test not yet implemented.")
         a<-b<-NA
      }
   }else{
      if(hl){
         a<-qsignrank(alpha/2,nn)
         b<-nn*(nn+1)/2-a
      }else{
         a<-qbinom(alpha/2,length(dataset),tau)
         b<-length(dataset)+1-qbinom(alpha/2,length(dataset),1-tau)
      }
   }
   maintt<-if(twosamp) "Median Difference" else paste("Quantile",tau,sep=" ")
   if(newplot){
      plot(rd+c(-1,1)*.2*diff(rd),range(statval),type="n",
         ylab=if(twosamp) "Statistic" else "Statistic",
         xlab=if(twosamp) "Location Difference" else "Location",
         sub=paste("Confidence",1-alpha,if(twosamp) "Hodges-Lehman" else "Median",
            "Interval uses order statistics",a,"and",b,sep=" "),
         main=paste("Construction of CI for",maint,maintt,sep=" "))
   }
   temp<-c(rd[1]-.2*diff(rd),dataset,rd[2]+.2*diff(rd))
   for(jj in seq(length(temp)-1)){
      sv<-sum(dataset<=temp[jj])
      segments(temp[jj],sv,temp[jj+1],sv)
   }
#  lines(mu,statval,col=2)
   abline(h=a/mx,lty=2)
   abline(h=(b-1)/mx,lty=2)
   axis(4,at=c(a,b-1),labels=c("a","b-1"))
# dataset is already sorted
   ci<-dataset[c(a,b)]
   for(i in 1:2) abline(v=ci[i],lty=3)
   legend(ci[2]+.01*diff(range(dataset)),b-1/2,
      lty=c(2,3),legend=c("Critical Values","CI Endpoints"))
   out<-if(returndata) list(ci=ci,dataset=dataset) else ci
   return(out)
}
