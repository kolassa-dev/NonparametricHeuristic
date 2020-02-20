#' Construct confidence interval for regression slope parameter by inverting correlation tests.
#'
#' @param xxx Explanatory vector
#' @param yyy Response vector
#'
#' @return A matrix of coverage levels.  Dimensions are distribution and confidence interval technique.
#'
#' @return Matrix of confidence intervals
#'
#' @examples
#' invertcor(rnorm(10),rnorm(10))
#'
#' @export
#' @importFrom quantreg rq
#' @importFrom stats cor
invertcor<-function(xxx,yyy){
   regp<-(-200:800)/400
   out<-array(NA,c(length(regp),3))
   for(j in seq(length(regp))){
      out[j,1]<-cor(xxx,yyy-regp[j]*xxx)
      out[j,2]<-cor(xxx,yyy-regp[j]*xxx,method="spearman")
      out[j,3]<-cor(xxx,yyy-regp[j]*xxx,method="kendall")
   }
   ci<-array(NA,c(3,2))
   dimnames(ci)<-list(c("Pearson","Spearman","Kendall"), c("Lower","Upper"))
   n<-length(xxx)
   sds<-sqrt(c(1/(n-1),1/(n-1),2*(2*n+5)/(9*n*(n-1))))
   for(i in 1:3) 
      ci[i,]<-approx(out[,i],regp,xout=c(1,-1)*1.96*sds[i])$y
   plot(range(regp),range(out),type="n",
      xlab="Regression Parameter",ylab="Statistic Value",
      main="Construction of Regression Estimates for BP Data",
      sub="Statistic Quantiles are approximated using Normal")
   for(i in 1:3){
      lines(regp,out[,i],lty=i)
      for(k in c(-1,1)) abline(h=k*1.96*sds[i],lty=i)
      for(k in 1:2) abline(v=ci[i,k],lty=i)    
   }
   legend(1.8,.85,legend=dimnames(ci)[[1]],lty=1:3)
   return(ci)
}
