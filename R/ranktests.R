#' Perform a rank-based one- or two-sample test using a Gaussian approximation.
#'
#' @param x First sample
#' @param y Second sample.  If default value of NULL, indicates one-sample test.
#' @param sv Type of scores, if character "ns", "vw", or "ss", or vector of scores if numeric.  Defaults to rank scores.
#'
#' @return A list of class htest
#'
#' @examples
#' ranktests(rnorm(10),rnorm(15),"ss")
#'
#' @export
#' @importFrom stats pnorm
ranktests<-function(x,y=NULL,sv=NULL){
   n<-length(y)
   m<-length(x)
   if(!is.null(y)){
      score<-fun.givescore(c(x,y),sv)
      test<-sum(score[seq(m)])
      b1<-n*m/(n+m)^2
      b2<- - b1/(n+m-1)
      vm<-(b1-b2)*sum(score^2)
   }else{
      score<-fun.givescore(abs(x),sv)
      test<-sum(score[x>0])
      vm<-sum(score^2)/4
   }
   z<-test/sqrt(vm)
   pv<-2*pnorm(-abs(z))
   out<-list(null.value=0,alternative="two-sided",method="Asymptotic rank score test",estimate=NA,data.name=NA,statistic=z,p.value=pv)
   class(out)<-"htest"
   return(out)
}
