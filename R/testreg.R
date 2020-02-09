#' Compare coverage of confidence intervals from various regression techniques.  Intervals are nominal coverage .1.
#'
#' @param dists vector of names of random number generators for sample.  Defaults to normal.
#' @param npergp Number of observations per sample; defaults to 10
#' @param mcsamp Number of Monte Carlo samples.
#'
#' @return A matrix of coverage levels.  Dimensions are distribution and confidence interval technique.
#'
#' @return A list with two components.  The first is a two-dimensional array of average interval lengths.  The second is coverage levels.  Dimensions of both are technique and distribution.
#'
#' @examples
#' library(VGAM)
#' fun.testreg(mcsamp=100,dists=c("rnorm","rcauchy","rlaplace","runif","rexp"))
#'
#' @export
#' @importFrom quantreg rq
fun.testreg<-function(dists="rnorm",npergp=10,mcsamp=10000){
   out<-array(NA,c(mcsamp,2,length(dists)))
   ciends<-array(NA,c(mcsamp,2,2,length(dists)))
   dimnames(out)<-list(NULL,c("L1","L2"),dists)
   dimnames(ciends)<-list(NULL,c("L1","L2"),c("l","u"),dists)
   for(k in seq(length(dists))){
      thisdist<-eval(parse(text=dists[k]))
      for(j in seq(dim(out)[1])){
         x<-seq(npergp)
         y<-thisdist(npergp)
         if(dists[k]=="runif") y<-y-.5
         if(dists[k]=="rexp") y<-y-1
# Confidence intervals for rq are calculated by rq.fit.br inside summary.rq .
# Pass this as noted in the summary.rq documentation by alpha= in the ...
# argument.
         cc<-summary(rq(y~x))$coefficients[2,2:3]
         ciends[j,1,,k]<-cc
         out[j,1,k]<-(cc[1]<0)&(cc[2]>0)
         ss<-summary(lm(y~x))$coefficients[2,]
         out[j,2,k]<-ss[4]>.1
         ciends[j,2,,k]<-ss[1]+c(-1,1)*qt(.95,8)*ss[2]
      }
   }
   ciendav<-apply(ciends,2:4,mean)
   return(list(cover=apply(out,2:3,mean),cil=ciendav[,2,]-ciendav[,1,]))
}
