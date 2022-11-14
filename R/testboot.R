#' Compare coverage of bootstrap confidence intervals.  
#'
#' @param fcn function to be bootstrapped.
#' @param sampsize Number of observations per sample; defaults to 10
#' @param mcsamp Number of Monte Carlo samples.
#' @param nsamp Deprecated version of mcsamp.
#' @param B Bootstrap size
#' @param dists vector of names of random number generators for sample.  Defaults to normal and Cauchy.  This fails for a single distribution.
#' @param truem True parameter value; should be of length 1 or length of distribution vector.
#' @param alpha 1-confidence level.
#'
#' @return A list with two components.  The first is a two-dimensional array of average interval lengths.  The second is coverage levels.  Dimensions of both are distribution and confidence interval technique.
#'
#' @examples
#' fun.testboot(function(x,indices) return(c(mean(x[indices]),1)),mcsamp=100,B=999)
#'
#' @export
#' @importFrom boot boot boot.ci
fun.testboot<-function(fcn,sampsize=10,nsamp=NA,mcsamp=10000,B=9999,dists=c("rnorm","rcauchy"),truem=NULL,alpha=.05){
  out<-NULL
  ok<-NULL
  if(is.null(truem)) truem<-rep(0,length(dists))
  if(!is.na(nsamp)){mcsamp=nsamp}
  for(j in seq(length(dists))){
     for(k in seq(mcsamp)){
        x<-eval(parse(text=dists[j]))(sampsize)
        ciout<-boot.ci(boot(x,fcn,B),conf=1-alpha)
        nci<-length(ciout)-3
        if(is.null(out)) out<-array(NA,c(length(dists),mcsamp,nci,2))
        if(nci==4) out[j,k,,]<-rbind(ciout$normal[2:3], ciout$basic[4:5], ciout$percent[4:5], ciout$bca[4:5])
        if(nci==5) out[j,k,,]<-rbind(ciout$normal[2:3], ciout$basic[4:5], ciout$percent[4:5], ciout$bca[4:5], ciout$student[4:5])
     }
  }
  ok<-array(TRUE,dim(out)[1:3])
  for(j in seq(length(dists))){
# R version 3.4.4 fails here if the first dimension of ok and out are 1.  Use multiple distributions.
     ok[j,,]<-(out[j,,,2]>=truem[j])&(out[j,,,1]<=truem[j])
  }
  cilen<-apply(out[,,,2]-out[,,,1],c(1,3),"mean")
  ok<-apply(ok,c(1,3),"mean")
  dimnames(ok)<-list(dists,
     c("Normal","Basic","Percentile","BCa","Studentized")[seq(nci)])
  dimnames(cilen)<-dimnames(ok)
  return(list(cilen=cilen,ok=ok))
}
