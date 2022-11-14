#' Calculate and plot the density of the studentized cauchy means, via simulation.
#'
#' @param nobs Sample size for studentized means.
#' @param nsamp Number of Monte Carlo samples.  Defaults to 10000.
#' @return A five-dimensional array of test levels (if altevalue equals hypoth) or power otherwise.  Dimensions are test, distribution, alternative, alternative value, and first dimension sample size.
#'
#' @examples
#' fun.studentizedcaucyplot(10,nsamp=100)
#'
#' @export
#' @importFrom stats rcauchy sd density median
fun.studentizedcaucyplot<-function(nobs,nsamp=10000){
  densities<-vector("list",length(nobs))
  ry<-NULL
  u<-array(rcauchy(max(nobs)*nsamp),c(max(nobs),nsamp))
  for(j in seq(length(nobs))){
     xbar<-apply(u[seq(nobs[j]),],2,mean)
     s<-apply(u[seq(nobs[j]),],2,sd)
     t<-xbar/(s/sqrt(nobs[j]))
     densities[[j]]<-density(c(t,-t))
     ry<-range(c(ry,densities[[j]]$y))
  }
  plot(range(densities[[1]]$x),ry,type="n",xlab="Data Value",
     main="Symmetrized Density of Studentized Cauchy")
  for(j in seq(length(nobs))) lines(densities[[j]],lty=j)
  legend(0,median(ry),lty=seq(length(nobs)),legend=nobs)
}
