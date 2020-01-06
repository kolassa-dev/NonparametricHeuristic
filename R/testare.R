#' Check whether powers are the same for sample sizes in the ratio indicated by asymptotic relative efficiency are similar, for sign and t tests, for Laplace variables.
#'
#' @param testare Putative relative efficiency, taking a value between 0 and 1.  Set to .5 for Laplace and t to sign test comparison.
#' @param nvec vector of larger sample sizes to try.
#' @param nsamp number of Monte Carlo samples.
#' @return A table of powers.
#'
#' @examples
#' testare(.5)
#'
#' @importFrom stats t.test binom.test
#' @importFrom VGAM rlaplace
#' @export
testare<-function(targetare=1,nvec=c(20,100,1000,10000), nsamp=10000){
   pw<-array(0,c(2,length(nvec)))
   dimnames(pw)<-list(c("t test","sign test"),as.character(nvec))
   for(k in seq(length(nvec))){
      for(j in seq(nsamp)){
         x<-rlaplace(nvec[k]*targetare)+3/sqrt(nvec[k])
         y<-rlaplace(nvec[k]*(1-targetare))+3/sqrt(nvec[k])
         pw[1,k]<-pw[1,k]+(t.test(c(x,y))$p.value<0.05)
         pw[2,k]<-pw[2,k]+(binom.test(sum(x>0),nvec[k]/2)$p.value<0.05)
      }
   }
   return(t(pw/nsamp))
}
