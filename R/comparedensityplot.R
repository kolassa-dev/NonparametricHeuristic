#' Plot various distributions: Cauchy, normal, and uniform.
#'
#' @examples
#' fun.comparedensityplot()
#'
#' @export
#' @importFrom graphics plot
#' @importFrom stats dt dnorm qnorm dunif
fun.comparedensityplot<-function(){
   ordinate<- (-400:400)/100
   densout<-array(NA,c(length(ordinate),3))
   densout[,1]<-dt(ordinate,1)
   densout[,2]<-dnorm(ordinate,sd=1/qnorm(.75))
   densout[,3]<-dunif(ordinate,min=-2,max=2)
   plot(range(ordinate),range(densout),
      main="Fig. 1: Comparison of Three Densities",
      sub="Distributions scaled to have IQR 2",type="n",
      xlab="ordinate",ylab="density")
   legend(-1.3,.11,lty=1:3,legend=c("Cauchy","Normal","Uniform"))
   for(j in 1:3) lines(ordinate,densout[,j],lty=j)
}
