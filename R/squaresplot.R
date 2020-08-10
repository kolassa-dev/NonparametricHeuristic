#' Draw contours of OLS regression as a function of the two parameters.
#'
#' @param x explanatory variable
#' @param y response variable
#' @param centerx Logical variable determining whether to center covariate
#' @param npts number of lattice points on each margin.
#'
#' @importFrom graphics plot.default contour.default points.default
#' @export
squaresplot<-function(x,y,centerx=FALSE,npts=100){
   if(centerx) x<-x-mean(x)
   cc<-summary(lm(y~x))$coefficients
   lims<-outer(cc[,1],c(1,1))+outer(cc[,2],c(-3,3))
   plot.default(lims[1,],lims[2,],type="n",
      xlab=dimnames(lims)[[1]][1], ylab=dimnames(lims)[[1]][2])
   cont<-array(NA,c(npts,npts))
   bvl<-vector(mode="list",2)
   for(kk in 1:2) bvl[[kk]]<-lims[kk,]%*%rbind(
      1-(0:(npts-1))/(npts-1),(0:(npts-1))/(npts-1))
   for(ii in 1:npts) for(jj in 1:npts)
      cont[ii,jj]<-mean((y-bvl[[1]][ii]-bvl[[2]][jj]*x)^2)
   contour.default(bvl[[1]],bvl[[2]],cont,add=TRUE)
   points.default(cc[1,1],cc[2,1])
}
