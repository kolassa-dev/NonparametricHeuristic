#' Compare power of various one- and two-sample tests via simulation.
#'
#' @param samp1sz the size of the first sample.
#' @param samp2sz the size of the second sample.  Defaults to NA, triggering one-sample tests.
#' @param dist vector of names of random number generators, determinating distributions to be compared.  Defaults to just normal.
#' @param ndim Dimension of distributions to be compared.  Defaults to one.
#' @param distb vector of names of random number generators for second sample.  Defaults to distribution for the first sample.
#' @param level level of tests.
#' @param hypoth vector null hypothesis values.  For uniform on [0,1] this should be taken as .5.
#' @param nsamp Number of Monte Carlo samples.
#' @param alternative alternative, of the form aceptable by t.test, z.test, binom.test, and wilcox.test.  Defaults to "two.sided".  Can be a vector.
#' @param altvalue vector of alternative hypothesis values.  Defaults to null value for calculation of actual test level.
#' @param flip logical indicator to allow asymetric distributions to be flipped.
#'
#' @return A five-dimensional array of test levels (if altevalue equals hypoth) or power otherwise.  Dimensions are test, distribution, alternative, alternative value, and first dimension sample size.
#'
#' @examples
#' fun.comparepower(10,15,nsamp=100)
#'
#' @details Tests are given by fun.onesampletests or fun.twosampletests, as appropriate.
#' @export
fun.comparepower<-function(samp1sz=10,samp2sz=NA,dist="rnorm",ndim=1,distb=NULL,level=.05, hypoth=0,nsamp=1000,alternative="two.sided",altvalue=NULL,flip=F){
#  if(!is.list(dist)) dist<-list(dist)
   if(is.null(distb)) distb<-dist
#  if(!is.list(distb)) distb<-list(distb)
   if(is.null(altvalue)) altvalue<-0
   count<-NULL
   if((length(hypoth)!=1)&&(length(hypoth)!=length(dist))){
       cat("Error: length of hypoth wrong\n")
   }else{
       if(length(hypoth)==1) hypoth<-rep(hypoth,length(dist))
   }
   for(m in seq(length(samp1sz))){
      for(i in seq(length(altvalue))){
         for(l in seq(length(alternative))){
            for(k in seq(length(dist))){
               thisdist<-eval(parse(text=dist[k]))
               for(j in 1:nsamp){
                  x<-thisdist(samp1sz[m]*ndim)+hypoth[k]
                  if(ndim>1){
                     x<-array(x,c(samp1sz[m],ndim))
                     if(flip){
                        flipme<-abs(x[,2])<1
                        x[flipme,2]<-x[flipme,1]
                        x[!flipme,2]<--x[!flipme,1]
                     }
                  }
                  if(!is.na(samp2sz)){
                     secdist<-eval(parse(text=distb[k]))
                     y<-secdist(samp2sz*ndim)+altvalue[i]
                     if(ndim>1) x<-array(x,c(samp2sz,ndim))
                     testout<-fun.twosampletests(x,y)
#                    browser()
#                    cat(testout);cat("\n"); cat(x); cat("\n") ; cat(y) ; cat("\nxxxx\n")
                  }else{
                     testout<-fun.onesampletests(x,alternative[l])
                  }
                  if(is.null(count)){
                     count<-array(0,c(length(testout),length(dist),length(alternative),length(altvalue),length(samp1sz)))
                     dnn<-if(flip) paste(dist,"flip") else dist
                     dimnames(count)<-list(names(testout),dist,alternative,as.character(round(altvalue,3)),as.character(samp1sz))
                  }
                  count[,k,l,i,m]<-count[,k,l,i,m]+(testout<level)
               }
            }
         }
      }
   }
   return(count/nsamp)
}
