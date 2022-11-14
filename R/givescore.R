#' Compare normal, van der Waerden, or Savage scores.
#'
#' @param data Vector of data to score
#' @param sv Type of scores, if character "ns", "vw", or "ss", or vector of scores if numeric.  Defaults to rank scores.
#'
#' @return A vector of scores.
#'
#' @examples
#' fun.givescore(seq(15))
#'
#' @details Tests are given by fun.onesampletests or fun.twosampletests, as appropriate.
#' @export
#' @importFrom stats qqnorm
fun.givescore<-function(data,sv=NULL){
   totsamp<-length(data)
   score<-NULL
   if(is.character(sv)){
      if(sv=="ns") score<-qqnorm(seq(totsamp),plot.it=F)$x
      if(sv=="vw") score<-qnorm(seq(totsamp)/(totsamp+1))
      if(sv=="ss") score<-cumsum(1/rev(seq(totsamp)))
      if(is.null(score)) cat("Error: character score does not match allowable options\n")
   }
   if(is.numeric(sv)) score<-sv
   if(is.null(score)) score<-seq(totsamp)
#  browser()
   score<-(score-mean(score))[rank(data,ties.method="first")]
   for(xx in unique(data)) score[data==xx]=mean(score[data==xx])
   return(score)
}
