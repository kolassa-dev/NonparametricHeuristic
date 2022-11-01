#' Calculate a vector of two-sample tests of the null hypothesis of equal center
#'
#' Calculates t-test, exact Wilcoxon, approximate Wilcoxon, normal scores, Savage scores, and Mood tests.
#'
#' @param x the data vector for the first sample.
#' @param y the data vector for the second sample.
#'
#' @return A vector of two-sided p-values
#'
#' @examples
#' fun.twosampletests(rnorm(10),rnorm(15))
#'
#' @details The exact Wilcoxon is performed only for small samples.
#' @importFrom stats t.test wilcox.test 
#' @export
fun.twosampletests<-function(x,y){
   ef<-(length(x)<20)&(length(y)<20)
   pv<-c(t.test(x,y)$p.value,
      wilcox.test(x,y,exact=ef,correct=F)$p.value,
      wilcox.test(x,y,exact=F,correct=F)$p.value,
      ranktests(x,y,"ns")$p.value,ranktests(x,y,"ss")$p.value,mood.median.test(x,y)$p.value
      )
   names(pv)<-c("T-test","Exact Wilcoxon","Approximate Wilcoxon", "Normal Scores", "Savage Scores","Mood")
   return(pv)
}
