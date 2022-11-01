#' Calculate a vector of one-sample tests of the null hypothesis of median 0
#'
#' Calculates z, sign, and wilcoxon rank sum tests of the null hypothesis of center 0
#'
#' @param x the data vector
#' @param alternative alternative, of the form aceptable by t.test, z.test, binom.test, and wilcox.test
#'
#' @return A vector of p-values
#'
#' @examples
#' fun.onesampletests(rnorm(10))
#'
#' @importFrom stats t.test binom.test wilcox.test 
#' @importFrom BSDA z.test
#' @importFrom ICSNP rank.ctest HotellingsT2
#' @export
fun.onesampletests<-function(x,alternative="two.sided"){
   if(is.vector(x)){
      pv<-c(t.test(x ,alternative=alternative)$p.value, 
         z.test((x>0)-sign(mean(x>0)-.5)/(2*length(x)),mu=.5,sigma.x=.5,
            alternative=alternative)$p.value,
         binom.test(sum(x>0),length(x),alternative=alternative)$p.value,
         wilcox.test(x,alternative=alternative)$p.value)
      names(pv)<-c("T","Sign","Exact Sign","Signed Rank")
   }else{
      test1<-try(rank.ctest(x),silent=TRUE)
      if(inherits(test1,"try-error")) test1<-list(p.value=1)
      test2<-try(rank.ctest(x,scores="sign"),silent=TRUE)
      if(inherits(test2,"try-error")) test2<-list(p.value=1)
#     browser()
      pv<-c(HotellingsT2(x)$p.value,test1$p.value,test2$p.value)
      names(pv)<-c("T2","Sign test","Sign rank test")
   }
   return(pv)
}
