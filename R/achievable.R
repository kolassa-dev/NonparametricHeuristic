#' List the achievable p-values for various sample sizes.
#' @export
#' @examples
#' fun.achievable()
fun.achievable<-function(){
   u<-array(NA,c(41,4))
   for(n in seq(dim(u)[1])) {
      v<-rep(NA,(n+1)/2)
      for(i in seq(length(v))) v[i]<-binom.test(i-1,n)$p.value
      u[n,2]<-sum(v<=.05)-1
      if(u[n,2]>-1) u[n,4]<-max(v[v<=.05])
      w<-v
      if(n>2){
         for(i in seq(length(w))) w[i]<-z.test(c(rep(0,n-i),
            .5,rep(1,i-1)),mu=.5,sigma.x=.5)$p.value
#        u[n,3]<-max(w[w<=.05])
         u[n,3]<-sum(w<=.05)-1
      }
      u[n,1]<-n
   }
   dimnames(u)<-list(NULL,c("n","Exact Critical Value",
      "Asymptotic Critical Value","Exact Size"))
#  v<-cbind(u[6:17,],u[18:29,],u[30:41,])
   v<-cbind(u[6:23,],u[24:41,])
   return(v)
}
