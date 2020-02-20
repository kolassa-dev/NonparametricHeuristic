#' Draw sample paths associated with the Renyi test based on the log rank statistic.
#'
#' @param time event/censoring times.  If null, hypothetical data set will be censored and plotted.
#' @param delta censoring indicator
#' @param group variable defining group.  
#' @param nsamp If time is null, number of random observations.
#' @param p if data are random, gives multiplier for alternative.  p=0 for null.
#' @return A list with components m for the statistic, x for ordered times, and y for statistic matrix.
#'
#' @examples
#' reyniplot(nsamp=20)
#'
#' @importFrom stats rbinom rexp
#' @export
reyniplot<-function(time=NULL,delta=NULL,group=NULL,nsamp=100,p=1){
   if(is.null(time)){
      group<-1+rbinom(nsamp,1,.5); time<-rexp(nsamp)*group^p; delta<-rbinom(nsamp,1,.5)
   }
   group<-as.numeric(as.factor(group))
   ut<-unique(sort(time))
   ug<-unique(group)
   ng<-length(ug)
   ev<-ce<-ar<-array(NA,c(length(ut),ng))
   ar[1,]<-table(group)
   for(jj in seq(length(ut))){
      thistime<-time==ut[jj]
      for(ii in ug){
         ev[jj,ii]<-sum(delta[thistime]*(group[thistime]==ii))
         ce[jj,ii]<-sum((1-delta[thistime])*(group[thistime]==ii))
      }
      if(jj<length(ut)) ar[jj+1,]<-ar[jj,]-ev[jj,]-ce[jj,]
   }
   evt<-apply(ev,1,sum); art<-apply(ar,1,sum)
   pv<-sum(ar[,1]*ar[,2]*evt*(art-evt)/(art^2*(art-1+(art==1))))
   temp<-outer(evt/art,c(1,1))
   stv<-apply(ev-outer(evt/art,c(1,1))*ar,2,cumsum)/sqrt(pv)
   m<-max(abs(stv[,1]))
   i<-match(m,abs(stv[,1]))/length(ut)
   plot(c(0,1),range(stv),type="n")
   for(j in 1:2) lines((0:length(ut))/length(ut),c(0,stv[,j]),lty=j)
   abline(v=i,lty=3)
   legend(0,max(stv),lty=1:2,legend=paste("Group",1:2))
   return(invisible(list(max=m,x=ut,y=stv)))
}
#temp<-survdiff(Surv(time,delta)~group)
#reyniplot(time,delta,group)
