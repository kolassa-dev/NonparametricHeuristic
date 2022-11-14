#' Plot the results of a per-value structure of confidence intervals.  Vestigial.
#'
#' @param tauv Vector of ordinates for the empirical cdf
#' @param cistr matrix with two rows, and as many colums as elements in tauv.  First row is lower end of CI, and second row is upper end of CI.
#' @export
#' @importFrom graphics segments
#'
fun.plotexactci<-function(tauv,cistr){
   for(j in 1:2){
      uvec<-unique(cistr[j,])
      for(u in uvec){
         rr<-range(tauv[cistr[j,]==u])
#        cat(j,u,rr,"\n")
         segments(u,rr[1],u,rr[2])
      }
   }
}
#' Intermdiate function for calculating Mood median test.  Odd case is calculated recursively. The even-sample size question is reduced either to Fisher's exact test, or to the Gaussian approximation.
#' @param u One corner of 2x2 contingency table
#' @param rt Row totals
#' @param ct Column totals
#' @param exact Logical determining whether Fisher's test or the approximate counterpart is used.
#' @param alternative alternative, defaults to "two.sided".
#' @return A list with one component, the p-value.
#' @importFrom stats fisher.test pnorm
mmti<-function(u,rt,ct,exact=FALSE,alternative="two.sided"){
   if(rt[2]==0){
      tt<-cbind(c(u,ct[1]-u),
         c(rt[1]-u,rt[3]-ct[1]+u))
      if(exact){
         p.value<-fisher.test(tt,alternative=alternative)$p.value
      }else{
         nn<-sum(rt)
         v<-prod(rt[-2])*prod(ct)/(nn^2*(nn-1))
         ucr<-u-.5*sign(u)
         z<-(ucr/2)/sqrt(v)
         if(alternative=="two.sided") p.value<-2*pnorm(-abs(z))
         if(alternative!="two.sided") p.value<-pnorm(-z)
      }
   }else{
      rt[2]<-0
      p.value<-
         mmti(u,rt,ct-c(1,0),exact=exact,alternative=alternative)$p.value*ct[1]/sum(ct)+
         mmti(u,rt,ct-c(0,1),exact=exact,alternative=alternative)$p.value*ct[2]/sum(ct)
   }
   return(list(p.value=p.value))
}
#' Perform Mood's median test
#' @param x First sample
#' @param y Second sample
#' @param exact logical, determing whether the exact version is done.
#' @param alternative Determines whether two-sided test, or which one-sided test, is done.
#' 
mood.median.test<-function(x,y,exact=F,alternative="two.sided"){
  mm<-median(c(x,y))
  a21<-sum(x>mm); a22<-sum(x<mm)
  a11<-sum(y>mm); a12<-sum(y<mm)
  a13<-sum(x==mm); a23<-sum(y==mm)
  return(mmti(a21-a22,c(a21+a11,a13+a23,a22+a12),
     c(a11+a12+a13,a21+a22+a23),exact=exact,
     alternative=alternative))
}
#' Compare coverage of bootstrap confidence intervals from regression.
#'
#' @param dists vector of names of random number generators for sample.  Defaults to normal.
#' @param regparam function to pass to the bootstrapping tool
#' @param sampsize Number of observations per sample; defaults to 10
#' @param mcsamp Number of Monte Carlo samples.
#' @param B Bootstrap size
#'
#' @return A matrix of coverage levels.  Dimensions are distribution and confidence interval technique.
#'
#' @examples
#' regparam<-function(data,indices,fittedv,residuals){
#'    y<-fittedv+residuals[indices]
#'    return(summary(lm(y~data[,1]))$coefficients[2,1:2]^(1:2))}
#' fun.testregboot(c("rnorm","rcauchy","rexp"),regparam,mcsamp=20,B=99)
#'
#' @export
#' @importFrom stats lm qt resid
fun.testregboot<-function(dists,regparam,sampsize=10, mcsamp=10000, B=9999){ 
   x<-1:sampsize
   out<-array(0,c(5,length(dists)))
   dimnames(out)<-list(
      c("normal","basic","student","percent","bca"),
      dists)
   for(dist in dists){
      for(j in 1:mcsamp){
         y<-x+eval(parse(text=dist))(sampsize)
         residuals<-resid(lm(y~x))
         fittedv<-y-residuals
         ciout<-boot.ci(boot(cbind(x,y),regparam,B,
            fittedv=fittedv,residuals=residuals))
         for(n in names(out)){
            r<-rev(ciout[[n]])
            out[n,dist]<-out[n,dist]+((r[1]>1)&(r[2]<1))
         }
      }
   }
   return(out/mcsamp)
}
#' Collection of location measures
#' @param y data vector
#' @returns mean median, and trimmed mean
#' @export
#' @importFrom stats median
location<-function(y){
    out<-c(mean(y),median(y),mean(y,trim=.1))
    names(out)<-c("mean","median","trim mean"); return(out)}
#' Collection of scale measures
#' @param y data vector
#' @returns standard deviation, interquartile range, and mean absolute deviation.
#' @export
#' @importFrom stats sd IQR mad
scale<-function(y){out<-c(sd(y),IQR(y),mad(y))
    names(out)<-c("sd","IQR","MAD"); return(out)}
sensitivity.plot<-function(y,sub,stats){
   ra<-range(y)
   xr<-mean(ra)+c(-1,1)*diff(ra)
   outlier<-xr[1]+(0:100)*diff(xr)/100
   base<-stats(y)
   sens<-array(NA,c(3,length(outlier)))
   dimnames(sens)[[1]]<-names(base)
   for(j in seq(length(outlier))){
      sens[,j]<-stats(c(y,outlier[j]))-base
   }
   plot(xr,range(sens),type="n",main="Sensitivity",sub=sub)
   inds<-seq(length(base))
   for(i in inds) lines(outlier,sens[i,],lty=i,col=i)
   legend(min(ra),max(sens),lty=inds,col=inds,
      legend=names(base))
}
#' Diagram for the construction of Feiller confidence intervals
#'
#' @param xv vector of X means
#' @param yv vector of Y means
#' @param z critical value for Z test.
#' @examples
#' fun.feiller(c(1,3,3),c(3,1.5,.3),1.96)
#' @export
#' @importFrom graphics par lines points abline legend mtext
fun.feiller<-function(xv,yv,z){
   par(mfrow=c(2,1))
   par(mar=c(0, 4, 4, 1) + 0.1)
   rho<-(-100:100)/40
   ff<-function(x,y,rho) return((x-rho*y)/sqrt(1+rho^2))
   gg<-function(x,y,rho,z) return((x-rho*y)^2-(1+rho^2)*z^2)
   gout<-fout<-array(NA,c(length(xv),length(rho)))
   epa<-array(NA,c(2,length(xv)))
   for(i in 1:length(xv)){
      epa[,i]<-(xv[i]/yv[i]+c(-1,1)*(z/yv[i])*sqrt(
         (xv[i]/yv[i])^2+1-(z/yv[i])^2))/(1-(z/yv[i])^2)
      fout[i,]<-ff(xv[i],yv[i],rho)
      gout[i,]<-gg(xv[i],yv[i],rho,z)
   }
   plot(range(rho),range(fout),type="n",xlab="",#xlab="rho",
#     sub=paste("n=1, variance=1, z=",1.96),
      ylab="Z Statistic"#,main="Ratio of Normal Means CI"
      )
   for(i in 1:length(xv)){
      lines(rho,fout[i,],lty=i)
      points(epa[,i],ff(xv[i],yv[i],epa[,i]))
   }
   abline(h=z^2)
   legend(mean(rho),z,lty=seq(length(xv)),legend=paste("x=",xv,"y=",yv))
   par(mar=c(4, 4, 0, 1) + 0.1)
   plot(range(rho),range(gout),type="n",xlab="",#xlab="rho",
#     sub=paste("n=1, variance=1, z=",1.96),
      ylab="Quadratic Form"#,main="Ratio of Normal Means CI"
      )
   abline(h=0)
   for(i in 1:length(xv)){
      lines(rho,gout[i,],lty=i)
      points(epa[,i],gg(xv[i],yv[i],epa[,i],z))
   }
   legend(mean(rho),0,lty=seq(length(xv)),legend=paste("x=",xv,"y=",yv))
   mtext("Ratio of Normal Means CI",outer=TRUE,side=3,line=-3)
   mtext(paste("n=1, variance=1, z=",1.96),outer=TRUE,side=1,line=-1)
   mtext("rho",outer=TRUE,side=1,line=-2)
   return(epa)
}
#' Siegel-Tukey rank function.
#' @param x the data vector
#' @return a vector with the same length as the input, with Siegel-Tukey Ranks.
#' @details
#' The complete test in one place is at https://raw.github.com/talgalili/R-code-snippets/master/siegel.tukey.r
#' @examples
#' siegel.tukey.ranks(c(rnorm(10),10*rnorm(10)))
#' @export
siegel.tukey.ranks<-function(x) {
   n<-length(x)
   sss<-rank(x)
   top<-FALSE
   count<-0
   rr<-rep(NA,n)
   lowest<-0
   highest<-n+1
   for(i in seq(n)){
      top<-(i-0)%%4>1
      if(top){
         highest<-highest-1
         rr[highest]<-i
      }else{
         lowest<-lowest+1
         rr[lowest]<-i
      }
   }
   out<-rr[rank(x,ties.method="first")]
   for(xx in unique(x)) out[x==xx]<-mean(out[x==xx])
   return(out)
}
#' Calculate the general multi-group score test.
#' @param x response vector
#' @param g group indicator.
#' @param s scores to be averaged over ties.
#' @return an object of class htest, with the chi-square approximate p-valu.e
#' @export
#' @importFrom stats approx aggregate pchisq
genmultscore<-function(x,g,s=NULL){
   if(is.null(s)) s<-seq(length(x))
# Average over tied values
   news<-sort(approx(seq(length(x)),s,rank(x))$y)
   abar<-mean(news)
   ahat<-mean(news^2)
   sigsq<-ahat-abar^2
   reps<-table(g)
   numtot<-sum(reps)
#  vv<-(numtot*diag(reps)-outer(reps,reps))/(numtot-1)
   rs<-aggregate(news,list(g=g),FUN=sum)$x
   H<-sum((rs-reps*abar)^2/reps)*(numtot-1)/(numtot*sigsq)
   pp<-length(reps)-1
   names(pp)<-"df"
   out<-list(parameter=pp,data.name=NA, 
      statistic=H,method="General Rank Score Multi-Group Test",
      p.value=1-pchisq(H,length(reps)-1))
   class(out)<-"htest"
   return(out)
}
#' Diagram the critical region for a three-group score test, as a function of the first two dimensions.
#' @param nrep Number in each group.
#' @param scores Scores.
#' @param docontour logical value controlling drawing contour.
#' @examples
#' showmultigroupscoretest(c(3,3,4))
#' @export
#' @importFrom MultNonParam nextp
#' @importFrom graphics contour
#' @importFrom stats qchisq
showmultigroupscoretest<-function(nrep,scores=NULL,docontour=TRUE){
   if(length(nrep)==3){
      g<-rep(c(1,2,3),nrep)
      oldg<-0
      out1<-NULL
      nn<-sum(nrep)
      if(is.null(scores)){
         scores<-(1:nn)
      }else{
         if(length(scores)==nn){
            docontour<-FALSE
         }else{
            message("Scores are the wrong length.  Resetting to ranks.")
            scores<-(1:nn)
         }
      }
      while(!all(oldg==g)){
         oldg<-g
         g<-nextp(g)
         me<-c(sum(scores[g==1]),sum(scores[g==2]))
         if(!all(oldg==g)){
            out1<-if(is.null(out1)) array(me,c(1,2)) else rbind(out1,me)
         }
      }
#     browser()
      out<-unique(as.data.frame(out1))
      out<-cbind(out,nn*(nn+1)/2-apply(out,1,sum))
      if(docontour){
         ers<-(nn+1)*nrep/2
#        kw<-apply(apply(apply(out,2,"-",ers)^2,2,"/",nrep),1,sum)*12/(nn*(nn+1))
         kwx<-min(out[,1])+-1:(diff(range(out[,1]))+1)
         kwy<-min(out[,2])+-1:(diff(range(out[,2]))+1)
         kwout<-array(NA,c(length(kwx),length(kwy)))
         for(ii in seq(length(kwx))) for(jj in seq(length(kwy))){
            xx<-kwx[ii]
            yy<-kwy[jj]
            zz<-nn*(nn+1)/2-xx-yy
            kwout[ii,jj]<-12*sum((c(xx,yy,zz)-ers)^2/nrep)/(nn*(nn+1))
         }
         contour(kwx,kwy,kwout, levels=qchisq(.95,2),
            xlab="Group 1 rank sum",ylab="Group 2 rank sum",
            main="Asymptotic Critical Region for Kruskall Wallis Test, level 0.05",
            sub=paste("Group sizes",paste(nrep,collapse=" "),sep=" "))
      }else{
         plot(range(out[,1]),range(out[,2]),type="n",
            xlab="Group 1 rank sum",ylab="Group 2 rank sum",
            main="Asymptotic Critical Region for Kruskall Wallis Test, level 0.05",
            sub=paste("Group sizes",paste(nrep,collapse=" "),sep=" "))
      }
      points(out[,1],out[,2])
   }else{
      message("You can only do this in three dimensions")
   }
}
#' Diagram the construction of the one-sample Hodges-Lehmann estimator
#' @param zz Data vector.  Defaults to a random Cauchy sample of an appropriate length.
#' @param nn Length of data vector.  Ignored if no real vector supplied, and defaults to 6.
#' @param alpha 1-confidence level.  Defaults to 0.05.
#' @return a vector with two components giving the confidence interval.
#' @examples
#' hodgeslehmannexample()
#' @export
#' @importFrom stats qsignrank
hodgeslehmannexample<-function(zz=NULL,nn=6,alpha=0.05){
   if(is.null(zz)){
      subt<-paste("Random Cauchy Example of size",nn)
      zz<-rcauchy(nn) 
   }else{
      subt<-paste("Data example of size",length(zz))
      nn<-length(zz)
   }
   wa<-outer(zz,zz,"+")/2
   wa<-sort(wa[!lower.tri(wa)])
   eps<-0.0001
   plot(c(min(wa)-1,max(wa)+1),c(0,nn*(nn+1)/2),type="n",sub=subt,
      xlab="Potential point of symmetry",ylab="Signed Rank Statistic",
      main="Construction of Median Estimator")
   for(ii in seq(nn*(nn+1)/2-1)) 
      segments(wa[ii],nn*(nn+1)/2-ii,wa[ii+1],nn*(nn+1)/2-ii)
   segments(wa[1]-1,nn*(nn+1)/2,wa[1],nn*(nn+1)/2)
   segments(wa[nn*(nn+1)/2]+1,0,wa[nn*(nn+1)/2],0)
   tl<-qsignrank(alpha/2,nn)
   tu<-nn*(nn+1)/2+1-tl
   axis(4,at=c(tl,tu-1),labels=c("tl","tu-1"))
   abline(h=tl,lty=2)
   abline(h=tu-1,lty=2)
   abline(v=wa[tl],lty=3)
   abline(v=wa[tu],lty=3)
   return(c(tl,tu))
}
#' Inverts the relationship giving standardized resituals from raw residuals as given, for example, in lmwork from library MASS .
#' @param stdres Standardized residuals
#' @param modelfit Fit of the linear model
#' @param hat Vector of hat values
#' @param stddev Sample standard deviation of residuals.
#' @return raw residuals
#' @importFrom stats lm.influence
stdres2resid<-function(stdres,modelfit=NULL,hat=NULL,stddev=NULL){
   if(is.null(hat)) hat<-lm.influence(modelfit,do.coef=FALSE)$hat
   if(is.null(stddev)) stddev<-sqrt(sum(modelfit$residuals^2)/modelfit$df.residuals)
   return(stdres*sqrt(1-hat)*stddev)
}
#' Test bias of the jackknife approach.
#' @param fcns functions to jackknife.  Must be given as a list.
#' @param dists Random number generator, or list of random number generators.
#' @param nsamp number of data sets to sample
#' @param npersamp number of observations per sample
#' @param nbig number of Monte Carlo examples to get target value.
#' @param others list of additional arguments for entries in components of fcns.
#' @return Expectation, jackknife bias, and true value (via large MC)
#' @examples
#' testjack(list(mean,median,mean),dists=list(rexp,rnorm),
#'    others=list(NULL,NULL,list(trim=0.25)),nsamp=5)
#' @importFrom bootstrap jackknife
#' @importFrom stats rexp
#' @export
testjack<-function(fcns,dists=rexp,nsamp=1000,npersamp=11,nbig=NA,others=NULL){
   if(!is.list(dists)) dists<-list(dists)
   if(is.na(nbig)) nbig<-npersamp*1000
   if(is.null(others)) others<-vector("list",length(fcns))
   out<-array(NA,c(length(dists),length(fcns),nsamp,2))
   meanout<-array(NA,c(length(dists),length(fcns),3))
   dimnames(meanout)<-list(
      as.character(substitute(dists))[-1],
      as.character(substitute(fcns))[-1],
      c("Expectation of Statistic","Expectation of Bias Estimate","Parameter"))
   for(ii in seq(length(dists))){
      for(k in seq(length(fcns))){
         fcn<-fcns[[k]]
         for(j in seq(dim(out)[3])){
            x<-dists[[ii]](npersamp)
            mylist<-list(x=x)
            if(!is.null(others[[k]])) mylist<-c(mylist,others[[k]])
            out[ii,k,j,1]<-do.call(fcn,mylist)
            mylist<-list(x=x,theta=fcn)
            if(!is.null(others[[k]])) mylist<-c(mylist,others[[k]])
#           browser()
            out[ii,k,j,2]<-do.call(jackknife,mylist)$jack.bias
         }
#        meanout[ii,k,1:2]<-apply(out[ii,k,,],3,mean)
         mylist<-list(x=dists[[ii]](nbig))
         if(!is.null(others[[k]])) mylist<-c(mylist,others[[k]])
         meanout[ii,k,3]<-do.call(fcn,mylist)
      }
   }
   meanout[,,1:2]<-apply(out,c(1,2,4),mean)
   return(meanout)
}
