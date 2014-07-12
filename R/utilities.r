
# one-sixth the skew of the Binomial score
scoreSkew6<-function(theta,n) (1-2*theta)/(6*sqrt(n*theta*(1-theta)))

### Non-vectorized for now
bcaFwd<-function(bootdat,thetaHat,n,probs,full=FALSE)
{
# BCa parameters
z0=qnorm(mean(bootdat<=thetaHat))
accel=scoreSkew6(thetaHat,n)
cat(c(z0,accel),'\n')

zalpha=qnorm(probs)
zbca=z0+(z0+zalpha)/(1-accel*(z0+zalpha))
cat(zbca,'\n')
endpoints=quantile(bootdat,prob=pnorm(zbca),type=6)
if(!full) return(endpoints)
return(list(enddpoints=endpoints,a=accel,z0=z0))
}

#shrinkIt<-function(
#if(is.null(prior)) prior=doseResponse(y=(1:m)/(m+1),x=xout,wt=rep(1/m,m))
