
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


#### bootstrap-t for Binomial, with variance shrinkage
bootTFwd<-function(bootdat,thetaHat,n,probs,full=FALSE,prior,nprior)
{
theta0Shrink=(thetaHat*n+prior*nprior)/(n+nprior)
thetaBootShrink=(bootdat*n+prior*nprior)/(n+nprior)
#cat(theta0Shrink,quantile(thetaBootShrink),'\n')
pivots=(thetaHat-bootdat)/sqrt(thetaBootShrink*(1-thetaBootShrink)/n)

endpoints=thetaHat+sqrt(theta0Shrink*(1-theta0Shrink)/n)*quantile(pivots,prob=probs,type=6)
if(!full) return(endpoints)
return(list(endpoints=endpoints,pivots=pivots))
}

