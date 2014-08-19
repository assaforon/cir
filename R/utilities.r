
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

#### Some standard CI calculation functions

wilsonCI<-function(phat,n,conf=0.9)
{
zalpha=qnorm(1-(1-conf)/2)
wid=zalpha*sqrt((phat*(1-phat)+zalpha^2/(4*n))
/n)
return(cbind(pmax(0,(phat+zalpha^2/(2*n)-wid)/(1+zalpha^2/n)),
	pmin(1,(phat+zalpha^2/(2*n)+wid)/(1+zalpha^2/n))))
}

agcouCI<-function(phat,n,conf=0.9)
{
zalpha=qnorm(1-(1-conf)/2)
ptilde=(phat+zalpha^2/(2*n))/(1+zalpha^2/n)
ntilde=n+zalpha^2
wid=zalpha*sqrt(ptilde*(1-ptilde)/ntilde)
return(cbind(pmax(0,ptilde-wid),
	pmin(1,ptilde+wid)))
}

jeffCI<-function(phat,n,conf=0.9,w1=0.5,w2=w1)
{
x=n*phat
clow=ifelse(phat==0,0,qbeta((1-conf)/2,x+w1,n-x+w2))
chigh=ifelse(phat==1,1,qbeta(1-(1-conf)/2,x+w1,n-x+w2))

return(cbind(clow,chigh))
}
