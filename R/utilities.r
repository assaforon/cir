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

#' Linearly interpolate/extrapolate from a single segment.

extrapol<-function(point1,point2,xout)
{
slope=(point2[2]-point1[2])/(point2[1]-point1[1])
return(point2[2]+(xout-point2[1])*slope)
}

