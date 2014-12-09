#### Some standard CI calculation functions

wilsonCI<-function(phat,n,conf=0.9)
{
zalpha=qnorm(1-(1-conf)/2)
wid=zalpha*sqrt((phat*(1-phat)+zalpha^2/(4*n))/n)
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

########### Intervals from Morris 1988

morrisUCL<-function(x,n,halfa=0.05,full=FALSE)
{
m=length(x)
if(length(n)!=m) stop("Mismatched lengths in Morris.\n")
# weird prep...
uout=rep(1,m)
G=list()
for (a in 1:m) G[[a]]<-function(x) 1

a=m
while(x[a]==n[a] && a>=1) a=a-1
uout[[a]]=wilsonCI(phat=x[a]/n[a],n=n[a],conf=1-2*halfa)[2]
if(a<=1) return(uout)
G[[a]]<-function(theta,g,y,n,k) pbinom(q=y[k],size=n[k],prob=theta)

for(b in (a-1):1)
{
	G[[b]]=function(theta,g,y,n,k) 
		pbinom(q=y[k]-1,size=n[k],prob=theta)+dbinom(x=y[k],size=n[k],prob=theta)*g[[k+1]](theta=theta,g=g,y=y,n=n,k=k+1)
	uout[b]=uniroot(function(theta,h,d,alpha,...) h[[d]](theta=theta,...)-alpha,interval=c(0,1),
		alpha=halfa,d=b,k=b,g=G,h=G,n=n,y=x)$root
}
if(full) return(list(u=uout,G=G))
return(uout)
}
	

	
########### Interpolation, Extrapolation, Parapolation....

#' Linearly interpolate/extrapolate from a single segment.

extrapol<-function(point1,point2,xout)
{
slope=(point2[2]-point1[2])/(point2[1]-point1[1])
return(point2[2]+(xout-point2[1])*slope)
}

#' Locally monotone quadratic interpolation between points on a plane

parapolate<-function(x,y,xout,upward,full=FALSE)
{
## validations
if(any(diff(x)<=0)) stop("x must be monotone strictly increasing.\n")
m=length(x)
if(m<2) stop ("Need at least 2 points.\n")
if(min(xout)<x[1] || max(xout)>x[m]) stop("This function does not extrapolate.\n")
if(length(y)!=m) stop("mismatched x,y lengths.\n")

yout=y[match(xout,x)]
if(!any(is.na(yout))) return(yout)  ### no interpolation to do

## These, to help make vectorized and interpretable formulae
x1=x[-m]
y1=y[-m]
x2=x[-1]
y2=y[-1]
denom=diff(x)^2
# A bit of boolean finesse is needed:
if(length(upward)==1) upward=rep(upward,m-1)
oneflat=xor(upward,diff(y)<0)
# Explanation: 'oneflat' determines whether parabola is 'flat' at point 1 or 2 of each segment

## parabola coefficients
a=(y1-y2)/denom
a[oneflat]=-a[oneflat]

b=ifelse(oneflat,-2*x1*a,-2*x2*a)
cee=y1-a*x1^2-b*x1  # per parabola definitions :)

intchoose=findInterval(xout,x)
candy=a[intchoose]*xout^2+b[intchoose]*xout+cee[intchoose]
yout[!(xout %in% x)]=candy[!(xout %in% x)]
if(!full) return(yout)
return(list(a=a,b=b,c=cee,outdat=data.frame(x=xout,y=yout)))
}

#k=length(xout)


