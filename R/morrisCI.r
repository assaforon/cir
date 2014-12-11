### Morris and Morris-style CI: functions and accessories

Gupper<-function(theta,y,n,k)
{
if(k==length(y)) return (pbinom(q=y[k],size=n[k],prob=theta))
pbinom(q=y[k]-1,size=n[k],prob=theta)+dbinom(x=y[k],size=n[k],prob=theta)*Gupper(theta=theta,y=y,n=n,k=k+1)
}	

Glower<-function(theta,y,n,k)
{
if(k==1) return (pbinom(q=y[k]-1,size=n[k],prob=theta,lower.tail=FALSE))
pbinom(q=y[k],size=n[k],prob=theta,lower.tail=FALSE)+dbinom(x=y[k],size=n[k],prob=theta)*Glower(theta=theta,y=y,n=n,k=k-1)
}	


##' Confidence intervals for ordered Binomial, based on Morris 1988

 
morrisUCL<-function(y,n,halfa=0.05)
{
m=length(y)
if(length(n)!=m) stop("Mismatched lengths in Morris.\n")
# weird prep...
uout=rep(1,m)
a=m
### At uppermost doses as long as phat=1, no need for algorithms
while(y[a]==n[a] && a>=1) a=a-1
if(a<1) return(uout)

for(b in a:1)
{
	uout[b]=uniroot(function(theta,h,d,alpha,...) h(theta=theta,...)-alpha,interval=c(0,1),
		alpha=halfa,k=b,h=Gupper,n=n,y=y)$root
}
return(uout)
}
#####

morrisLCL<-function(y,n,halfa=0.05)
{
m=length(y)
if(length(n)!=m) stop("Mismatched lengths in Morris.\n")
# weird prep...
uout=rep(0,m)
a=1
### At lowermost doses as long as phat=0, no need for algorithms
while(y[a]==0 && a<=m) a=a+1
if(a>m) return(uout)

for(b in a:m)
{
	uout[b]=uniroot(function(theta,h,d,alpha,...) h(theta=theta,...)-alpha,interval=c(0,1),
		alpha=halfa,k=b,h=Glower,n=n,y=y)$root
}
return(uout)
}
#####

morrisCI<-function(y,n,phat=y/n,conf=0.9,narrower=TRUE,alternate=wilsonCI,...)
{
tailp=(1-conf)/2
lcl=morrisLCL(y=y,n=n,halfa=tailp)
ucl=morrisUCL(y=y,n=n,halfa=tailp)

if(narrower)
{
	altout=alternate(phat=phat,n=n,conf=conf,...)
	lcl=pmax(lcl,altout[,1])
	ucl=pmin(ucl,altout[,2])
}

return(cbind(lcl,ucl))
}

	
#uout[[a]]=wilsonCI(phat=x[a]/n[a],n=n[a],conf=1-2*halfa)[2]
#if(a<=1) return(uout)

	