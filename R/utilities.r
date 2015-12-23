#### Some standard CI calculation functions
#' Standard unordered-Binomial confidence interval utilitlies.
#' 
#' Standard small-sample Binomial confidence interval utilities, using the methods of Wilson, Agresti-Coull and Jeffrys.
#' 
#' These functions implement the basic (uncorrected) three intervals which are seen by the consensus of literature as the "safest" off-the-shelf formulae. None of them account for ordering or monotonicity; therefore the \code{cir} package default is \code{\link{morrisCI}} which does account for that, with the 3 unordered formulae used for optional narrowing of the interval at individual points.

#' @aliases agcouCI jeffCI
#' @seealso \code{\link{isotInterval}} for more details about how forward CIs are calculated, \code{\link{quickInverse}} for inverse (dose-finding) intervals.
#' @export

#' @return A two-column matrix with the same number of rows as \code{length(phat)}, containing the calculated lower and upper bounds, respectively.
#' 
#' @param phat numeric vector, point estimates for which an interval is sought
#' @param n integer vector of same length, of pointwise sample sizes
#' @param conf numeric in (0,1), the confidence level
#' @param w1,w2 numeric, weights used in \code{jeffCI} only
#' @param ... pass-through for compatibility with a variety of calling functions
#' 
wilsonCI<-function(phat,n,conf=0.9,...)
{
zalpha=qnorm(1-(1-conf)/2)
wid=zalpha*sqrt((phat*(1-phat)+zalpha^2/(4*n))/n)
return(cbind(pmax(0,(phat+zalpha^2/(2*n)-wid)/(1+zalpha^2/n)),
	pmin(1,(phat+zalpha^2/(2*n)+wid)/(1+zalpha^2/n))))
}

#' @rdname wilsonCI
#' @export
agcouCI<-function(phat,n,conf=0.9,...)
{
zalpha=qnorm(1-(1-conf)/2)
ptilde=(phat+zalpha^2/(2*n))/(1+zalpha^2/n)
ntilde=n+zalpha^2
wid=zalpha*sqrt(ptilde*(1-ptilde)/ntilde)
return(cbind(pmax(0,ptilde-wid),
	pmin(1,ptilde+wid)))
}

#' @rdname wilsonCI
#' @export
jeffCI<-function(phat,n,conf=0.9,w1=0.5,w2=w1,...)
{
x=n*phat
clow=ifelse(phat==0,0,qbeta((1-conf)/2,x+w1,n-x+w2))
chigh=ifelse(phat==1,1,qbeta(1-(1-conf)/2,x+w1,n-x+w2))

return(cbind(clow,chigh))
}


########### Interpolation, Extrapolation, Parapolation....

# Linearly interpolate/extrapolate from a single segment.

extrapol<-function(point1,point2,xout)
{
slope=(point2[2]-point1[2])/(point2[1]-point1[1])
return(point2[2]+(xout-point2[1])*slope)
}

# Locally monotone quadratic interpolation between points on a plane

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

slope<-function(x,y,xout,allowZero=FALSE,full=FALSE)
{
### Validation (might be mostly redundant if using doseResponse as input)

if (any(xout>max(x) | xout<min(x))) stop("No extrapolation allowed in 'slopes'.\n")
m=length(x)
if(length(y)!=m) stop("Mismatched lengths in 'slopes'.\n")
xdiffs=diff(x)
ydiffs=diff(y)
if (any(xdiffs<=0 | ydiffs<0)) stop("Monotonicity violation in 'slopes'.\n")
slopes=ydiffs/xdiffs
sslopes=c(0,slopes,0)

interval=findInterval(xout,x)
## The trivial ones
candidate=slopes[interval]
## ones falling on design points
design=which(xout %in% x)
if (length(design)>0) {
	for(a in design) candidate[a]=(sslopes[interval[a]]+sslopes[interval[a]+1])/2
}
candidate0=candidate

## tougher nut: zero slope
if(!allowZero && any(candidate==0))
{
	xstep=mean(xdiffs)/2
	for (a in which(candidate==0))
	{	
		b=0
		while(candidate[a]==0)
		{
			b=b+1
			xends=c(max(x[1],xout[a]-b*xstep),min(x[m],xout[a]+b*xstep))
			yends=approx(x,y,xout=xends)$y
			candidate[a]=diff(yends)/diff(xends)
		}
	}
}
if(!full) return (candidate)
return(list(scrappy=candidate,clean=candidate0,rawslopes=slopes))		
}
