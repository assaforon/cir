#### Some standard CI calculation functions
#' Standard unordered-Binomial confidence interval utilitlies.
#' 
#' Standard small-sample Binomial confidence interval utilities, using the methods of Wilson, Agresti-Coull and Jeffrys.
#' 
#' These functions implement the basic (uncorrected) three intervals which are seen by the consensus of literature as the "safest" off-the-shelf formulae. None of them account for ordering or monotonicity; therefore the \code{cir} package default is \code{\link{morrisCI}} which does account for that, with the 3 unordered formulae used for optional narrowing of the interval at individual points.

#' @aliases agcouCI jeffCI
#' @seealso \code{\link{isotInterval}} for more details about how forward CIs are calculated, \code{\link{quickInverse}} for inverse (dose-finding) intervals.
#' @export

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


