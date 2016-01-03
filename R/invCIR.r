##' Inverse (dose-finding) estimate of a target x value (e.g., a percentile)
#'
#'
#' Inverse ("dose-finding") point estimation of a dose (x) for a specified target y value (e.g., a response rate), using centered-isotonic-regression (\code{invCIR}) or a generic forward-estimation algorithm (\code{doseFind}).
#'
#'
#' The function works by calling \code{estfun} for forward estimation of the x-y relationship, then using \code{\link{approx}} with the x and y roles reversed for inverse estimation. The \code{extrapolate} option sets the \code{rule} argumet for this second call: 
#' \itemize{
#' \item {}{\code{extrapolate=TRUE} translates to \code{rule=2}, which actually means that the x value on the edge of the estimated y range will be assigned.}
#' \item{}{\code{extrapolate=FALSE} (default) translates to \code{rule=1}, which means an \code{NA} will be returned for any target y value lying outside the estimated y range.}
#' }
#' Note also that the function is set up to work with a vector of targets.
 
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}

#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). 
#' @param wt weights (if not included in y).
#' @param target A vector of target response rate(s), for which the percentile dose estimate is needed.
#' @param full logical, is a more complete output desired (relevant only for doseFind)? if \code{FALSE} (default), only a point estimate of the dose (x) for the provided target rate is returned
#' @param extrapolate logical: should extrapolation beyond the range of estimated y values be allowed? Default \code{FALSE}.
#' @param dec (relevant only for doseFind) logical, is the true function is assumed to be monotone decreasing? Default \code{FALSE}.
#' @param estfun the name of the dose-response estimation function (relevant only for doseFind). Default \code{\link{cirPAVA}}.
#' @param ...	Other arguments passed on to \code{\link{doseResponse}} and \code{estfun}.

#' @return under default, returns point estimate(s) of the dose (x) for the provided target rate(s). With \code{full=TRUE}, returns a list with
#' \itemize{
#' \item {xout } {  The said point estimate of x}
#' \item {input  }  {  a \code{doseResponse} object summarizing the input data}
#' \item {cir  }  {  a \code{doseResponse} object which is the \code{alg} output of the forward-estimation function}
#' }

#' @seealso \code{\link{oldPAVA}},\code{\link{cirPAVA}}. If you'd like point and interval estimates together, use \code{\link{quickInverse}}.

#' @export

doseFind<-function(y,x=NULL,wt=NULL,estfun=cirPAVA,target=NULL,full=FALSE,dec=FALSE,extrapolate=FALSE,...) {


### converting to doseResponse object 
### Basically it's a numeric data frame with x,y,weight, and x increasing
if(is.null(target)) stop("Must provide target to dose-find for.\n")
dr=doseResponse(y=y,x=x,wt=wt,...)
if (any(is.na(dr))) stop ("Missing values are not allowed.\n")  

# We start via forward estimation

pavout<-estfun(y=dr,full=TRUE,dec=dec,...)

newx<-pavout$alg$x
newy<-pavout$alg$y
newn=pavout$alg.wt

### The estimate is generated by using 'approx' with x and y interchanged

tout=approx(x=newy,y=newx,xout=target,ties="ordered",rule=1)$y
# targets outside y range are provisionally NA'ed... now if so desired, extrapolate
if(any(is.na(tout)) && extrapolate)
{
	lowTargets=which(target<min(newy))
	if(length(lowTargets)>0) tout[lowTargets]=extrapol(c(newy[1],newx[1]),c(newy[2],newx[2]),xout=target[lowTargets])
	hiTargets=which(target>max(newy))
	m=length(newx)
	if(length(hiTargets)>0) tout[hiTargets]=extrapol(c(newy[m-1],newx[m-1]),c(newy[m],newx[m]),xout=target[hiTargets])
	tout[!is.finite(tout)]=NA
}
	
#### output

if (!full)  return(tout) 

return (list(targest=tout,input=dr,fwd=pavout$alg,fwdDesign=pavout$output))
}

#' Point and Interval Inverse Estimation ("Dose-Finding"), using CIR and IR
#'
#'
#' Convenience wrapper for point and interval estimation of the "dose" that would generate a \code{target} "response" value, using CIR and IR.
#' 
#' 
#' The inverse point estimate is calculated in a straightforward manner from a forward estimate, using \code{\link{doseFind}}. For the inverse interval, \code{\link{quickIsotone}} is called for forward estimation on a high-resolution grid of $x$ values (the default interval curves are piecewise-parabolic), then the closest points to the line \code{y=target} are found. Note that the forward UCL is used to find the inverse LCL, and vice versa. 
#' 
#' @return A data frame with 4 variables: \code{x} either the input x values, or \code{outx} of specified; \code{y} the point estimates; and the lower and upper confidence bounds.


#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). 
#' @param wt weights (if not included in y).
#' @param target A vector of target response rate(s), for which the percentile dose estimate is needed.
#' #' @param cir logical, is centered-isotonic-regression (CIR) to be used? If \code{FALSE}, traditional isotonic regression is used. Default \code{TRUE}.
#' @param intfun the function to be used for interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param resolution numeric: how fine should the grid for the inverse-interval approximation be? Default 100, which seems to be quite enough. See 'Details'.
#' @param xbounds numeric vector of 2, lower and upper bounds for the confidence intervals, beyond which the interval "doesn't make sense". Under the default (\code{NULL}), the function will set these to one spacing level outside the boundaries of \code{x}.
#' @param extrapolate logical: should extrapolation beyond the range of estimated y values be allowed? Default \code{FALSE}. Note this affects only the point estimate; interval boundaries are extrapolated in any case.
#' @param seqDesign logical, should intervals be further widened using a simple adjustment for the data having been obtained via a sequential (adaptive) design? Default \code{FALSE} due to futility.
#' @param ...	Other arguments passed on to \code{\link{doseFind}} and \code{\link{quickIsotone}}, and onwards from there.

#' @return A data frame with 4 elements:
#' \itemize{
#' \item {\code{target}  }  { The user-provided target values of y, at which x is estimated}
#' \item {\code{point} } {  The point estimates of x}
#' \item {\code{lowerPPconf,upperPPconf}  }  { the interval-boundary estimates for a 'PP'=\code{100*conf} confidence interval}
#' }

#' @seealso \code{\link{quickIsotone}},\code{\link{doseFind}},\code{\link{isotInterval}}
#' @export

quickInverse<-function(y,x=NULL,wt=NULL,target,cir = TRUE, intfun = morrisCI,delta=TRUE,conf = 0.9,resolution=100,xbounds=NULL,extrapolate=FALSE,seqDesign=FALSE,...)
{

#### Point estimate first
dr=doseResponse(y,x,wt)
m=length(dr$x)
if(cir) estfun<-cirPAVA else estfun<-oldPAVA
pestimate=doseFind(y=dr,estfun=estfun,target=target,full=TRUE,extrapolate=extrapolate,...) 
foundPts=pestimate$targest[!is.na(pestimate$targest)]

dout=data.frame(target=target,point=pestimate$targest,low=-Inf,high=Inf)

if(delta) { ## new default, delta-method ("local") intervals
	dout[,3:4]=deltaInverse(y=dr,target=target,cir=cir,intfun=intfun,conf=conf,seqDesign=seqDesign,...)
} else {
#### CI, using "global" interpolation, a more stable method

	## Calculate forward CIs at high-rez grid
	higrid=seq(dr$x[1],dr$x[m],length.out=1+resolution*(m-1))
	fwdCIgrid=quickIsotone(dr,outx=higrid,conf=conf,intfun=intfun,seqDesign=seqDesign,...)
	fwdCIdesign=fwdCIgrid[match(dr$x,fwdCIgrid$x),]

	### Setting extrapolation boundaries for CI. If not given by user, we establish "logical" boundaries for CIs, at one spacing level out

	if(is.null(xbounds)) xbounds=rep(NA,2)
	xmin=dr$x[1]-(dr$x[2]-dr$x[1])
	xmax=dr$x[m]+(dr$x[m]-dr$x[m-1])
	xbounds[1]=min(xmin,xbounds[1],na.rm=TRUE)
	xbounds[2]=max(xmax,xbounds[2],na.rm=TRUE)
	maxLower=max(c(fwdCIgrid[,3],extrapol(point1=c(dr$x[1],fwdCIdesign[1,3]),point2=c(dr$x[m],fwdCIdesign[m,3]),xout=xbounds[2])))
	minUpper=min(c(fwdCIgrid[,4],extrapol(point1=c(dr$x[1],fwdCIdesign[1,4]),point2=c(dr$x[m],fwdCIdesign[m,4]),xout=xbounds[1])))

	### CI by finding 'right points' on grid
	for (a in 1:length(dout$target))
	{	
		if(dout$target[a]<=max(fwdCIgrid[,3]))  # upper inverse interval taken from fwd LCL 
		{
			dout$high[a]=fwdCIgrid[min(which(fwdCIgrid[,3]>=dout$target[a])),1]
		} else if(dout$target[a]<=maxLower)  {
			dout$high[a]=extrapol(point1=c(fwdCIdesign[1,3],fwdCIdesign$x[1]),
				point2=c(fwdCIdesign[m,3],fwdCIdesign$x[m]),xout=dout$target[a])
		}
		if(dout$target[a]>=min(fwdCIgrid[,4]))  # and vice versa
		{
			dout$low[a]=fwdCIgrid[max(which(fwdCIgrid[,4]<=dout$target[a])),1]
		} else if(dout$target[a]>=minUpper)  {
			dout$low[a]=extrapol(point1=c(fwdCIdesign[1,4],fwdCIdesign$x[1]),
				point2=c(fwdCIdesign[m,4],fwdCIdesign$x[m]),xout=dout$target[a])
		}
	}
}

#dout$low[!is.na(pestimate$targest)]=ciLow
#dout$high[!is.na(pestimate$targest)]=ciHigh

names(dout)[3:4]=paste(c("lower","upper"),round(100*conf),"conf",sep="")
return(dout)
}






