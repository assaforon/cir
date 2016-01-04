##' Returns analytic interval estimates, given isotonic-regression (centered or not) point estimates
#'
#'
#' For confidence intervals at design points ($x$ values with obesrvations), this function calls \code{intfun} to do the work. In addition, CIs for any $x$ value are calculated. Naively, one would do a linear interpolation just like with the point estimate. However, between observations the uncertainty is greater, and it can also be generally assumed that the curve is nonlinear. Therefore, a piecewise-convex (for the LCL) or piecewise-concave (for the UCL) confidence band seems in order, analogous to the curvilinear bands of linear regression. 
#' 
#' The solution used was chosen because it does not require tuning parameters, while maintaining monotonicity. Each segment is a parabola whose slope becomes zero at one endpoint. The requirement of maintaining convexity/concavity for LCL/UCL determines that endpoint. This is implemented via the internal function \code{parapolate}
#' 
#' @note All provided algorithm and formulae are for Binomial data only. For other data, write your own \code{intfun}, returning a two-column matrix. Also, interval estimate for dose-finding (inverse estimation) is only available via the \code{\link{quickInverse}} function. That function calls \code{isotInterval} and then reads the inverse interval "horizontally". See details there.
#'
#' @seealso \code{\link{quickIsotone}},\code{\link{quickInverse}},\code{\link{morrisCI}},
#'
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}
#' @export

#' @return a data frame with two variables \code{ciLow, ciHigh} containing the estimated lower and upper confidence bounds, respectively.
#' 
#' @param isotPoint a \code{\link{doseResponse}} object with the x values and isotonic point estimates at design points.
#' @param outx vector of x values for which estimates will be made. If \code{NULL} (default), this will be set to the set of unique values in isotPoint$x argument (or equivalently in y$x).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param intfun the function to be used for interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param sequential logical, should a rough accounting for the added randomness due to the use of a sequential dose-finding design? Default \code{FALSE} due to little benefit.
#' @param parabola logical, should the interpolation between design points follow a parabola (\code{TRUE}, default) or a straight line? See details.
#' @param ... additional arguments passed on to \code{intfun}

isotInterval<-function(isotPoint,outx=isotPoint$x,conf=0.9,intfun=morrisCI,sequential=FALSE,parabola=TRUE,...)
{
## Validation
if(conf<=0 || conf>=1) stop("Confidence must be between 0 and 1.\n")
if(!is.doseResponse(isotPoint)) stop("Point-estimate data must be in doseResponse format.\n")
if(min(outx)<min(isotPoint$x) || max(outx)>max(isotPoint$x)) stop("Cannot predict outside data boundaries.\n")

ycount=round(isotPoint$weight*isotPoint$y)
designInt=intfun(phat=isotPoint$y,n=isotPoint$weight,y=ycount,conf=conf,...)

if(sequential) ## correction for sequential designs
{
	pi_j=isotPoint$weight/sum(isotPoint$weight)
	corfac=sqrt(1+(1-pi_j)/isotPoint$weight)

	newidths=(designInt-isotPoint$y)*as.vector(corfac)
	designInt=newidths+isotPoint$y
	designInt[designInt<0]=0
	designInt[designInt>1]=1
}
#if(all(outx %in% isotPoint$x)) return(designInt[match(outx,isotPoint$x),])

if(parabola) lcl=parapolate(isotPoint$x,designInt[,1],xout=outx,upward=TRUE) else
	lcl=approx(isotPoint$x,designInt[,1],xout=outx)$y
if(parabola) ucl=parapolate(isotPoint$x,designInt[,2],xout=outx,upward=FALSE) else
	ucl=approx(isotPoint$x,designInt[,2],xout=outx)$y

return(data.frame(ciLow=lcl,ciHigh=ucl))
}

####################### Inverse

#' Calculate inverse (dose-finding) intervals, using local inversion and the Delta method

#' Calculate left-bound to right-bound intervals for the dose point estimates, using local slopes at design points (places where observations exist) to invert the forward lower-upper bounds.

#' The Delta method in this application boils down to dividing the distance to the forward (vertical) bounds by the slope. Slope estimates are subcontracted to \code{\link{slope}}.  

#' @return two-column matrix with the left and right bounds, respectively

#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). 
#' @param wt weights (if not included in y).
#' @param target A vector of target response rate(s), for which the interval is needed. If \code{NULL} (default), interval will be returned for the point estimates at design points (e.g., if the forward point estimate at $x_1$ is 0.2, then the first returned interval is for the 20th percentile).
#' @param cir logical, is centered-isotonic-regression (CIR) to be used? If \code{FALSE}, traditional isotonic regression is used. Default \code{TRUE}.
#' @param intfun the function to be used for interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.

#' @export

deltaInverse<-function(y,x=NULL,wt=NULL,target=NULL,cir = TRUE, intfun = morrisCI, conf = 0.9,seqDesign=FALSE,parabola=FALSE,...)
{
dr=doseResponse(y,x,wt)
k=length(target)
# We start by constructing inverse intervals based on design-point estimates
forward=quickIsotone(dr,outx=NULL,conf=conf,intfun=intfun,seqDesign=seqDesign,cir=cir,...)
forward$y=round(forward$y,10) ### avoid rounding errors from PAVA
yvals=sort(unique(forward$y))
#cat(yvals)
if(length(yvals)==1 || var(yvals)<.Machine$double.eps*1e3) return(cbind(rep(NA,k),rep(NA,k))) ## degenerate case, completely flat
fslopes=slope(dr$x,forward$y)

# inverse widths raw
rwidths=(forward$y-forward$lower)/fslopes
lwidths=(forward$y-forward$upper)/fslopes

rbounds=rev(cummin(rev(tapply(dr$x+rwidths,forward$y,max))))
lbounds=cummax(tapply(dr$x+lwidths,forward$y,min))
designCIs=cbind(lbounds,rbounds)[match(forward$y,yvals),]

if (is.null(target)) return(designCIs)

if(parabola) {  ### direction chosen by trial & error...
	good=(target>=min(yvals) & target<=max(yvals))
	lout=rep(NA,k)
	rout=lout
	lout[good]=parapolate(yvals,lbounds,xout=target[good],upward=TRUE)
	rout[good]=parapolate(yvals,rbounds,xout=target[good],upward=FALSE)
} else {
	lout=approx(yvals,lbounds,xout=target,rule=1)$y
	rout=approx(yvals,rbounds,xout=target,rule=1)$y
}
return(cbind(lout,rout))
}



