##' Returns analytical interval estimates, given isotonic-regression (centered or not) point estimates
#'
#'
#' For confidence intervals at design points ($x$ values with obesrvations), this function calls \code{intfun} to do the work. In addition, CIs for any $x$ value are calculated using linear interpolation between design points (note that for CIR, this differs from the interpolation of point estimates which is carried out between shrinkage points, as explained in \code{\link{quickIsotone}})
#' 
#' @note All provided algorithm and formulae are for Binomial data only. For other data, write your own \code{intfun}, returning a two-column matrix. The interval estimation method is presented and discussed by Oron and Flournoy (2017).
#'
#' @warning As of mid-2021, we've found that interval coverage for extreme percentiles with adaptive designs is relatively lacking: by 10-15% for ED90/ED10 type targets, and deteriorating further towards the edges. 
#'
#' @seealso \code{\link{quickIsotone}},\code{\link{quickInverse}},\code{\link{morrisCI}},
#' @example inst/examples/fwdCiExamples.r
#'
##' @author Assaf P. Oron \code{<aoron.at.idmod.org>}
#' @export
#' @references Oron, A.P. and Flournoy, N., 2017. Centered Isotonic Regression: Point and Interval Estimation for Dose-Response Studies. Statistics in Biopharmaceutical Research 3, 258-267.

#' @return a data frame with two variables \code{ciLow, ciHigh} containing the estimated lower and upper confidence bounds, respectively.
#' 
#' @param isotPoint The output of an estimation function such as \code{\link{cirPAVA}}  with the option \code{full=TRUE}. Should be a list of 3 \code{\link{doseResponse}} objects named \code{input, output, shrinkage}.
#' @param outx vector of x values for which estimates will be made. If \code{NULL} (default), this will be set to the set of unique values in isotPoint$x argument (or equivalently in y$x).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param intfun the function to be used for interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param ... additional arguments passed on to \code{intfun}

isotInterval<-function(isotPoint,outx=isotPoint$output$x,conf=0.9,intfun=morrisCI,...)
{
## Validation
if(conf<=0 || conf>=1) stop("Confidence must be between 0 and 1.\n")
#if(!is.doseResponse(isotPoint)) stop("Point-estimate data must be in doseResponse format.\n")
if(min(outx)<min(isotPoint$output$x) || max(outx)>max(isotPoint$output$x)) stop("Cannot predict outside data boundaries.\n")

ycount=round(isotPoint$shrinkage$weight*isotPoint$shrinkage$y)
rawInt=intfun(phat=isotPoint$shrinkage$y,n=isotPoint$shrinkage$weight,y=ycount,conf=conf,...)

#if(all(outx %in% isotPoint$x)) return(designInt[match(outx,isotPoint$x),])

if(length(unique(rawInt[,1]))==1 || length(unique(rawInt[,2]))==1)
{ # degenerate case: only one y value
	lcl=rep(NA,length(outx))
	ucl=rep(NA,length(outx))
} else {
	lcl=approx(isotPoint$shrinkage$x,rawInt[,1],xout=outx)$y
	ucl=approx(isotPoint$shrinkage$x,rawInt[,2],xout=outx)$y
}
return(data.frame(ciLow=lcl,ciHigh=ucl))
}

####################### Inverse

#' Calculate inverse (dose-finding) intervals, using local inversion and the Delta method
#'
#'
#' Calculate left-bound to right-bound intervals for the dose point estimates, using local slopes at design points (places where observations exist) to invert the forward lower-upper bounds.
#'
#'
#' The Delta method in this application boils down to dividing the distance to the forward (vertical) bounds, by the slope, to get the left/right interval width. Slope estimates are performed by \code{\link{slope}}. An alternative method (dubbed "global") is hard-coded into \code{\link{quickInverse}}. 

#' @return two-column matrix with the left and right bounds, respectively

#' @param isotPoint The output of an estimation function such as \code{\link{cirPAVA},\link{doseFind}},  with the option \code{full=TRUE}. Should be at least a list of 3 \code{\link{doseResponse}} objects named \code{input, output, shrinkage}.
#' @param target A vector of target response rate(s), for which the interval is needed. If \code{NULL} (default), interval will be returned for the point estimates at design points (e.g., if the forward point estimate at $x_1$ is 0.2, then the first returned interval is for the 20th percentile).
#' @param intfun the function to be used for initial (forward) interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param starget The shrinkage target. Defaults to \code{target[1]}.
#' @param ... additional arguments passed on to \code{\link{quickIsotone}}

#' @seealso \code{\link{quickIsotone}},\code{\link{quickInverse}},\code{\link{isotInterval}},
#' \code{\link{slope}}; \code{\link{DRshrink}} for the shrinkage fix.
#' @example inst/examples/invCiExamples.r


#' @export

deltaInverse<-function(isotPoint,target=NULL,intfun = morrisCI, conf = 0.9,starget=target[1],...)
{
k=length(target)
#isotPoint$shrinkage$y=round(isotPoint$shrinkage$y,10) ### avoid rounding errors from PAVA
yvals=sort(unique(isotPoint$shrinkage$y))
#cat(yvals)
if(is.null(yvals) || length(yvals)<=1 || var(yvals)<.Machine$double.eps*1e3) return(cbind(rep(NA,k),rep(NA,k))) ## degenerate case, completely flat

cestimate=isotInterval(isotPoint,conf=conf,intfun=intfun,outx=isotPoint$shrinkage$x,...)
fslopes=slope(isotPoint$shrinkage$x,isotPoint$shrinkage$y)

# inverse widths raw
rwidths=(isotPoint$shrinkage$y-cestimate$ciLow)/fslopes
lwidths=(isotPoint$shrinkage$y-cestimate$ciHigh)/fslopes
# Adding the widths to the mean curve, self-consistently
rbounds=rev(cummin(rev(tapply(isotPoint$shrinkage$x+rwidths,isotPoint$shrinkage$y,max))))
lbounds=cummax(tapply(isotPoint$shrinkage$x+lwidths,isotPoint$shrinkage$y,min))

### Returning
# Note we use approx() with rule=1, forcing NAs when specified target is outside bounds

if(length(unique(lbounds))==1 || length(unique(rbounds))==1)
{ # degenerate case: only one y value. No interval can be calculated
	nout=ifelse(is.null(target),nrow(isotPoint$output),length(target))
	lout=rep(NA,nout)
	rout=rep(NA,nout)
	return(cbind(lout,rout))
}	
if (is.null(target)) 
{ ## No target specified, returning CIs at design points
	lout = approx(isotPoint$shrinkage$y,lbounds,isotPoint$output$y,rule=1)$y
	uout = approx(isotPoint$shrinkage$y,rbounds,isotPoint$output$y,rule=1)$y
	return(cbind(lout,uout))
}
# Otherwise: target was specified
lout=approx(yvals,lbounds,xout=target,rule=1)$y
rout=approx(yvals,rbounds,xout=target,rule=1)$y

return(cbind(lout,rout))
}



