##' Returns analytic interval estimates, given isotonic-regression (centered or not) point estimates
#'
#'
#' Forward interval estimation of a monotone response (y) as a function of dose (x), given a set of point estimates produced by  \code{\link{cirPAVA}} or \code{\link{oldPAVA}}.
#'
#' @seealso \code{\link{quickIsotone}},\code{\link{quickInverse}},\code{\link{morrisCI}},
#'
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}
#' @export

#' @param isotPoint a \code{\link{doseResponse}} object with the x values and isotonic point estimates at design points.
#' @param outx vector of x values for which estimates will be made. If \code{NULL} (default), this will be set to the set of unique values in isotPoint$x argument (or equivalently in y$x).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param intfun the function to be used for interval estimation. Default \code{\link{morrisCI}} (see help on that function for additional options).
#' @param logical, should a rough accounting for the added randomness due to the use of a sequential dose-finding design? Default \code{FALSE} due to little benefit.
#' @param logical, should the interpolation between design points follow a parabola (\code{TRUE}, default) or a straight line? See details.
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




