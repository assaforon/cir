##' Returns analytic interval estimates, given isotonic-regression (centered or not) point estimates
#'
#'
#' Forward interval estimation of a monotone response (y) as a function of dose (x), given a set of point estimates produced by \code{\link{pava}} or \code{\link{cirPAVA}}.
#'
#'
#'
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}

#' @param isotPoint a \code{\link{doseResponse}} object with the x values and isotonic point estimates at design points.
#' @param outx vector of x values for which estimates will be made. If \code{NULL} (default), this will be set to the set of unique values in isotPoint$x argument (or equivalently in y$x).
#' @param intfun the function to be used for interval estimation. Default \code{\link{wilsonCI}} (see help on that function for additional options).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.

isotInterval<-function(isotPoint,outx=isotPoint$x,conf=0.9,intfun=wilsonCI)
{
## Validation
if(conf<=0 || conf>=1) stop("Confidence must be between 0 and 1.\n")
if(!is.doseResponse(isotPoint)) stop("Point-estimate data must be in doseResponse format.\n")
if(min(outx)<min(isotPoint$x) || max(outx)>max(isotPoint$x)) stop("Cannot predict outside design boundaries.\n")

designInt=intfun(phat=isotPoint$y,n=isotPoint$weight,conf=conf)

#if(all(outx %in% isotPoint$x)) return(designInt[match(outx,isotPoint$x),])

return(data.frame(ciLow=approx(isotPoint$x,designInt[,1],xout=outx)$y
	,ciHigh=approx(isotPoint$x,designInt[,2],xout=outx)$y))
}

#'
#' One-Step Forward point and interval estimation via IR or CIR
#'
#' One-Step Forward point and confidence-interval estimation of a monotone response (y) as a function of dose (x). Input format is rather flexible.
#'
#'
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}

#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). Note that the PAV algorithm doesn't really use them. 
#' @param wt weights (if not included in y).
#' @param outx vector of x values for which predictions will be made. If \code{NULL} (default), this will be set to the set of unique values in the x argument (or equivalently in y$x).
#' @param dec logical, is the true function is assumed to be monotone decreasing? Default \code{FALSE}.
#' @param cir logical, is centered-isotonic-regression (CIR) to be used? If \code{FALSE}, traditional isotonic regression is used. Default \code{TRUE}.
#' @param intfun the function to be used for interval estimation. Default \code{\link{wilsonCI}} (see help on that function for additional options).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param ...	Other arguments passed on to the estimating function.

quickIsotone<-function (y,x=NULL,wt=rep(1,length(y)),outx=NULL,dec=FALSE,cir=TRUE,intfun=wilsonCI,conf=0.9,...) 
{
dr=doseResponse(y=y,x=x,wt=wt,...)
if(is.null(outx)) outx=dr$x

if (cir) {
	pestimate=cirPAVA(y=dr,dec=dec,full=TRUE)
} else pestimate=oldPAVA(y=dr,dec=dec,full=TRUE)

cestimate=isotInterval(pestimate$output,conf=conf,intfun=intfun,outx=outx)

if(all(outx %in% dr$x)) 
{
	dout=cbind(pestimate$output[,1:2],cestimate)
} else
{
	dout=data.frame(x=outx,y=approx(pestimate$output$x,pestimate$output$y,xout=outx)$y,low=cestimate[,1],hi=cestimate[,2])
}
names(dout)[3:4]=paste(c("lower","upper"),round(100*conf),"conf",sep="")
return(dout)
}




