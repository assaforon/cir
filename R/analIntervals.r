##' Returns analytic interval estimates, given isotonic-regression (centered or not) point estimates
#'
#'
#' Forward interval estimation of a monotone response (y) as a function of dose (x), given a set of point estimates produced by \code{\link{pava}} or \code{\link{cirPAVA}}.
#'
#'
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}

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

quickIsotone<-function (y,x=NULL,wt=rep(1,length(y)),outx=NULL,dec=FALSE,cir=TRUE,intfun=wilsonCI,conf=0.9,...) 
{
dr=doseResponse(y=y,x=x,wt=wt,...)
if(is.null(outx)) outx=dr$x

if (cir) {
	pestimate=cirPAVA(y=dr,dec=dec,full=TRUE)
} else pestimate=drPAVA(y=dr,dec=dec,full=TRUE)

cestimate=isotInterval(pestimate$output,conf=conf,intfun=intfun,outx=outx)

if(all(outx %in% dr$x)) 
{
	dout=cbind(pestimate$output[,1:2],cestimate)
} else
{
	dout=data.frame(x=outx,y=approx(pestimate$output$x,pestimate$output$y,xout=outx)$y,low=cestimate[,1],hi=cestimate[,2])
}
names(dout)[3:4]=paste(c("Lower","Upper"),round(100*conf),"conf",sep="")
return(dout)
}




