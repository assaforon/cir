##' Returns centered-isotonic-regression estimate
#'
#'
#' Nonparametric forward point estimation of a monotone response (y) as a function of dose (x), using the centered-isotonic-regression (CIR) algorithm.
#'
#'
#' This is the underlying "engine" function implementing CIR. For a quick and somewhat more user-friendly wrapper, use \code{\link{quickIsotone}}. CIR is a variation of isotonic regression (IR) that shrinks IR's constant ("flat") intervals to single points and interpolates between these points, generating a curve that is strictly monotone everywhere except (possibly) near the boundaries.
#'
#' Flat intervals in the raw input data, are handled with care. Under the default setting (\code{strict=FALSE,interiorStrict=TRUE}), flat intervals are treated as monotonicity violations, unless the $y$ value is on the boundary of its allowed range (default $[0,1]$, appropriate for binary-response data). On that boundary, flat intervals are left unchanged.

##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}


#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). 
#' @param wt weights (if not included in y).
#' @param outx vector of x values for which predictions will be made. If \code{NULL} (default), this will be set to the set of unique values in the x argument (or equivalently in y$x). Non-NULL inputs are relevant only if \code{full=TRUE}.
#' @param full: logical, is a more complete output desired? if \code{FALSE} (default), only a vector of point estimates for y at the provided dose levels is returned
#' @param dec logical, is the true function is assumed to be monotone decreasing? Default \code{FALSE}.
#' @param strict logical, should CIR enforce strict monotonicity by "fixing" flat intervals as well? Default \code{FALSE}.
#' @param interiorStrict logical, should CIR enforce strict monotonicity, but only for y values inside of \code{ybounds}?  Default \code{TRUE}. Choosing \code{FALSE} will be overridden if \code{strict=TRUE}, and a warning will be given.
#' @param ybounds numeric vector of length 2, relevant only under the default setting of \code{strict=FALSE,interiorStrict=TRUE}. Default \code{0:1}. See 'Details'.
#' @param ...	Other arguments passed on to the constructor functions that pre-process the input.

#' @return under default, returns a vector of y estimates at unique x values. With \code{full=TRUE}, returns a list of 3 \code{\link{doseResponse}} objects name \code{output,input,alg} for the output data at dose levels, the input data, and the function as fit at algorithm-generated points, respectively.

#' @seealso \code{\link{oldPAVA}},\code{\link{quickIsotone}}
#' @export
#' 
cirPAVA <-function (y,x=NULL,wt=NULL,outx=NULL,full=FALSE,dec=FALSE,strict=FALSE,interiorStrict=TRUE,ybounds=0:1,...) {

### converting to doseResponse object 
### Basically it's a numeric data frame with x,y,weight, and x increasing

dr=doseResponse(y=y,x=x,wt=wt,...)
if (any(is.na(dr))) stop ("Missing values are not allowed.\n")  

### Predictions will be delivered for x=outx
if(is.null(outx)) outx=dr$x

### validations
if(min(outx)<min(dr$x) || max(outx)>max(dr$x)) stop("Cannot predict outside design boundaries.\n")
if(strict && !interiorStrict) warning("strict=TRUE overrides interiorStrict=FALSE.\n")

m <- dim(dr)[1]
if (m <= 1) {  ## degenerate case: only one dose level
if (!full) return (dr$y)
else return(list(output=dr,input=dr,alg=dr))
}


dr0=dr ## clean copy of input data
### Decreasing monotone case: simple fix
if (dec) dr$y = -dr$y


#### Core algorithm

repeat {

# Find adjacent violators
# Definition of violators comes in 1 of 3 'flavors', see help
	viol <- (as.vector(diff(dr$y)) < 0)
	if(interiorStrict) 
	{
		addviol <- (as.vector(diff(dr$y)) == 0)
		addviol[addviol & (dr$y[-m] %in% ybounds) & (dr$y[-1] %in% ybounds)]=FALSE
		viol<-(viol | addviol)
	}
	if(strict) viol <- (as.vector(diff(dr$y)) <= 0)

    if (!(any(viol))) break
    i <- min( (1:(m-1))[viol]) # Pool first pair of violators
    dr$y[i] <- (dr$y[i]*dr$weight[i]+dr$y[i+1]*dr$weight[i+1]) / (dr$weight[i]+dr$weight[i+1])
    dr$x[i] <- (dr$x[i]*dr$weight[i]+dr$x[i+1]*dr$weight[i+1]) / (dr$weight[i]+dr$weight[i+1])  # new x is calculated
    dr$weight[i]<-dr$weight[i]+dr$weight[i+1]  # weights are combined

# Deleting the i-1-th element and updating n
    dr<-dr[-(i+1),]
    m <- dim(dr)[1]
    if (m <= 1) break
  }

# Now, if relevant we re-interpolate to original x boundaries, stored in z
# This will give constant y values for x values falling outside new range of x
# Which is identical to the PAVA solution on those ranges

## Extending back to original boundaries if needed
if(dr$x[1]>dr0$x[1] ) 
{
	dr=rbind(doseResponse(x=dr0$x[1],y=dr$y[1],wt=0),dr)
}
if(max(dr$x)<max(dr0$x)) ## Extending back to original boundaries if needed
{
	dr=rbind(dr,doseResponse(x=max(dr0$x),y=max(dr$y),wt=0))
}

if (dec) dr$y = -dr$y

outy=approx(dr$x,dr$y,outx,rule=2)$y

if (!full) {
	return(outy)
} else {
	
	if(all(outx %in% dr0$x)) {
		dr1=dr0
		dr1$y=outy
		dr1=dr1[match(outx,dr1$x),]
	} else dr1=doseResponse(y=outy,x=outx,wt=rep(0,length(outy)))
	return(list(output=dr1,input=dr0,alg=dr))   }
}

#'
#' One-Step Forward point and interval estimation via CIR or IR
#'
#' One-Step Forward point and confidence-interval estimation of a monotone response (y) as a function of dose (x), using centered-isotonic-regression (CIR, default) or isotonic regression. Input format is rather flexible.

#' This function calls either \code{\link{cirPAVA}} or \code{\link{oldPAVA}} for the point estimate, then \code{\link{isotInterval}} for the confidence interval. Vector input is allowed, but the preferred input format is a \code{\link{doseResponse}} object.

#' An analogous function for dose-finding (inverse estimation) is \code{\link{quickInverse}}.
#'
#'
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}
#' @export

#' @seealso \code{\link{cirPAVA}},\code{\link{oldPAVA}},\code{\link{quickInverse}},\code{\link{doseResponse}}

#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). Note that the PAV algorithm doesn't really use them. 
#' @param wt weights (if not included in y).
#' @param outx vector of x values for which predictions will be made. If \code{NULL} (default), this will be set to the set of unique values in the x argument (or equivalently in y$x).
#' @param dec logical, is the true function is assumed to be monotone decreasing? Default \code{FALSE}.
#' @param cir logical, is centered-isotonic-regression (CIR) to be used? If \code{FALSE}, traditional isotonic regression is used. Default \code{TRUE}.
#' @param intfun the function to be used for interval estimation. Default \code{\link{wilsonCI}} (see help on that function for additional options).
#' @param conf numeric, the interval's confidence level as a fraction in (0,1). Default 0.9.
#' @param ...	Other arguments passed on to the estimating function.

quickIsotone<-function (y,x=NULL,wt=rep(1,length(y)),outx=NULL,dec=FALSE,cir=TRUE,intfun=morrisCI,conf=0.9,seqDesign=FALSE,...) 
{
dr=doseResponse(y=y,x=x,wt=wt,...)
if(is.null(outx)) outx=dr$x

if (cir) {
	pestimate=cirPAVA(y=dr,dec=dec,full=TRUE,...)
} else pestimate=oldPAVA(y=dr,dec=dec,full=TRUE,...)

cestimate=isotInterval(pestimate$output,conf=conf,intfun=intfun,outx=outx,sequential=seqDesign,...)

if(all(outx %in% dr$x)) 
{
	dout=cbind(pestimate$output[match(outx,pestimate$output$x),1:2],cestimate)
} else
{
	dout=data.frame(x=outx,y=approx(pestimate$output$x,pestimate$output$y,xout=outx)$y,low=cestimate[,1],hi=cestimate[,2])
}
names(dout)[3:4]=paste(c("lower","upper"),round(100*conf),"conf",sep="")
return(dout)
}

