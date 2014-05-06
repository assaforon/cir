##' Returns standard isotonic-regression estimate
#'
#'
#' Nonparametric forward point estimation of a monotone response (y), using the standard isotonic-regression pool-adjacent-violators algorithm (PAVA). Core code from Raubertas (1994) with many modifications.
#'
#'
#'  Compute the isotonic regression of a numeric vector 'y', with
#'  weights 'wt', with respect to simple order. The core algorithm is still the one
#' coded by R.F. Raubertas, dated 02 Sep 1994. However, the input and output modules have been
#' modified to allow more flexible formats in either direction.
#' note that unlike centered-isotonic-regression (CIR, see \code{\link{cirPAVA}}), this algorithm does not use the dose (x) values at all. 
 
#  
##' @author C.R. Raubertas, Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}


#' @param y  can be either of the following: y values (response rates), a 2-column matrix with positive/negative response counts by dose, a \code{\link{DRtrace}} object or a \code{\link{doseResponse}} object. 
#' @param x dose levels (if not included in y). Note that the PAV algorithm doesn't really use them. 
#' @param wt weights (if not included in y).
#' @param full: logical, is a more complete output desired? if \code{FALSE} (default), only a vector of point estimates for y at the provided dose levels is returned
#' @param dec logical, is the true function is assumed to be monotone decreasing? Default \code{FALSE}.
#' @param ...	Other arguments passed on to the constructor functions that pre-process the input.

#' @return under default, returns a vector of y estimates at unique x values. With \code{full=TRUE}, returns a list of 3 \code{\link{doseResponse}} objects named \code{output,input,alg} for the output data at dose levels, the input data, and the function as fit at algorithm-generated points, respectively. For this function, the first and thrid objects are identical.

pava<-function (y,x=NULL,wt=rep(1,length(x)),full=FALSE,dec=FALSE,...) {

### converting to doseResponse object 
### Basically it's a numeric data frame with x,y,weight, and x increasing

dr=doseResponse(y,x,wt,...)
if (any(is.na(dr))) stop ("Missing values are not allowed.\n")  

m <- dim(dr)[1]
if (m <= 1) {  ## degenerate case: only one dose level
if (!full) return (dr$y)
else return(list(output=dr,input=dr,alg=dr))
}

dr0=dr ## clean copy of input data
### Decreasing monotone case: simple fix
if (dec) dr$y = -dr$y

#### Core algorithm

lvlsets <- (1:m)
 
repeat {
      viol <- (as.vector(diff(dr$y)) < 0)  # Find adjacent violators
      if (!(any(viol))) break
 
      i <- min( (1:(m-1))[viol])        # Pool first pair of violators
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i+1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      dr$y[ilvl] <- sum(dr$y[ilvl]*dr$weight[ilvl]) / sum(dr$weight[ilvl])
      lvlsets[ilvl] <- lvl1
}

if (dec) dr$y = -dr$y

if (!full) {
	return(dr$y)
} else {   ### somewhat redundant structure to match cirPAVA output
	return(list(output=dr,input=dr0,alg=dr))   }
}
