
##' Returns centered-isotonic-regression y values at original x points 
##' @author Assaf P. Oron
#

# y: y values (responses). Can be a vector or a yes-no table (for binary responses)
#    it is Given as first argument, both for compatibility with 'pava'
# and to enable parallel running via 'apply' type routines
# x: treatments. Need to be pre-sorted in increasing order, with order matching y's
# wt: weights. Will be overwritten in case of a yes-no input for y

# full: if FALSE, only point estimates at x values are returned
# dec: Whether the true function is assumed to be
#    decreasing. Defaults to FALSE

cirPAVA <-function (y,x=NULL,wt=rep(1,length(x)),full=FALSE,dec=FALSE,...) {

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

repeat {

# Find adjacent violators
	viol <- (as.vector(diff(dr$y)) < 0)

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

if (dr$x[m]<max(dr0$x)) {
    dr<-rbind(dr,c(max(dr0$x),dr$y[m],0)) 
    m=m+1
}
if (dr$x[1]>min(dr0$x)) {
	dr<-rbind(c(min(dr0$x),dr$y[1],0),dr)  
	m=m+1
}
if (dec) dr$y = -dr$y

if (!full) {
	return(approx(dr$x,dr$y,dr0$x,rule=2)$y)
} else {
	
	dr1=dr0
	dr1$y=approx(dr$x,dr$y,dr1$x,rule=2)$y
	return(list(output=dr1,input=dr0,alg=dr))
}
}
