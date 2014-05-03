
# Generic dose
# Assaf Oron, 10/2007
# Returns Point & Interval Estimates of Target percentile

### ARGUMENTS:
# y: Yes-no table of binary responses. Can contain rows of zeros
# xseq: Treatment values matched to the responses
# target: The target response rate (between 0 and 1)
# xbounds,ybounds: used for interpolation in case (estimated)
# target falls xoutside CIR xoutput boundaries.
# full: complete xoutput or estimates only
# cioption: which method to use for interval estimation
# can choose from 'poisson', 't' or binomial (any other string)
# plist: percentile list for interval estimation. Is a vector
# of any length

invCIR<-function(y,x=NULL,wt=rep(1,length(x)),target,...) doseFind(y,x,wt,estfun=cirPAVA,target=target,full=FALSE,...) 

doseFind<-function(y,x=NULL,wt=rep(1,length(x)),estfun=cirPAVA,target,full=FALSE,dec=FALSE,...) {


### converting to doseResponse object 
### Basically it's a numeric data frame with x,y,weight, and x increasing

dr=doseResponse(y,x,wt,...)
if (any(is.na(dr))) stop ("Missing values are not allowed.\n")  


### Point Estimate  ###########################

# We start via forward estimation, just calling 'cir.pava'

pavout<-estfun(y=dr,full=TRUE,dec=dec)


newx<-pavout$alg$x
newy<-pavout$alg$y
newn=pavout$alg.wt

### Error control if the interpolation target is xoutside boundary
if (min(newy)>target || max(newy)<target) return(NA)
### Otherwise, the estimate is generated by using 'approx' with x and y interchanged

xout=approx(x=newy,y=newx,xout=target,ties="ordered",rule=2)$y

#### xoutput

if (!full)  return(xout) 

return (list(xout=xout,input=dr,cir=pavout$alg))
}
