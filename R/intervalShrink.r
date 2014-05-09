#' Shrinkage bootstrap confidence intervals for isotonic regressions
#'
#'
#' Performs shrinkage percentile bootstrap interval estimation, simultaneously for all estimated y values (forward) or percentiles (inverse). Bootstrap samples are size-matched to the real data (i.e., stratified by value of x). When generating these samples, we use the isotonic point estimates shrunk towards a generic or user-defined prior, to account for e.g., "variability sources outside the model".

#' @param dat The input data. Has to be a \code{\link{DRtrace}} or \code{\link{doseResponse}} object.
#' @param fun The estimation function
#' @param xout the vector of x values on which prediction is needed. If NULL (default) will be set as the sorted unique values of the input x values.
#' @param B The number of bootstrap replications (default 999). The true data are added to the denominator as well.
#' @param width the confidence-interval width in probability units. Default 0.9.
#' @param prior either NULL, or a \code{\link{doseResponse}} object with the shrinkage prior.
#' @param detailed logical: should the output be detailed, or just the CI boundaries? (default FALSE).
#' @param verbose logical: should start and end times be echoed out to assess runtime? Default TRUE.
#' @param ... parameters passed on to 'fun'.

#' @return a matrix with B rows, and as length(xout) columns, containing the full


#' @seealso \code{\link{pava}},\code{\link{cirPAVA}}
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}


shrinkCirBoot<-function(dat,fun=cirPAVA,xout=NULL,B=999,width=0.9,prior=NULL,detailed=FALSE,verbose=TRUE,...)
{
if(!is.DRtrace(dat)) stop("Input data has to be a DRtrace or doseResponse object.\n")
if(!is.doseResponse(dat))	dat=doseResponse(dat)
if(width<=2/B || width>=1-2/B) stop("Unrealistic input interval width.\n")

if(verbose) cat(date(),' ')

basest=fun(dat,...)
m=length(basest)
#print(basest)
if(is.null(prior)) prior=doseResponse(y=(1:m)/(m+1),x=xout,wt=rep(1/m,m))

# validation on user-provided prior
if(!is.doseResponse(prior)) stop("User-defined prior must be a doseResponse object.\n")
if(!all(dat$x %in% prior$x)) stop("User-defined prior must include all requested x values.\n")

shrunkest=dat
shrunkest$y=(dat$weight*basest+prior$weight*prior$y)/(dat$weight+prior$weight)
#print(shrunkest)

### Bootstrap data generated here
counts=apply(shrunkest,1,function(x,n) rbinom(n,size=x[3],prob=x[2]),n=B)

if(is.null(xout)) xout=sort(unique(dat$x)) #x values for prediction
k=length(xout)

dout=matrix(NA,nrow=B+1,ncol=k)
colnames(dout)=xout
### True estimate is last row
dout[B+1,]=approx(dat$x,basest,xout=xout,rule=2)$y
### The bootstrap loop is for estimation
for (a in 1:B) {

	tmpdat=dat
	tmpdat$y=counts[a,]/dat$weight
	dout[a,]=fun(tmpdat,outx=xout,...)
}
if(verbose) cat(date(),'\n')

tails=(1-width)/2
ciout=apply(dout,2,quantile,prob=c(tails,1-tails),type=6)
if(!detailed) return(ciout)

return(list(interval=ciout,samples=dout))
}
