#' Bootstrap confidence intervals for isotonic regressions
#'
#'
#' Performs percentile bootstrap interval estimation, simultaneously for all estimated y values (forward) or percentiles (inverse).

#' @param dat The input data. Has to be a \code{\link{DRtrace}} or \code{\link{doseResponse}} object.
#' @param fun The estimation function
#' @param xout the vector of x values on which prediction is needed. If NULL (default) will be set as the sorted unique values of the input x values.
#' @param B The number of bootstrap replications (default 999). The true data are added to the denominator as well.
#' @param width the confidence-interval witdh in probability units.
#' @param stratify logical: should the bootstrapping be stratified by x value to maintain the original design? Defaults to FALSE. See 'Details' for more.
#' @param detailed logical: should the output be detailed, or just the CI boundaries? (default FALSE).
#' @param verbose logical: should start and end times be echoed out to assess runtime? Default TRUE.
#' @param ... parameters passed on to 'fun'.

#' @return a matrix with B rows, and as length(xout) columns, containing the full


#' @seealso \code{\link{pava}},\code{\link{cirPAVA}}
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}


plainCirBoot<-function(dat,fun=cirPAVA,xout=NULL,B=999,width=0.95,stratify=FALSE,detailed=FALSE,verbose=TRUE,...)
{
if(!is.DRtrace(dat)) stop("Input data has to be a DRtrace or doseResponse object.\n")
if(stratify) {
	dat=doseResponse(dat)
	if(verbose)	cat("Stratified (Binomial) bootstrap...\n")
} else if(is.doseResponse(dat)) stop("Non-stratified only possible with DRtrace input.\n")
if(width<=2/B || width>=1-2/B) stop("Unrealistic input interval width.\n")

if(verbose) cat(date(),' ')

n=dim(dat)[1]
if(!stratify) {
	indices=replicate(B,sample(n,replace=TRUE))
} else counts=apply(dat,1,function(x,n) rbinom(n,size=x[3],prob=x[2]),n=B)

if(is.null(xout)) xout=sort(unique(dat$x)) #x values for prediction
truest=fun(dat,outx=xout,...)
m=length(truest)

dout=matrix(NA,nrow=B+1,ncol=m)
colnames(dout)=xout
dout[B+1,]=truest
for (a in 1:B) {

	if(!stratify) {
		tmpdat=dat[indices[,a],]
	} else {
		tmpdat=dat
		tmpdat$y=counts[a,]/dat$weight
	}
	dout[a,]=fun(tmpdat,outx=xout,...)
}
if(verbose) cat(date(),'\n')

tails=(1-width)/2
ciout=apply(dout,2,quantile,prob=c(tails,1-tails),type=6)
if(!detailed) return(ciout)

return(list(interval=ciout,samples=dout))
}
