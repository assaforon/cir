
#' Utility 

bootversion<-function(dat,indices,fun,...) fun(dat[indices,],...)


plainCirBoot<-function(dat,fun,B=999,width=0.95,detailed=FALSE,...)
{
if(!is.DRtrace(dat)) stop("Input data has to be a DRtrace object.\n")
if(is.doseResponse(dat)) stop("Input data has to be a DRtrace object, but not a doseResponse object.\n")

n=dim(dat)[1]
indices=replicate(B,sample(n,replace=TRUE))
cat(dim(indices))
truest=fun(dat,...)
cat('\n',truest)
m=length(truest)

dout=matrix(NA,nrow=B+1,ncol=m)
dout[B+1,]=truest
for (a in 1:B) dout[a,]=fun(dat[indices[,a],],...)

return(dout)
}
