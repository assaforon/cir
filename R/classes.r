#' Functions to check class

#' @param dr the (purported) dose-response object being checked
#' @return logical value indicating whether the object meets class definition

is.DRtrace<-function(dr)
{
if(!inherits(dr,"DRtrace")) return(FALSE)
if(!inherits(dr,"data.frame")) return(FALSE)
if(!all(c("x","y","weight") %in% names(dr))) return(FALSE)
if(!all(sapply(dr[,c("x","y","weight")],is.numeric))) return(FALSE)
if(any(dr$weight<0)) return(FALSE)
return(TRUE)
}

is.doseResponse<-function(dr)
{
if(!is.DRtrace(dr)) return(FALSE)
if(!inherits(dr,"doseResponse")) return(FALSE)
if(any(duplicated(dr$x))) return(FALSE)
if (any(dr$x!=sort(dr$x))) return(FALSE)
return(TRUE)
}

#########################

#' Constructor functions for DRtrace and doseResponse classes

DRtrace<-function(y,x=NULL,wt=NULL,noyes=FALSE)
{
if(is.DRtrace(y)) return(y)

ll=dim(y)

if (length(ll)>2 || (length(ll)==2 && ll[2]!=2)) stop ("For conversion to dose-response, y can only be a vector or yes-no table.\n")

if (length(ll)==2 && ll[2]==2) { # converting a yes-no table

	xvals<-x
	if(is.null(xvals)){
		warning("No dose (x) values provided; using natural numbers.\n")
		xvals=1:ll[1]
	}
	if(length(xvals) != ll[1]) stop("Mismatched lengths.\n")
	
	yntab=y
	if(noyes) yntab=yntab[,2:1]  # if the table is no-yes rather than yes-no, we reverse it here.
    y<-c(rep(1,sum(yntab[,1])),rep(0,sum(yntab[,2])))
    x<-c(rep(xvals,yntab[,1]),rep(xvals,yntab[,2]))
    wt<-rep(1,length(y))  ## in case of yes-no table we ignore incoming weights 
	warning("Raw data is a yes-no table; therefore, observation order is arbitrary.\n")
}

n=length(y)
if(is.null(wt)) wt=rep(1,n)
if(length(x)!=n || length(wt)!=n) stop("Mismatched lengths.\n")

tout<-data.frame(x=x,y=y,weight=wt)
attr(tout,'class')<-c('DRtrace','data.frame')
return(tout)
}


#############

doseResponse<-function(y,x=NULL,wt=NULL,...)
{
if(is.doseResponse(y)) return(y)

z<-suppressWarnings(DRtrace(y,x,wt,...))

tout<-data.frame(x=sort(unique(z$x)),y=tapply(z$y,z$x,mean),weight=tapply(z$weight,z$x,sum))
attr(tout,'class')<-c('doseResponse','DRtrace','data.frame')
return(tout)
}








