neundatLevel=doseResponse(x=c(1,2.5,5,10,20,25),y=c(rep(0,4),2/9,1),wt=c(3,4,5,4,9,2))

neundatDose=doseResponse(x=c(1,2.5,5,10,20,25),y=c(rep(0,4),2/9,1),wt=c(3,4,5,4,9,2))
plot(neundat,main="Neuenschwander et al. (2008) Final Dose-Toxicity",ylim=c(0,1),xlab="Dose (mg/sq.m./wk)",ylab="Toxicity Response Curve (F)",cex.main=1.5)

neunIR=oldPAVA(neundat)
neunCIR0=cirPAVA(neundat,full=TRUE,strict=TRUE)
lines(neundat$x,neunIR,col=4,type='b',pch=19)
lines(neunCIR0$alg$x,neunCIR0$alg$y,type='b',col=2,lwd=2)
