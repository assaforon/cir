# Interesting run (#664) from a simulated up-and-down ensemble:
# (x will be auto-generated as dose levels 1:5)
dat=doseResponse(y=c(1/7,1/8,1/2,1/4,4/17), wt=c(7,24,20,12,17))
# The experiment's goal is to find the 30th percentile
quick1=quickIsotone(dat, adaptiveShrink = TRUE, target = 0.3)
# For inverse confidence intervals "the long way", 
#    we need a full CIR output object:
fwd1=cirPAVA(dat, full=TRUE, adaptiveShrink = TRUE, target = 0.3)
# Inverse intervals at doses
invDelta=deltaInverse(fwd1, target = NULL)

### Showing the data and the estimates
par(mar=c(3,3,4,1), mgp=c(2,.5,0), tcl=-0.25)
# Following command uses plot.doseResponse()
plot(dat, ylim=c(0.05,0.55), 
     las=1, xlim=c(0,6), main="Inverse-Estimation CIs") 

# The true response function; true target is where it crosses the y=0.3 line
lines(seq(0,7,0.1), pweibull(seq(0,7,0.1), shape=1.1615, scale=8.4839), col=4)
abline(h=0.3, col=2, lty=3) ### The experiment's official target

# Forward CIs; the "global" inverse interval just draws horizontal lines between them
# To get these "global" intervals calculated for you at specific targets, choose 'delta=FALSE' 
#      when calling quickInverse()
lines(quick1$lower90conf,lty=2,col=3) 
lines(quick1$upper90conf,lty=2,col=3) 
# Note that only the upper forward bounds intersect the horizontal line at y=0.3.
#   Therefore, via the "global" approach there won't be a finite CI for the target estimate.

# Now, the default "local" inverse interval, which is finite for the range of estimated y values.
# In particular, it is finite (albeit very wide) for y=0.3.
lines(invDelta[,1],quick1$y,lty=2, lwd=2)
lines(invDelta[,2],quick1$y,lty=2, lwd=2)

legend('topleft',pch=c(NA,'X',NA,NA),lty=c(1,NA,2,2), col=c(4,1,1,3), lwd=c(1,0,2,1),
       legend = c('True Curve', 'Observations', 'Local Interval (default)',
                  'Forward/Global Interval'), bty='n')

