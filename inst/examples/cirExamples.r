# Interesting run from a simulated up-and-down ensemble:
# (x will be auto-generated as 1:5)
dat=doseResponse(y=c(6/25,4/17,11/17,1/3,1),wt=c(25,34,17,3,1))
# CIR using the default 'quick' function that also provides 90% CIs
quick1=quickIsotone(dat)
quick1
# Use 'estfun' argument to operate the same function with old PAVA as the estimator
quick0=quickIsotone(dat,estfun=oldPAVA)
quick0

