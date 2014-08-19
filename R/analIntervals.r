##' Returns analytic interval estimates, given isotonic-regression (centered or not) point estimates
#'
#'
#' Forward interval estimation of a monotone response (y) as a function of dose (x), given a set of point estimates produced by \code{\link{pava}} or \code{\link{cirPAVA}}.
#'
#'
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}

isotInterval<-function(isotPoint,xout,conf=0.9)