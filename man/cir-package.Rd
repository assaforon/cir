\name{cir-package}
\alias{cir-package}
\alias{cir}
\docType{package}
\title{
Isotonic Regression and Centered Isotonic Regression
}
\description{
Centered Isotonic Regression (CIR) is an empirical shrinkage improvement to isotonic regression (IR), requiring no tuning parameters. Numerical studies suggest around 2x reduction in mean-square error at design points (i.e., doses where observations were taken), and about 30 percent reduction elsewhere. 

This package implements CIR focusing on Binomial dose-response data, but compatible with other data types. Interval estimates are also supplied for both CIR and IR. Also included are classes and function for managing and plotting the trajectories or summaries of dose-response experiments.
}

\details{

The main functions to become acquainted with:

\tabular{ll}{
 \code{\link{quickIsotone}}\tab plug-and-play function providing CIR (or IR) point and interval estimates at all doses (x values) requested (default, all doses present in the input data).\cr
 \code{\link{quickInverse}}\tab the mirror-image function, providing point and interval estimates of all percentiles requested (i.e., dose-finding).\cr
 \code{\link{cirPAVA}}\tab more flexible forward point-estimate CIR function. There is also \code{\link{oldPAVA}} for traditional IR, arguably the most flexible input-format IR function currently available.\cr
  \code{\link{doseFind}}\tab likewise for inverse estimation.\cr
  \code{\link{isotInterval}}\tab interval estimate function. \cr
}


\tabular{ll}{
Package: \tab cir\cr
Type: \tab Package\cr
Version: \tab 1.1.0\cr
Date: \tab 2014-08-13\cr
License: \tab GPL 2\cr
}
}
\author{
Assaf Oron

Maintainer: Assaf Oron <assaf.oron@seattlechildrens.org>
}
\references{
}

\keyword{ package }
\seealso{

}
\examples{

}
