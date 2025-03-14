% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rowcalwm.R
\name{rowcalwm}
\alias{rowcalwm}
\title{Calculate Weighted Mean Date for a Radiocarbon Determination}
\usage{
rowcalwm(date, sigma, ...)
}
\arguments{
\item{date}{A radiocarbon determination in years BP.}

\item{sigma}{Sigma values corresponding to the radiocarbon determination.}

\item{...}{Additional arguments to be passed to the \code{rowcal} function.}
}
\value{
The weighted mean date for the radiocarbon determinations.
}
\description{
This function calculates the weighted mean date for a given radiocarbon determination, based on the rowcal function.
}
\details{
This function internally calls the \code{rowcal} function to compute some intermediate results.
It then calculates the weighted mean date based on these results.
}
\examples{
rowcalwm(5100, 35) # For 5100±35 BP

}
\seealso{
\code{\link{rowcal}} \code{\link{rowcalmedian}} \code{\link{rowcalmode}}
}
