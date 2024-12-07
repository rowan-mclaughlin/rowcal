% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_density.R
\name{bootstrap_density}
\alias{bootstrap_density}
\title{Bootstrap a confidence interval for a set of Calendar Dates}
\usage{
bootstrap_density(DATA = CLIP()[, 1], N = 100, bw = 30)
}
\arguments{
\item{DATA}{A vector of calendar dates. BY default the function takes input from the clipboard.}

\item{N}{Number of bootstrap samples.}

\item{bw}{Bandwidth parameter for density calculation.}
}
\value{
A matrix of type 'MCd' containing the \code{MCdensity}-like object.
}
\description{
This function generates a \code{MCdensity}-like object from a list of calendar dates without uncertainty using bootstrapping (sampling with replacement).
}
\examples{
# Simulate 50 calendar dates centered on AD 1066 with a sigma of 50
dates<-rnorm(50, 1066, 50)
# calculate and plot bootstrap density model
denmod<-bootstrap_density(dates)
plot(denmod)

}
\seealso{
\code{\link{MCdensity}} \link{'mixdensity'} \link{'phasedensity'} \code{\link{density}} \code{\link{plot.MCd}}
}
