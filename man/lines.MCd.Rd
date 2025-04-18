% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lines.MCd.R
\name{lines.MCd}
\alias{lines.MCd}
\title{Add Outline of Density Model to Current Plot}
\usage{
\method{lines}{MCd}(x, scalefactor = 1, ...)
}
\arguments{
\item{x}{Density model data, typically output from functions like MCdensity or phasedensity.}

\item{scalefactor}{Scaling factor for adjusting the density model.}

\item{...}{Additional graphical parameters to be passed to lines function.}
}
\value{
The density model outline is added to the current plot.
}
\description{
This function adds an outline of a density model to the current plot.
}
\examples{
# Plot density model and add outline
plot(MCdensity_output)
lines.MCd(MCdensity_output)

}
