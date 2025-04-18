% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCsimulate.R
\name{MCsimulate}
\alias{MCsimulate}
\title{Simulate Radiocarbon Measurements from Distribution in Calendar Time}
\usage{
MCsimulate(mod, Nd = 1000, type = "KDE", ...)
}
\arguments{
\item{mod}{Model object representing the distribution in calendar time; this must be a data frame or list with components \code{x} and \code{y}}

\item{Nd}{Number of simulated radiocarbon measurements to generate. Defaults to 1000.}

\item{type}{Type of simulation to perform. Options are 'KDE' for Kernel Density Estimation or 'SPD' for Summed Probability Distribution. Defaults to 'KDE'.}

\item{...}{Additional parameters to be passed to the underlying KDE or SPD function.}
}
\value{
A simulated dataset of radiocarbon measurements, summarised using either KDE or SPD.
}
\description{
This function simulates radiocarbon measurements from a given distribution in calendar time.
}
\examples{
# Simulate 100 radiocarbon measurements following an expotential growth rate
# of 0.1, in calendar time from 2000 BC to AD 500:
sim_data <- MCsimulate(list(x=seq(-2000,500,10), y=exp((0:250)*0.1)), Nd=100)
plot(sim_data)


}
