% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rowcal_datelist.R
\name{rowcal_datelist}
\alias{rowcal_datelist}
\title{Calibrate a List of Radiocarbon Dates}
\usage{
rowcal_datelist(dl = CLIP(), prob = 0.95, copy = TRUE, ...)
}
\arguments{
\item{dl}{A data frame containing radiocarbon determination dates. If the data frame has two columns, it is assumed that the first column represents the radiocarbon dates, and the second column represents their corresponding uncertainties. If it has three columns, the third column represents the calibration curve (default is 'intcal'). If nothing is provided the function will attempt to read from the clipboard.}

\item{prob}{The probability range for the calibrated dates. Default is 0.95.}

\item{copy}{Logical. If TRUE, the output is copied to the clipboard (Windows and macOS). Default is TRUE.}

\item{...}{Additional arguments to be passed to the 'rowcal' function.}
}
\value{
A string with the two-sigma range of calibrated dates on each line.
}
\description{
This function calibrates a list of radiocarbon determination dates and returns text indicating the two-sigma range of calibrated dates.
}
\details{
The format includes the calibrated date range with respect to the common era (BC/AD).
}
\examples{
# Calibrate a list of radiocarbon dates with default parameters
dl<-BIRE[BIRE$Where=='Britain' & BIRE$ccode=='F',2:3] # British burnt mounds
rowcal_datelist(dl, copy=FALSE)

# Calibrate a list of radiocarbon dates with custom parameters
rowcal_datelist(dl, prob=0.68, copy=FALSE)

}
