% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tides.R
\name{IT}
\alias{IT}
\title{Inundation time}
\usage{
IT(h, h0, hoffset = 0, dtMax = 15, unit = "mins")
}
\arguments{
\item{h}{Water level time series. data frame with time and h column}

\item{h0}{Reference level, either single valued or vector with same length as h}

\item{hoffset}{Offset level to cope with small fluctuations due to rain, ripples. h <= h0 + hoffset is considered dry; h> h0+hoffset is considered wet}

\item{dtMax}{Maximum time interval in continuous water level series. Larger time intervals are considered gaps}

\item{unit}{Unit of dtMax.}
}
\value{
a list containing:
\itemize{
   \item{IT }{Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of inundation}
   \item{DT }{Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of dry time}
}
}
\description{
Calculate inundation times, i.e. time intervals for which water level h > h0. Care must be taken when there are gaps (long time periods for which there is no data )in the time series. Either the erroneous values have to removed manually, or a wrapper making use of the function gapsts can be used.
}
\author{
Tom Cox <tom.cox@uantwerp.be>
}
\keyword{utilities}
