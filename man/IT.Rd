\name{IT}
\alias{IT}
\title{Inundation time}
\description{Calculate inundation times, i.e. time intervals for which water level h > h0. Care must be taken when there are gaps (long time periods for which there is no data )in the time series. Either the erroneous values have to removed manually, or a wrapper making use of the function gapsts can be used.}
\usage{IT(h, h0, h0marg = 0.3, dtMax, unit = "mins")}
\arguments{
  \item{h }{Water level time series. data frame with time and h column}
  \item{h0 }{Reference level, either single valued or vector with same length as h}
  \item{h0marg }{Margin on reference level, to cope with small fluctuations in the Water level time series}
  \item{dtMax }{Maximum time interval in continuous water level series}
  \item{unit }{Unit of dtMax}
}
\value{
  a list containing:
  \item{IT }{Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of inundation}
  \item{DT }{Data frame with start time (t1), end time (t2) and duration (dt, unit = unit) of dry time}
}
\author{Tom Cox <tom.cox@ua.ac.be>}
\keyword{utilities}

