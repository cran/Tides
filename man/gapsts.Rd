\name{gapsts}
\alias{gapsts}
\title{Gaps in time series}
\description{Determine gaps (time intervals bigger than a certain maximum dtMax) in a time series}
\usage{gapsts(ts, dtMax,unit = "mins")}
\arguments{
  \item{ts }{Array of times, not necessarily of class \code{POSIXt}. Time sequences with differences of subsequent values < dtMax are considered a continuous series. Time differences > dtMax are considered gaps between two consecutive series.}
  \item{dtMax }{Maximum time interval that is not considered as a gap.}
  \item{unit }{Unit of dtMax. only used when ts is of class \code{POSIXt}}
}
\value{
  Data frame with the initial time, end time and time difference (unit = unit) of each interval > dtMax
}
\author{Tom Cox <tom.cox@uantwerp.be>}
\keyword{utilities}

