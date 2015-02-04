\name{extrema}
\alias{extrema}
\title{Extrema}
\description{Calculate maxima and minima of quasi periodic time series. The function works good when the time series is 'smooth enough'. The function will fail when multiple local extrema exist, with a magnitude similar to the extremum extremorum and in an interval larger than [t-T2, t+T2] around the extremum extremorum at time t. When this is the case, this might be solved by adding extra wrappers around the function.}
\usage{extrema(h, h0, T2 = 5*60*60, hoffset = 0, filtconst = 1)}
\arguments{
  \item{h }{Water level time series. Data frame with time and h column}
  \item{h0 }{Reference level, either single valued vector with dimension corresponding to h}
  \item{T2 }{'Lower' bound on half the quasi period, but higher than expected stagnant phase; default = 5h}
  \item{hoffset}{Offset level, to prevent spurious maxima generation due to small fluctuations}
  \item{filtconst }{Filtering constant for smoothing the time series with \code{filter(x, rep(1/filtconst,filtconst))}}
}
\value{
  a list containing:
  \item{HL }{Data frame with time, water level and reference level of the extrema. Column HL denotes H(igh) (=maximum) or L(ow) (=minimum) water levels}
  \item{h }{Original water level data frame with additional columns HL, denoting H(igh) or L(ow) water phase and N(umber) of tidal cycle.}
}
\author{Tom Cox <tom.cox@untwerpbe>}
\keyword{utilities}

