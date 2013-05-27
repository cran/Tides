\name{IF}
\alias{IF}
\title{Inundation frequence}
\description{Calculate inundation frequence at height h0 from series of high water levels H}
\usage{IF(H, h0,N = length(H[,1]))}
\arguments{
  \item{H }{High water levels. Data Frame with column h}
  \item{h0 }{Reference level for which IF has to be calculated}
  \item{N }{number of cycles in time series, equals the number of high water levels when these are complete (= default value)}
}
\value{
  Inundation frequence [\%] at reference level h0
}
\author{Tom Cox <tom.cox@ua.ac.be>}
\keyword{utilities}

