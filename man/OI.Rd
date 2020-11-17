\name{Optimal Interpolation}
\alias{OI}
\title{Optimal Interpolation}
\description{
  Filled the ieth data with Optimal Interpolation.
}
\usage{
  OI(Est_target, Est_reference, i, cord= T, criterion = 'SBE' )

}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{criterion}{criterion is which criteria to be used to select the fisrt guest value, SBE (default) for the station with maximal correlation, or NN for the nearest station. If Nearest cord must be given}

  \item{cord}{matrix or data.frame, first column Longitude second column Latitude, first row coordinates of the target, then in order the coordinates of Est_reference}

}

\value{
  \code{OI},returns ieth data interpolate with the Optimal Interpolation method discribe in E ISCHEID 1999.

}
\references{
 - Eischeid, J. K., Pasteris, P. A., Diaz, H. F., Plantico, M. S., & Lott, N. J. (2000). Creating a serially complete, national daily time series of temperature and precipitation for the western United States. Journal of Applied Meteorology, 39(9), 1580-1591. https://doi.org/10.1175/1520-0450(2000)039<1580:CASCND>2.0.CO;2

}

\examples{
  # Make a Target serie with a seasonal sicle, and a missing value
  set.seed(524)
  Est_target <- sin(1:100 * pi/6) * 50 + 70 + rnorm(100, 0, 10)
  i <- 30
  Est_target[i] <- NA

  # Make reference station data with 7 stations:
  Est_reference <- matrix(
  data = c(sin(1:100 * pi/6) * 100 + 300 + rnorm(100, 0, 50),
  sin(1:100 * pi/6 +1.5) * 15 + 35 + rnorm(100, 0, 9),
  sin(1:100 * pi/6 + 2) * 8 + 10 + rnorm(100, 0, 1),
  sin(1:100 * pi/6 + 1) * 30 + 200 + rnorm(100, 0, 10),
  sin(1:100 * pi/6 - 1) * 80 + 120 + rnorm(100, 0, 15),
  sin(1:100 * pi/6 + 0.5) * 40 + 67 + rnorm(100, 0, 10),
  sin(1:100 * pi/6 -0.5) * 10 + 14 + rnorm(100, 0, 5)),
  ncol = 7, nrow = 100)
  Est_reference[Est_reference < 0] <- 0

  cord <- matrix( c(rnorm(8,-27,1),rnorm(8,-61,1)),ncol = 2, nrow = 8)

  # Estimate the ieth value
  print(
    OI(Est_target, Est_reference, i)
  )
}
