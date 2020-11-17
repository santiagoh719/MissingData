\name{TSLR}
\alias{TSLR_ranked}
\alias{TSLR}
\title{Theil-Sen missing value filled}
\description{
  Filled the ieth data with a Theil-Sen regression.
}
\usage{
  TSLR(Est_target, Est_reference,i, criterion='Max_Cor',cord=T,distans=F)

  TSLR_ranked(Est_target, Est_reference,i, criterion='Max_Cor',cord=T,distans=F)

}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{criterion}{criterion is which criteria to be used to select reference estation, Max_Cor (default) for the station with maximal correlation, or Nearest for the nearest station. If Nearest either cord or distans must be given}

  \item{cord}{matrix or data.frame, first column Longitude second column Latitude, first row coordinates of the target, then in order the coordinates of Est_reference}

  \item{distans}{numeric vector,  in order the distans from target to Est_reference}

}

\value{
  \code{TSLR},returns ieth data with a Theil-Sen regression of the data.

  \code{TSLR_ranked} returns the ieth data with a Theil-Sen regression of the ranked data.

}
\references{
- Theil, H., 1992. A rank-invariant method of linear and polynomial regression analysis. In: Henri Theil’s contributions to economics and econometrics. Springer, pp. 345–381. URL https://doi.org/10.1007/978-94-011-2546-8_20

- Sen, P. K., 1968. Estimates of the regression coefficient based on kendall’s tau. Journal of the American statistical association 63 (324), 1379–1389.
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
    TSLR_ranked(Est_target, Est_reference,i)
  )
  print(
      TSLR(Est_target, Est_reference,i, criterion='Max_Cor')
  )
}
