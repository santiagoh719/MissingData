\name{One value replace}
\alias{AA}
\alias{AM}
\alias{NN}
\alias{SBE}
\alias{Clim_AA}

\title{Replace with one specific value}
\description{
  Replace the ieth data with one specific value with simple methods.
}
\usage{
  AA(Est_reference, i)

  AM(Est_reference, i)

  NN(Est_reference, i, cord = T , distans = F)

  SBE(Est_target, Est_reference, i)

  Clim_AA(Est_target, i)
}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{cord}{matrix or data.frame, first column Longitude second column Latitude, first row coordinates of the target, then in order the coordinates of Est_reference}

  \item{distans}{numeric vector,  in order the distans from target to Est_reference}

}

\value{
  \code{AA}, \code{AM}, \code{NN} and \code{SBE} give the ieth value of target.
}
\details{
  \code{AA} compute the ieth value of target with the aritmetic mean of the ieth value of Est_reference.
  \code{AM} compute the ieth value of target with the median of the ieth value of Est_reference.
  \code{NN} compute the ieth value of target with the nearest station of Est_reference.
  \code{SBE} compute the ieth value of target with the station of Est_reference that show hieghest correlation with target.
  \code{Clim_AA} compute the temporal mean of the target station. The data must be monthly value order by month with no time gap.
}

\references{
 - Sattari, M. T., Rezazadeh-Joudi, A., & Kusiak, A. (2017). Assessment of different methods for estimation of missing data in precipitation studies. Hydrology Research, 48(4), 1032-1044. https://doi.org/10.2166/nh.2016.364

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
    AA(Est_reference, i)
  )
  print(
    AM(Est_reference, i)
  )
  print(
    NN(Est_reference, i, cord)
  )
  print(
    SBE(Est_target, Est_reference, i)
  )
  print(
    Clim_AA(Est_target, i)
  )
}
