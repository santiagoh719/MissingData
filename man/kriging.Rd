\name{Kriging}
\alias{Ord_kriging}
\alias{KED_Alt}
\alias{KED}
\title{Kriging missing value filled}
\description{
  Fill the ieth data with kriging.
}
\usage{

  Ord_kriging( Est_reference, i, cord, method = NULL)

  KED_Alt( Est_reference, i, cord, method = NULL)

  KED(Est_reference,Est_reference_var, var_target, i, cord, method = NULL)
}
\arguments{
  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

 \item{cord}{matrix or data.frame, first column Longitude second column Latitude, for KED_Alt a third column with altitude is necesary, first row coordinates of the target, then in order the coordinates of Est_reference}

  \item{method}{character, model to fit the variogram, possibles: "Exp", "Mat", "Gau","Sph". Is not given will try with all four and select the better one.}

  \item{Est_reference_var}{vector of another variable in the stations to use for the kriging model}

  \item{var_target}{numeric, value of the other variable used in the target}

}

\value{
  \code{Ord_kriging}, \code{KED_Alt},\code{KED}, returns the ieth value of target with a kriging interpolation.
}
\details{
  \code{Ord_kriging} compute ordinary kriging.
  \code{KED_Alt} and \code{KED} compute kriging with external drift. \code{KED_Alt}  compute kriging with external drift using as secondary variable the altitude.
}

\references{
  - Erxleben, J., Elder, K., & Davis, R. (2002). Comparison of spatial interpolation methods for estimating snow distribution in the Colorado Rocky Mountains. Hydrological Processes, 16(18), 3627-3649. https://doi.org/10.1002/hyp.1239
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

  cord2 <- cbind(cord,rnorm(8,1200,400))
  # Estimate the ieth value
  print(
    Ord_kriging( Est_reference, i, cord, method = 'Gau')
  )
  print(
    KED_Alt( Est_reference, i, cord2, method = 'Gau')
  )
}
