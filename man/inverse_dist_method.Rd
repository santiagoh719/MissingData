\name{IDW}
\alias{IDW}
\alias{IDW_simple}
\alias{RNNWM}
\alias{ADW}

\title{Inverse distance weigth interpolation methods}
\description{
  Fill the ieth data with an inverse distance weigth method.
}
\usage{
  IDW(Est_reference, i, k, cord = T , distans = F)

  IDW_simple(Est_reference, i, cord = T , distans = F)

  RNNWM(Est_target, Est_reference, i, k, reg= 1)

  ADW(Est_reference, i, cord, CDD, m = 4)
}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{k}{numeric, exponent of the inverse distance weigth}

  \item{cord}{matrix or data.frame, first column Longitude second column Latitude, first row coordinates of the target, then in order the coordinates of Est_reference}

  \item{distans}{numeric vector,  in order the distans from target to Est_reference}

  \item{reg}{integer value, large of the paterns that will be look for, default 1 acept 1 or 2}

  \item{m}{numeric, exponent of the inverse distance weigth}

  \item{CDD}{numeric, search radius of reference stations (correlation decay distance)}

}

\value{
  \code{IDW}, returns the ieth value of target with an interpolation using a inverse distance weigth model
  \code{IDW_simple}, returns the ieth value of target with an interpolation using a inverse distance weigth model with k=1
  \code{RNNWM}, returns the ieth value of target with an interpolation using a inverse distance weigth model, with a different mesearu of distance.
}
\details{
    \code{IDW} compute the ieth value of target with an interpolation using a inverse distance weigth model. \code{IDW_simple} is as \code{IDW} with k=1.
  \code{RNNWM}, compute the ieth value of target with an interpolation using a inverse distance weigth model but with adifferent mesearu of distance.It use a surrogate measure for distance base on Teegavarapu and Chandramouli 2005.
}
\references{
  - New, M., Hulme, M., & Jones, P. (2000). Representing twentieth-century space–time climate variability. Part II: Development of 1901–96 monthly grids of terrestrial surface climate. Journal of climate, 13(13), 2217-2238. https://doi.org/10.1175/1520-0442(2000)013<2217:RTCSTC>2.0.CO;2

  - Teegavarapu, R. S., & Chandramouli, V. (2005). Improved weighting methods, deterministic and stochastic data-driven models for estimation of missing precipitation records. Journal of hydrology, 312(1-4), 191-206. https://doi.org/10.1016/j.jhydrol.2005.02.015
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
      IDW(Est_reference, i, k = 3, cord)
      )
print(
  IDW_simple(Est_reference, i, cord)
)
print(
  RNNWM(Est_target, Est_reference, i, k = 2, reg= 1)

)
print(
  ADW(Est_reference, i, cord, CDD = 300, m = 4)
)

}
