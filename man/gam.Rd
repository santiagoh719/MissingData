\name{gam_missing}
\alias{gam_missing}
\title{GAM missing value filled}
\description{
  \code{gam_missing} filled the ieth data with a GAM model.
}
\usage{
  gam_missing(Est_target, Est_reference,i)

}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}
}

\value{
  \code{gam_missing}, returns the ieth value of target with an interpolation using a GAM model.
}
\details{
  \code{gam_missing} filled the ieth data with a GAM model.
}
\references{
  - Hastie, T., & Tibshirani, R. (1987). Generalized additive models: some applications. Journal of the American Statistical Association, 82(398), 371-386. https://www.tandfonline.com/doi/abs/10.1080/01621459.1987.10478440
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
    gam_missing(Est_target, Est_reference, i)
  )

}
