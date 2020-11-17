\name{MLR}
\alias{MLR_ranked}
\alias{MLR}
\alias{RMLR_ranked}
\alias{RMLR}
\title{Multiple Linear Regression}
\description{
  Filled the ieth data with Multiple Linear Regression.
}
\usage{
  MLR(Est_target, Est_reference,i, criterion='AIC',alpha=0.05)

  MLR_ranked(Est_target, Est_reference,i, criterion='AIC',alpha=0.05)

  RMLR(Est_target, Est_reference,i, M_Method= 'MM')

  RMLR_ranked(Est_target, Est_reference,i, M_Method= 'MM')

}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{criterion}{criterion is which criteria to be used to select which stations to use from Est_reference, possible:  AIC (default), BIC and Cp (mallows Cp)}

  \item{alpha}{significanse level to select the model}

  \item{M_Method}{currently either M-estimation or MM-estimation or (for the formula method only) find the model frame. MM-estimation is M-estimation with Tukey's biweight initialized by a specific S-estimator. }

}

\value{
  \code{MLR},returns ieth data with a multilinear regression of the data.

  \code{MLR_ranked},returns ieth data with a multilinear regression of the ranked data.

  \code{RMLR},returns ieth data with a multilinear robust regression of the data.

  \code{RMLR_ranked},returns ieth data with a multilinear robust regression of the ranked data.

}
\details{
\code{RMLR} and \code{RMLR_ranked} use the SignifReg function from package SignifReg.
\code{MLR} and \code{MLR_ranked} use the rlm function from package MASS.

}

\references{
 - Simolo, C., Brunetti, M., Maugeri, M., & Nanni, T. (2010). Improving estimation of missing values in daily precipitation series by a probability density function‐preserving approach. International Journal of Climatology, 30(10), 1564-1576.  https://doi.org/10.1002/joc.1992
 - Presti, R. L., Barca, E., & Passarella, G. (2010). A methodology for treating missing data applied to daily rainfall data in the Candelaro River Basin (Italy). Environmental monitoring and assessment, 160(1-4), 1. https://doi.org/10.1007/s10661-008-0653-3
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

  # Estimate the ieth value
  print(
    RMLR_ranked(Est_target, Est_reference,i)
  )
  print(
    RMLR(Est_target, Est_reference,i)
  )

}
