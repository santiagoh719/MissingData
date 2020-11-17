\name{Weigth methods}
\alias{UK}
\alias{NR}
\alias{NR_1952}
\alias{CWM}

\title{Replace using weigth means}
\description{
  Replace the ieth data using weigth means or single estimator weigthed.
}
\usage{
  UK(Est_target, Est_reference, i, W = NULL)

  NR(Est_target, Est_reference, i)

  NR_1952(Est_target, Est_reference, i, tt='month', W = NULL)

  CWM(Est_target, Est_reference, i)

}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{tt}{tt is periodicity of sample data, acep 'month' (default), 'day', 'hour' and '6hour'}

  \item{W}{weight to used in the method, see details}

}

\value{
  \code{UK}, \code{NR}, \code{NR_1952} and \code{CWM} give the ieth value of target.
}
\details{
  \code{UK} compute the ieth value of target with United Kingdom traditional method. Select the station with hieghtes correlation with target, and multiple it by the reason of the mean values of target and it self.
  \code{NR_1952} compute the ieth value of target with the Normal Ratio method. Compute the mean of the ieth values of Est_reference weight by the ratio of the acumalated annual mean of the reference and target.
  \code{NR} compute the ieth value of target with the Normal Ratio method modified by Young 1992. Compute the mean of the ieth values of Est_reference weight by the  t-student estimator of correlation.
  \code{SBE} compute the mean of the ieth values of Est_reference weight by the correlation.
  \code{W} in \code{UK} method is the mean precipitation of the target station divided the mean precipitation of the reference station, meanwhile in \code{NR_1952} is similar but the mean is of the total annual precipitation. In both cases it is not mandatory to give this value, but if the imputs are anomalies it is recomended.
}

\references{
 - Sattari, M. T., Rezazadeh-Joudi, A., & Kusiak, A. (2017). Assessment of different methods for estimation of missing data in precipitation studies. Hydrology Research, 48(4), 1032-1044. https://doi.org/10.2166/nh.2016.364

  - Teegavarapu, R. S., & Chandramouli, V. (2005). Improved weighting methods, deterministic and stochastic data-driven models for estimation of missing precipitation records. Journal of hydrology, 312(1-4), 191-206. https://doi.org/10.1016/j.jhydrol.2005.02.015

 - Young, K. C. (1992). A three-way model for interpolating for monthly precipitation values. Monthly Weather Review, 120(11), 2561-2569. https://doi.org/10.1175/1520-0493(1992)120<2561:ATWMFI>2.0.CO;2

 - Paulhus, J. L., & Kohler, M. A. (1952). Interpolation of missing precipitation records. Monthly Weather Review, 80(8), 129-133. http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.394.6283&rep=rep1&type=pdf
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
    NR_1952(Est_target, Est_reference, i)
  )
  print(
    NR(Est_target, Est_reference, i)
  )
  print(
    CWM(Est_target, Est_reference, i)
  )
  print(
    UK(Est_target, Est_reference, i)
  )
}
