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

