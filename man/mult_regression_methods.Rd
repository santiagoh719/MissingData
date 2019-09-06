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

