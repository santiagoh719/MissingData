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

