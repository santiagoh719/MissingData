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

