\name{One value replace}
\alias{AA}
\alias{AM}
\alias{NN}
\alias{SBE}

\title{Replace with one specific value}
\description{
  Replace the ieth data with one specific value with simple methods.
}
\usage{
  AA(Est_reference, i)

  AM(Est_reference, i)

  NN(Est_reference, i, cord = T , distans = F)

  SBE(Est_target, Est_reference, i)

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
}
