\name{Optimal Interpolation}
\alias{OI}
\title{Optimal Interpolation}
\description{
  Filled the ieth data with Optimal Interpolation.
}
\usage{
  OI(Est_target, Est_reference, i, cord= T, criterion = 'Cor' )

}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{criterion}{criterion is which criteria to be used to select the fisrt guest value, Max_Cor (default) for the station with maximal correlation, or Nearest for the nearest station. If Nearest cord must be given}

  \item{cord}{matrix or data.frame, first column Longitude second column Latitude, first row coordinates of the target, then in order the coordinates of Est_reference}

}

\value{
  \code{OI},returns ieth data interpolate with the Optimal Interpolation method discribe in E ISCHEID 1999.

}

