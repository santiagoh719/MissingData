\name{TSLR}
\alias{TSLR_ranked}
\alias{TSLR}
\title{Theil-Sen missing value filled}
\description{
  Filled the ieth data with a Theil-Sen regression.
}
\usage{
  TSLR(Est_target, Est_reference,i, criterion='Max_Cor',cord=T,distans=F)

  TSLR_ranked(Est_target, Est_reference,i, criterion='Max_Cor',cord=T,distans=F)

}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{criterion}{criterion is which criteria to be used to select reference estation, Max_Cor (default) for the station with maximal correlation, or Nearest for the nearest station. If Nearest either cord or distans must be given}

  \item{cord}{matrix or data.frame, first column Longitude second column Latitude, first row coordinates of the target, then in order the coordinates of Est_reference}

  \item{distans}{numeric vector,  in order the distans from target to Est_reference}

}

\value{
  \code{TSLR},returns ieth data with a Theil-Sen regression of the data.

  \code{TSLR_ranked} returns the ieth data with a Theil-Sen regression of the ranked data.

}

