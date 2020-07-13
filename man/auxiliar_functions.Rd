\name{patern_distans_change}
\alias{patern_distans_change}
\title{Patern identify and distanse base on paterns}
\description{
  \code{patern_search} look for patterns inside a vector. \code{patern_distans_change} compute a surrogate measure for distance base on Teegavarapu and Chandramouli 2005.
}
\usage{
  patern_distans_change(Est_target, Est_reference, reg= 1)
}
\arguments{
  \item{Est_target}{numeric vector, station data for which the distans will be calculated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{reg}{integer value, large of the paterns that will be look for, default 1 acept 1 or 2}

}



\value{
  \code{patern_distans_change}, returns a vector with a surrogate measure for distance base of the increased/decreased paterns.
    \code{patern_search}, returns a vector with the indexes of the patern in the vector.

}
\details{
\code{patern_distans_change} compute a surrogate measure for distance base on Teegavarapu and Chandramouli 2005.
\code{patern_search} look for patterns inside a vector.
}

