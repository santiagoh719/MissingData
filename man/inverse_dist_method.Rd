\name{IDW}
\alias{IDW}
\alias{IDW_simple}
\alias{RNNWM}

\title{Inverse distance weigth interpolation methods}
\description{
  Fill the ieth data with an inverse distance weigth method.
}
\usage{
  IDW(Est_reference, i, k, cord = T , distans = F)

  IDW_simple(Est_reference, i, cord = T , distans = F)

  RNNWM(Est_target, Est_reference, i, k, reg= 1)

}
\arguments{
  \item{Est_target}{numeric vector, data in which the ieth value will be interpolated}

  \item{Est_reference}{matrix or data.frame, each column is a different station}

  \item{i}{integer, index of the value to interpolated in target}

  \item{k}{numeric, exponent of the inverse distance weigth}

  \item{cord}{matrix or data.frame, first column Longitude second column Latitude, first row coordinates of the target, then in order the coordinates of Est_reference}

  \item{distans}{numeric vector,  in order the distans from target to Est_reference}

  \item{reg}{integer value, large of the paterns that will be look for, default 1 acept 1 or 2}
}

\value{
  \code{IDW}, returns the ieth value of target with an interpolation using a inverse distance weigth model
  \code{IDW_simple}, returns the ieth value of target with an interpolation using a inverse distance weigth model with k=1
  \code{RNNWM}, returns the ieth value of target with an interpolation using a inverse distance weigth model, with a different mesearu of distance.
}
\details{
    \code{IDW} compute the ieth value of target with an interpolation using a inverse distance weigth model. \code{IDW_simple} is as \code{IDW} with k=1.
  \code{RNNWM}, compute the ieth value of target with an interpolation using a inverse distance weigth model but with adifferent mesearu of distance.It use a surrogate measure for distance base on Teegavarapu and Chandramouli 2005.
}

