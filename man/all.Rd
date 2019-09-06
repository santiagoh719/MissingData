\name{interpolate_all}
\alias{interpolate_all}
\title{Interpolate all NAs from a serie}
\description{
  \code{interpolate_all} filled all NAs from a serie with a given function.
}
\usage{
  interpolate_all(Est_target, fun,...)

}
\arguments{
  \item{Est_target}{numeric vector, data in which the NA values will be interpolated}

  \item{fun}{function to use for interpolate the data, any function of this package can be used}

  \item{...}{other arguments used by fun}

}

\value{
  \code{interpolate_all}, returns the Est_target vector without NAs values or with some NAs values that could not be interpolated with the given method.
}
\details{
  \code{interpolate_all} filled all NAs from a serie with a given function. It can hapend that not all NAs values can be interpolated with the given function, in that case it will be printed that, and the return will be Est_target with the interpolation that could be done.
}

