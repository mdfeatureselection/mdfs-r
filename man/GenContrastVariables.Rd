% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{GenContrastVariables}
\alias{GenContrastVariables}
\title{Generate contrast variables from data}
\usage{
GenContrastVariables(data, n.contrast = max(ncol(data), 30))
}
\arguments{
\item{data}{data organized in matrix with separate variables in columns}

\item{n.contrast}{number of constrast variables (defaults to max of 1/10 of variables number and 30)}
}
\value{
A list with the following key names:
 \itemize{
  \item \code{indices} -- vector of indices of input variables used to construct contrast variables
  \item \code{x} -- data with constrast variables appended to it
  \item \code{mask} -- vector of booleans making it easy to select just contrast variables
 }
}
\description{
Generate contrast variables from data
}
\examples{
GenContrastVariables(madelon$data)
}
