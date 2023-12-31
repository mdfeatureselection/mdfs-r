% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tuples.R
\name{ComputeInterestingTuples}
\alias{ComputeInterestingTuples}
\title{Interesting tuples}
\usage{
ComputeInterestingTuples(
  data,
  decision = NULL,
  dimensions = 2,
  divisions = 1,
  discretizations = 1,
  seed = NULL,
  range = NULL,
  pc.xi = 0.25,
  ig.thr = 0,
  I.lower = NULL,
  interesting.vars = vector(mode = "integer"),
  require.all.vars = FALSE,
  return.matrix = FALSE,
  stat_mode = "MI",
  average = FALSE
)
}
\arguments{
\item{data}{input data where columns are variables and rows are observations (all numeric)}

\item{decision}{decision variable as a binary sequence of length equal to number of observations}

\item{dimensions}{number of dimensions (a positive integer; 5 max)}

\item{divisions}{number of divisions (from 1 to 15)}

\item{discretizations}{number of discretizations}

\item{seed}{seed for PRNG used during discretizations (\code{NULL} for random)}

\item{range}{discretization range (from 0.0 to 1.0; \code{NULL} selects probable optimal number)}

\item{pc.xi}{parameter xi used to compute pseudocounts (the default is recommended not to be changed)}

\item{ig.thr}{IG threshold above which the tuple is interesting (0 and negative mean no filtering)}

\item{I.lower}{IG values computed for lower dimension (1D for 2D, etc.)}

\item{interesting.vars}{variables for which to check the IGs (none = all)}

\item{require.all.vars}{boolean whether to require tuple to consist of only interesting.vars}

\item{return.matrix}{boolean whether to return a matrix instead of a list (ignored if not using the optimised method variant)}

\item{stat_mode}{character, one of: "MI" (mutual information, the default; becomes information gain when \code{decision} is given), "H" (entropy; becomes conditional entropy when \code{decision} is given), "VI" (variation of information; becomes target information difference when \code{decision} is given); decides on the value computed}

\item{average}{boolean whether to average over discretisations instead of maximising (the default)}
}
\value{
A \code{\link{data.frame}} or \code{\link{NULL}} (following a warning) if no tuples are found.

 The following columns are present in the \code{\link{data.frame}}:
 \itemize{
   \item \code{Var} -- interesting variable index
   \item \code{Tuple.1, Tuple.2, ...} -- corresponding tuple (up to \code{dimensions} columns)
   \item \code{IG} -- information gain achieved by \code{var} in \code{Tuple.*}
 }

 Additionally attribute named \code{run.params} with run parameters is set on the result.
}
\description{
Interesting tuples
}
\details{
If running in 2D and no filtering is applied, this function is able to run in an
optimised fashion. It is recommended to avoid filtering in 2D if only it is
feasible.

This function calculates what \code{stat_mode} dictates.
When \code{decision} is omitted, the \code{stat_mode} is calculated on the descriptive variables.
When \code{decision} is given, the \code{stat_mode} is calculated on the decision variable, conditional on the other variables.
Translate "IG" to that value in the rest of this function's description.
}
\examples{
\donttest{
ig.1d <- ComputeMaxInfoGains(madelon$data, madelon$decision, dimensions = 1, divisions = 1,
                             range = 0, seed = 0)
ComputeInterestingTuples(madelon$data, madelon$decision, dimensions = 2, divisions = 1,
                         range = 0, seed = 0, ig.thr = 100, I.lower = ig.1d$IG)
}
}
