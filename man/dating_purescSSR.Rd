% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{dating_purescSSR}
\alias{dating_purescSSR}
\title{Compute global break dates for pure structural change model}
\usage{
dating_purescSSR(y, z, m, h)
}
\arguments{
\item{y}{A (\code{T x 1}) vector with endogenous variable.}

\item{z}{A (\code{T x q}) matrix with explanatory variables subject to change.}

\item{m}{An integer determining the number of breaks to find.}

\item{h}{An integer determining the minimum length of a regime.}
}
\description{
This is the main procedure which calculates the break points that globally
 minimizes the SSR. It returns optimal dates and associated SSR for all numbers of breaks less than or equal to m.
}
\details{
Note: This code is an adaptation of the one originally written by Yohei 
Yamamoto and Pierre Perron for MATLAB. Original code files can be found on 
Pierre Perron's website: https://blogs.bu.edu/perron/codes/
}
\references{
Bai, Jushan & Pierre Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes," \emph{Econometrica}, vol 66, 47-78.

Bai, Jushan & Pierre Perron (2003), "Computation and Analysis of Multiple Structural Change Models," \emph{Journal of Applied Econometrics}, 18, 1-22.
}
