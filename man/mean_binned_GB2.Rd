% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ineqQuantile.R
\name{mean_binned_GB2}
\alias{mean_binned_GB2}
\title{Means of the ranges according to a GB2 distribution}
\usage{
mean_binned_GB2(p, q, mu = 1, sigma = 1, nu = 1, tau = 0.5)
}
\arguments{
\item{p}{numeric, vector of probabilities}

\item{q}{numeric, vector of quantiles}

\item{mu}{numeric, vector of location parameter values}

\item{sigma}{numeric, vector of scale parameter values}

\item{nu}{numeric, vector of skewness parameter values}

\item{tau}{numeric, vector of kurtosis parameter values}
}
\description{
Means of the ranges according to a GB2 distribution
}
\references{
Belz (2019), \emph{Estimating Inequality Measures from Quantile Data} \url{https://halshs.archives-ouvertes.fr/halshs-02320110}
}
\seealso{
\code{\link{compute_LC}}
}
