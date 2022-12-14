% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cusum.R
\name{Z_tilde}
\alias{Z_tilde}
\alias{X_tilde}
\title{Conditional expectation of latent process}
\usage{
Z_tilde(Y, mu, cutpoints)

X_tilde(Y, mu, cutpoints)
}
\arguments{
\item{Y}{An integer vector of observed categories.}

\item{mu}{Estimated latent mean vector, same length as \code{Y}.}

\item{cutpoints}{A vector of cutpoints \eqn{(c_0, c_1, c_2, \dots, c_{K})}}
}
\description{
Conditional expectation of latent process
}
