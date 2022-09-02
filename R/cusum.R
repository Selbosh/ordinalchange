#' Conditional expectation of latent process
#' 
#' @rdname condexp
#' 
#' @param Y An integer vector of observed categories.
#' @param mu Estimated latent mean vector, same length as \code{Y}.
#' @param cutpoints A vector of cutpoints \eqn{(c_0, c_1, c_2, \dots, c_{K})}
Z_tilde <- function(Y, mu, cutpoints) {
    stopifnot(length(Y) == length(mu))
    stopifnot(max(Y) < length(cutpoints))
    stopifnot(min(Y) >= 1)
    mu + (dnorm(cutpoints[Y] - mu) - dnorm(cutpoints[Y + 1])) /
         (pnorm(cutpoints[Y + 1] - mu) - pnorm(cutpoints[Y] - mu))
}

#' @rdname condexp
X_tilde <- function(Y, mu, cutpoints) {
    stopifnot(length(Y) == length(mu))
    stopifnot(max(Y) < length(cutpoints))
    stopifnot(min(Y) >= 1)
    (dnorm(cutpoints[Y] - mu) - dnorm(cutpoints[Y + 1])) /
         (pnorm(cutpoints[Y + 1] - mu) - pnorm(cutpoints[Y] - mu))
}

#' Bartlett estimator
#' @param X A vector of estimates for the latent process. For example \eqn{\tilde{X}}.
bartlett_eta <- function(X) {
    n <- length(X)
    q <- floor(n ^ (1 / 3))
    X_bar <- mean(X)
    eta2 <- var(X) / (n - 1) * n
    for (ell in seq_len(q)) {
        eta2 + 2 * (1 - ell / (q + 1)) * sum((head(X, -ell) - X_bar) * (tail(X, -ell) - X_bar)) / n
    }
    return(eta2)
}

#' Detect changepoints using CUSUM test statistics
#' 
#' This function calculates a CUSUM test statistics at various time points and returns the highest one found.
#' This value can then be compared with (empirical) percentiles of the distribution of the test statistic,
#' to decide whether or not a changepoint has been found.
#' 
#' @return 
#' A vector containing the position of the best changepoint (\eqn{\tau}) and the respective value of the test statistic (\eqn{M}).
#' 
#' @param Y Integer vector of responses.
#' @param L Period of the sinusoidal component of the latent process.
#' @param K Number of ordinal categories. \code{Y} is assumed to take on values \eqn{1, \dots, K}.
#' @param t Candidate positions for the changepoints. Defaults to \eqn{2, \dots, n - 1}.
#' 
#' @source Section 3 of \url{https://doi.org/10.1002/env.2752}
#' 
#' @rdname cusum
#' @export
changepoint_detection <- function(Y, L, K = max(Y), t = 2:(length(Y) - 1),
                                  verbose = FALSE) {
  n <- length(Y)
  maximum <- -Inf
  best_changepoint <- NA
  for (tau in t) {
    if (verbose)
        message('Fitting model for timepoint tau = ', tau)
    model <- seq_estimation(Y, tau = tau, K = K, L = L)
    theta <- model$theta
    rho <- model$rho$par
    mu <- ordinalchange:::mu(alpha0 = theta[K-1],
                             alpha1 = theta[K],
                             A = theta[K + 1],
                             D = theta[K + 2],
                             Delta = theta[K + 3],
                             Y, tau, L = L)
    X_tilde <- ordinalchange:::X_tilde(Y, mu, cutpoints = c(-Inf, 0, theta[1:(K-2)], Inf))
    eta <- sqrt(ordinalchange:::bartlett_eta(X_tilde))
    test_stat <- abs(sum(X_tilde[1:tau]) - tau / n * sum(X_tilde)) / sqrt(n)
    if (test_stat > maximum) {
        if (verbose)
            message('\tBetter changepoint found: tau = ', tau,
                    ' with M = ', test_stat)
        maximum <- test_stat
        best_changepoint <- tau
    }
  }
  c(tau = best_changepoint, M = maximum)
}