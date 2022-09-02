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

#' CUSUM test statistics
#' 
#' @source Section 3 of \url{https://doi.org/10.1002/env.2752}
#' 
#' @rdname cusum
M_X <- function(Y, mu, cutpoints) {
    X_tilde <- X_tilde(Y, mu, cutpoints)
    eta <- sqrt(bartlett_eta(X_tilde))
    n <- length(Y)
    maximum <- -Inf
    best_tau <- 0
    for (tau in 1:n) {
        test <- abs(sum(X_tilde - tau / n * sum(X_tilde)) / sqrt(n))
        if (test > maximum) {
            maximum <- test
            best_tau <- tau
        }
    }
    c(M = maximum / eta, tau = best_tau)
}

#' Detect changepoints using CUSUM test statistics
#' 
#' 
