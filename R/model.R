#' Estimate the AOP(1) model from Li and Lu (2022)
#' 
#' The parameter vector is \eqn{\theta = (c_2, \dots, c_{K-1}, \alpha_0, \alpha_1, A, D, \Delta)}.
#' Starting values for \eqn{\alpha_1, A, D, \Delta} are set to 0.
#' Starting values for \eqn{\alpha_0} and thresholds \eqn{c_2,\dots,c_{K-1}} are \eqn{\alpha_0 = -\Phi^s that \eqn{Y_t} lies in category \eqn{j}.
#' {-1}(p_1)} and \eqn{c_k=-\Phi^{-1}(\sum_{j=1}^k p_j) + \alpha_0}, where \eqn{p_j} is the sample proportion of time
#' @param Y Integer vector.
#' @param K Integer number of response categories, denoted \eqn{k = 1, \dots, K}.
#' @param iterations Number of iterations for Newton's method.
#' 
#' @source \url{https://doi.org/10.1002/env.2752}
#' 
#' @export
model_estimation <- function(Y, K = max(Y), tau = 10, T = 7, ...) {
    # Input checking.
    stopifnot(all(Y > 0))
    stopifnot(K > 1)
    stopifnot(iterations > 0)
    
    # Initialization.
    cutpoint <- vector(length = K - 2)
    alpha0 <- -qnorm(mean(Y == 1))
    for (k in 2:(K - 1))
      cutpoint[k - 1] <- qnorm(mean(Y <= k)) - alpha0
    alpha1 <- A <- D <- Delta <- 0
    init_theta <- c(cutpoint, alpha0, alpha1, A, D, Delta)


    # Newton's method.
    result <- nlm(Q, init_theta, Y = Y, tau = tau, T = T)
}

mu <- function(alpha0, alpha1, A, D, Delta, Y, tau, T) {
    t <- seq_along(Y)
    n <- length(Y)
    s <- A * cos(2 * pi * t / T) + D * sin(2 * pi * t / T)
    alpha0 + alpha1 * t / n + s + Delta * (t >= tau)
}

#' Sum of squared discrepancies between observations and their expected values
#' 
#' @param theta A vector of parameters.
#' @param Y A vector of ordered categorical responses.
#' @param tau A candidate changepoint.
#' @param T The known period of the seasonal component.
#' 
#' @export
Q <- function(theta, Y, tau, T) {
    # Unpack the parameters
    K <- length(theta) - 5 + 2
    cutpoints <- c(0, head(theta, K - 2), +Inf)
    alpha0 <- theta[K - 1]
    alpha1 <- theta[K]
    A <- theta[K + 1]
    D <- theta[K + 2]
    Delta <- theta[K + 3]

    # Compute level probabilities and expected values.
    n <- length(Y)
    mu <- mu(alpha0, alpha1, A, D, Delta, Y, tau, T)
    E_Y <- xi <- omega <- matrix(nrow = n, ncol = K)
    E_Y[, 1] <- pnorm(0 - mu) - 0 # c0 = -Inf
    xi[, 1] <- dnorm(0 - mu) - 0  # c0 = -Inf
    omega[, 1] <- (Y == 1) - E_Y[, 1]
    for (k in 2:K) {
      E_Y[, k] <- pnorm(cutpoints[k] - mu) - pnorm(cutpoints[k - 1] - mu)
      xi[, k] <- dnorm(cutpoints[k] - mu) - dnorm(cutpoints[k - 1] - mu)
      omega[, k] <- (Y == k) - E_Y[, k]
    }

    # Compute the sum of squared differences.
    result <- sum((outer(Y, 1:K, '==') - E_Y)^2)
    # Compute the gradient.
    attr(result, 'gradient') <- gradient_theta(theta, omega, xi, mu)
    return(result)
}

gradient_theta <- function(theta, omega, xi, mu) {
    # Unpack the parameters
    K <- length(theta) - 5 + 2
    cutpoints <- c(0, head(theta, K - 2))
    alpha0 <- theta[K - 1]
    alpha0 <- theta[K]
    A <- theta[K + 1]
    D <- theta[K + 2]
    Delta <- theta[K + 3]
    gradient(cutpoints, alpha0, alpha1, A, D, Delta, omega, xi, mu)
}

gradient <- function(cutpoints, alpha0, alpha1, A, D, Delta, omega, xi, mu) {
    d_cutpoints <- vector(length = length(cutpoints)) # d_c1 is redundant
    for (k in 2:(K - 1)) {
        d_cutpoints[k] <- 2 * sum((omega[, k + 1] - omega[, k]) * dnorm(cutpoints[k] - mu))
    }

    t <- seq_along(Y)
    n <- length(Y)
    d_alpha0 <- 2 * sum(omega * xi * 1)
    d_alpha1 <- 2 * sum(omega * xi * t / n)
    d_A <- 2 * sum(omega * xi * cos(2 * pi / T * t))
    d_D <- 2 * sum(omega * xi * sin(2 * pi / T * t))
    d_Delta <- 2 * sum (omega * xi * (t >= tau))

    c(d_cutpoints[-1], d_alpha0, d_alpha1, d_A, d_D, d_Delta)
}

# hessian <- function(theta) {
# }